!--------------------------------------------------------------------------------------------------!
!  FORTNET: A Behler-Parrinello-Neural-Network Implementation                                      !
!  Copyright (C) 2020 - 2021  T. W. van der Heide                                                  !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Implements Atom-centered Symmetry Functions to generate input features.
module fnet_acsf

  use h5lt
  use hdf5

  use dftbp_accuracy, only : dp
  use dftbp_message, only : error
  use dftbp_constants, only: pi
  use dftbp_charmanip, only : i2c, tolower
  use dftbp_typegeometry, only : TGeometry
  use dftbp_dynneighlist, only : TDynNeighList, TDynNeighList_init, TNeighIterator,&
      & TNeighIterator_init

  use fnet_nestedtypes, only : TIntArray1D, TRealArray2D, TRealArray4D, TEnv
  use fnet_intmanip, only : getUniqueInt, getNumberOfUniqueInt
  use fnet_hdf5fx, only : h5ltfx_read_dataset_double_f, h5ltfxmake_dataset_double_f

#:if WITH_MPI
  use fnet_mpifx
  use fnet_parallel, only : getStartAndEndIndex
#:endif

  implicit none

  private

  public :: TAcsf, TAcsf_init, TMultiAcsfVals, TMultiAcsfPrimeVals
  public :: TGFunction, TGFunction_init, TGFunctions

  type :: TGFunction

    !> type of function
    character(len=2) :: type

    !> cutoff radius, defining the sphere to search for neighboring atoms
    real(dp) :: rCut

    !> kappa scaling
    real(dp) :: kappa

    !> center of Gaussian
    real(dp) :: rs

    !> width of Gaussian
    real(dp) :: eta

    !> lambda parameter of angular part
    real(dp) :: lambda

    !> parameter of angular resolution
    real(dp) :: xi

    !> dataset index for atomIDs
    integer :: atomId

    !> atomic number for the species-resolved scheme
    integer :: atomicNumbers(2)

    !> true, if the instance is a radial function (G1, G2, G3)
    logical :: tRadial

    !> true, if the instance is an angular function (G4, G5)
    logical :: tAngular

  end type TGFunction


  type :: TGFunctions

    !> wrapper around multiple G-functions
    type(TGFunction), allocatable :: func(:)

  contains

    procedure :: append => TGFunctions_append
    procedure :: fromAutoScheme => TGFunctions_fromAutoScheme

  end type TGFunctions


  type :: TMultiAcsfVals

    !> array of ACSF value instances
    type(TRealArray2D), allocatable :: vals(:)

  end type TMultiAcsfVals


  type :: TMultiAcsfPrimeVals

    !> array of ACSF derivative value instances
    type(TRealArray4D), allocatable :: vals(:)

  end type TMultiAcsfPrimeVals


  type :: TAcsf

    !> functions emerging from an automatic generation scheme
    type(TGFunctions) :: gFunctions

    !> wrapper around multiple ACSF value instances
    type(TMultiAcsfVals) :: vals

    !> wrapper around multiple ACSF derivative value instances
    type(TMultiAcsfPrimeVals) :: valsPrime

    !> storage container of means and variances to calculate z-score
    real(dp), allocatable :: zPrec(:,:)

    !> true, if z-score standardization should be applied
    logical :: tZscore

  contains

    procedure :: getMeansAndVariances => TAcsf_getMeansAndVariances
    procedure :: applyZscore => TAcsf_applyZscore
    procedure :: applyZscorePrime => TAcsf_applyZscorePrime
    procedure :: calculate => TAcsf_calculate
    procedure :: calculatePrime => TAcsf_calculatePrime
    procedure :: fromFile => TAcsf_fromFile
    procedure :: toFile => TAcsf_toFile

  #:if WITH_MPI
    procedure :: syncConfig => TAcsf_syncConfig
  #:endif

  end type TAcsf

  !> Chunk size to use when obtaining neighbours dynamically via an iterator
  !> (used to artificially restrict memory usage to a certain amount)
  integer, parameter :: iterChunkSize = 1000


contains

  !> Initialises an ACSF instance.
  pure subroutine TAcsf_init(this, functions, tZscore)

    !> representation of ACSF mappings
    type(TAcsf), intent(out) :: this

    !> wrapper around multiple G-functions
    type(TGFunctions), intent(in) :: functions

    !> true, if z-score standardization should be applied
    logical, intent(in), optional :: tZscore

    if (present(tZscore)) then
      this%tZscore = tZscore
    else
      this%tZscore = .false.
    end if

    this%gFunctions = functions

  end subroutine TAcsf_init


  !> Initialises a single ACSF G-function (G1, G2, G3, G4 or G5).
  subroutine TGFunction_init(this, type, rCut, kappa, rs, eta, lambda, xi, atomId)

    !> representation of an ACSF function
    type(TGFunction), intent(out) :: this

    !> type of G-function (G1, G2, G3, G4, G5)
    character(len=2), intent(in) :: type

    !> cutoff radius
    real(dp), intent(in) :: rCut

    !> kappa scaling
    real(dp), intent(in), optional :: kappa

    !> center of Gaussian
    real(dp), intent(in), optional :: rs

    !> width of Gaussian
    real(dp), intent(in), optional :: eta

    !> lambda parameter of angular part
    real(dp), intent(in), optional :: lambda

    !> parameter of angular resolution
    real(dp), intent(in), optional :: xi

    !> optional atom identifier index
    integer, intent(in), optional :: atomId

    !> error messages
    character(len=:), allocatable :: msg1, msg2

    msg1 = 'Superflous argument(s) for selected function type provided.'
    msg2 = 'Missing argument(s) for selected function type.'

    this%type = tolower(type)
    this%rCut = rCut

    if (present(atomId)) then
      this%atomId = atomId
    else
      this%atomId = 0
    end if

    ! provide dummy values for now
    this%kappa = 0.0_dp
    this%rs = 0.0_dp
    this%eta = 0.0_dp
    this%lambda = 0.0_dp
    this%xi = 0.0_dp
    this%atomicNumbers(:) = [0, 0]
    this%tRadial = .false.
    this%tAngular = .false.

    select case (this%type)

    case ('g1')
      if (present(kappa) .or. present(rs) .or. present(eta) .or. present(lambda) .or. present(xi))&
          & then
        call error(msg1)
      end if
      this%tRadial = .true.
      this%tAngular = .false.
    case ('g2')
      if (.not. (present(rs) .and. present(eta))) then
        call error(msg2)
      elseif (present(kappa) .or. present(lambda) .or. present(xi)) then
        call error(msg1)
      elseif (present(rs) .and. present(eta)) then
        this%rs = rs
        this%eta = eta
      end if
      this%tRadial = .true.
      this%tAngular = .false.
    case ('g3')
      if (.not. present(kappa)) then
        call error(msg2)
      elseif (present(rs) .or. present(eta) .or. present(lambda) .or. present(xi)) then
        call error(msg1)
      elseif (present(kappa)) then
        this%kappa = kappa
      end if
      this%tRadial = .true.
      this%tAngular = .false.
    case ('g4', 'g5')
      if (.not. (present(eta) .or. present(lambda) .or. present(xi))) then
        call error(msg2)
      elseif (present(kappa) .or. present(rs)) then
        call error(msg1)
      elseif (present(eta) .or. present(lambda) .or. present(xi)) then
        this%eta = eta
        this%lambda = lambda
        this%xi = xi
      end if
      this%tRadial = .false.
      this%tAngular = .true.
    case default
      call error('Invalid function type, aborting initialization.')

    end select

  end subroutine TGFunction_init


  !> Initialises multiple G-functions according to an automatic parameter generation scheme.
  subroutine TGFunctions_fromAutoScheme(this, rCut, nRadial, nAngular, atomId)

    !> representation of ACSF functions
    class(TGFunctions), intent(out) :: this

    !> cutoff radius
    real(dp), intent(in) :: rCut

    !> number of radial symmetry functions
    integer, intent(in) :: nRadial

    !> number of angular symmetry functions
    integer, intent(in) :: nAngular

    !> optional atom identifier index
    integer, intent(in), optional :: atomId

    !> stepsize for center parameter
    real(dp) :: rsStep

    !> width parameter for G2 and G5 function
    real(dp) :: g2eta, g5eta

    !> peak positions for G2 function
    real(dp), allocatable :: g2rs(:)

    !> lambda parameters for G5 function
    real(dp), allocatable :: g5lambda(:)

    !> xi parameters for G5 function
    real(dp), allocatable :: g5xi(:)

    !> auxiliary variable
    real(dp) :: xi

    !> auxiliary variables
    integer :: ii, jj, ind, iRadial, iAngular, identifier

    if (present(atomId)) then
      identifier = atomId
    else
      identifier = 0
    end if

    rsStep = rCut / (nRadial - 1)
    g2eta = 5.0_dp * log(10.0_dp) / (2.0_dp * rsStep)**2

    allocate(g2rs(nRadial))
    g2rs(1) = 0.0_dp

    if (nRadial > 1) then
      g2rs(2:) = [(ii * rsStep, ii = 1, nRadial - 1)]
    end if

    g5eta = 2.0_dp * log(10.0_dp) / (rCut)**2

    allocate(g5lambda(nAngular))
    allocate(g5xi(nAngular))

    ind = 1

    do ii = 0, ceiling(nAngular / 2.0_dp - 1.0_dp)
      if (nAngular <= 2) then
        xi = 1.0_dp
      else
        xi = 1.0_dp + ii * 30.0_dp / (nAngular - 2.0_dp)
      end if
      do jj = 1, -1, -2
        g5lambda(ind) = jj
        g5xi(ind) = xi
        if (ind >= nAngular) exit
        ind = ind + 1
      end do
    end do

    allocate(this%func(nRadial + nAngular))

    do iRadial = 1, nRadial
      call TGFunction_init(this%func(iRadial), 'G2', rCut, eta=g2eta, rs=g2rs(iRadial),&
          & atomId=identifier)
    end do

    do iAngular = 1, nAngular
      call TGFunction_init(this%func(iAngular + nRadial), 'G5', rCut, xi=g5xi(iAngular), eta=g5eta,&
          & lambda=g5lambda(iAngular), atomId=identifier)
    end do

  end subroutine TGFunctions_fromAutoScheme


  !> Appends G-functions to the current instance.
  pure subroutine TGFunctions_append(this, functions)

    !> representation of ACSF functions
    class(TGFunctions), intent(inout) :: this

    !> representation of ACSF functions to append
    type(TGFunctions), intent(in) :: functions

    !> temporary storage
    type(TGFunctions) :: tmp

    if (.not. allocated(functions%func)) then
      return
    end if

    if (allocated(this%func)) then
      call move_alloc(this%func, tmp%func)
      allocate(this%func(size(tmp%func) + size(functions%func)))
      this%func(1:size(tmp%func)) = tmp%func
      this%func(size(tmp%func)+1:) = functions%func
    else
      this%func = functions%func
    end if

  end subroutine TGFunctions_append


  !> Initialises a wrapper around multiple ACSF value instances.
  pure subroutine TMultiAcsfVals_init(this, geos, nAcsfVals)

    !> wrapper instance around multiple ACSF value types
    type(TMultiAcsfVals), intent(out) :: this

    !> system geometry container
    type(TGeometry), intent(in) :: geos(:)

    !> total number of ACSF mappings per atom
    integer, intent(in) :: nAcsfVals

    !> auxiliary variable
    integer :: iGeo

    allocate(this%vals(size(geos)))

    do iGeo = 1, size(geos)
      allocate(this%vals(iGeo)%array(nAcsfVals, geos(iGeo)%nAtom))
      this%vals(iGeo)%array(:,:) = 0.0_dp
    end do

  end subroutine TMultiAcsfVals_init


  !> Initialises a wrapper around multiple ACSF derivative value instances.
  pure subroutine TMultiAcsfPrimeVals_init(this, geos, nAcsfVals)

    !> wrapper instance around multiple ACSF derivative value types
    type(TMultiAcsfPrimeVals), intent(out) :: this

    !> system geometry container
    type(TGeometry), intent(in) :: geos(:)

    !> total number of ACSF mappings per atom
    integer, intent(in) :: nAcsfVals

    !> auxiliary variable
    integer :: iSys

    allocate(this%vals(size(geos)))

    do iSys = 1, size(geos)
      allocate(this%vals(iSys)%array(3, nAcsfVals, geos(iSys)%nAtom, geos(iSys)%nAtom))
      this%vals(iSys)%array(:,:,:,:) = 0.0_dp
    end do

  end subroutine TMultiAcsfPrimeVals_Init


  !> Calculates the means and variances for z-score preconditioning.
  subroutine TAcsf_getMeansAndVariances(this, weights)

    !> representation of ACSF mappings
    class(TAcsf), intent(inout) :: this

    !> weighting of each corresponding datapoint
    integer, intent(in) :: weights(:)

    !> auxiliary variables
    integer :: iGeo, iAtom, nTotAtoms

    if (allocated(this%zPrec)) then
      call error('Container for acsf means and variances is already allocated.')
    else
      allocate(this%zPrec(size(this%gFunctions%func), 2))
      this%zPrec(:,:) = 0.0_dp
    end if

    nTotAtoms = 0

    ! calculate means of the different inputs
    do iGeo = 1, size(this%vals%vals)
      do iAtom = 1, size(this%vals%vals(iGeo)%array, dim=2)
        nTotAtoms = nTotAtoms + weights(iGeo)
        this%zPrec(:, 1) = this%zPrec(:, 1) + real(weights(iGeo), dp)&
            & * this%vals%vals(iGeo)%array(:, iAtom)
      end do
    end do

    this%zPrec(:, 1) = this%zPrec(:, 1) / real(nTotAtoms, dp)

    ! calculate variances of the different inputs
    do iGeo = 1, size(this%vals%vals)
      do iAtom = 1, size(this%vals%vals(iGeo)%array, dim=2)
        this%zPrec(:, 2) = this%zPrec(:, 2) + real(weights(iGeo), dp)&
            & * (this%vals%vals(iGeo)%array(:, iAtom) - this%zPrec(:, 1))**2
      end do
    end do

    this%zPrec(:, 2) = sqrt(this%zPrec(:, 2) / real(nTotAtoms, dp))

  end subroutine TAcsf_getMeansAndVariances


  !> Applies a z-score preconditioning to the ACSF features.
  subroutine TAcsf_applyZscore(this)

    !> representation of ACSF mappings
    class(TAcsf), intent(inout) :: this

    !> auxiliary variables
    integer :: iGeo, iAcsf

    if (.not. allocated(this%zPrec)) then
      call error('Cannot apply z-score standardization without means and variances.')
    end if

    ! apply z-score standardization to ACSF mappings
    do iGeo = 1, size(this%vals%vals)
      do iAcsf = 1, size(this%vals%vals(iGeo)%array, dim=1)
        if (this%zPrec(iAcsf, 2) < 1e-08_dp) cycle
        this%vals%vals(iGeo)%array(iAcsf, :) = (this%vals%vals(iGeo)%array(iAcsf, :) -&
            & this%zPrec(iAcsf, 1)) / this%zPrec(iAcsf, 2)
      end do
    end do

  end subroutine TAcsf_applyZscore


  !> Account for possible z-score preconditioning in the ACSF derivatives.
  subroutine TAcsf_applyZscorePrime(this)

    !> representation of ACSF mappings
    class(TAcsf), intent(inout) :: this

    !> auxiliary variables
    integer :: iGeo, iAcsf

    if (.not. allocated(this%zPrec)) then
      call error('Cannot apply z-score standardization without means and variances.')
    end if

    ! apply z-score standardization to ACSF derivatives
    do iGeo = 1, size(this%valsPrime%vals)
      do iAcsf = 1, size(this%valsPrime%vals(iGeo)%array, dim=2)
        if (this%zPrec(iAcsf, 2) < 1e-08_dp) cycle
        this%valsPrime%vals(iGeo)%array(:, iAcsf, :, :)&
            & = this%valsPrime%vals(iGeo)%array(:, iAcsf, :, :) / this%zPrec(iAcsf, 2)
      end do
    end do

  end subroutine TAcsf_applyZscorePrime


  !> Calculates ACSF mappings for given geometries.
  subroutine TAcsf_calculate(this, geos, env, localAtToAtNum, extFeaturesInp, weights, zPrec)

    !> representation of ACSF mappings
    class(TAcsf), intent(inout) :: this

    !> system geometry container
    type(TGeometry), intent(in) :: geos(:)

    !> if compiled with mpi enabled, contains mpi communicator
    type(TEnv), intent(in) :: env

    !> index mapping local atom --> atomic number
    type(TIntArray1D), intent(in) :: localAtToAtNum(:)

    !> optional atom dependent scaling parameters for cutoff function
    type(TRealArray2D), intent(in), optional :: extFeaturesInp(:)

    !> optional weighting of each corresponding datapoint
    integer, intent(in), optional :: weights(:)

    !> storage container of means and variances to calculate z-score
    real(dp), intent(in), optional :: zPrec(:,:)

    !> atom dependent scaling parameters for cutoff function
    type(TRealArray2D), allocatable :: extFeatures(:)

    !> temporary storage of ACSF values of each node
    type(TMultiAcsfVals) :: tmpVals

    !> weighting of each corresponding datapoint
    integer, allocatable :: weighting(:)

    !> true, if current process is the lead
    logical :: tLead

    !> auxiliary variables
    integer :: iSys, iStart, iEnd

    allocate(weighting(size(geos)))

    if (present(weights)) then
      weighting(:) = weights
    else
      weighting(:) = 1
    end if

    if (present(zPrec)) then
      this%zPrec = zPrec
    end if

    ! prevent for accessing an unallocated array
    if (present(extFeaturesInp)) then
      extFeatures = extFeaturesInp
    else
      allocate(extFeatures(size(geos)))
    end if

  #:if WITH_MPI
    tLead = env%globalMpiComm%lead
    call getStartAndEndIndex(size(geos), env%globalMpiComm%size, env%globalMpiComm%rank, iStart,&
        & iEnd)
  #:else
    tLead = .true.
    iStart = 1
    iEnd = size(geos)
  #:endif

    call TMultiAcsfVals_init(tmpVals, geos, size(this%gFunctions%func))
    call TMultiAcsfVals_init(this%vals, geos, size(this%gFunctions%func))

    lpSystem: do iSys = iStart, iEnd
      if (.not. allocated(extFeatures(iSys)%array)) allocate(extFeatures(iSys)%array(0,0))
      call iGeoAcsf(tmpVals%vals(iSys), geos(iSys), this%gFunctions%func,&
          & localAtToAtNum(iSys)%array, extFeatures(iSys)%array)
    end do lpSystem

  #:if WITH_MPI
    ! sync ACSF mappings between MPI nodes
    do iSys = 1, size(geos)
      call mpifx_allreduce(env%globalMpiComm, tmpVals%vals(iSys)%array, this%vals%vals(iSys)%array,&
          & MPI_SUM)
    end do
  #:else
    this%vals = tmpVals
  #:endif

    if (tLead .and. this%tZscore) then
      if (.not. allocated(this%zPrec)) then
        call this%getMeansAndVariances(weighting)
      end if
      call this%applyZscore()
    end if

  #:if WITH_MPI
    do iSys = 1, size(this%vals%vals)
      call mpifx_bcast(env%globalMpiComm, this%vals%vals(iSys)%array)
    end do
  #:endif

  end subroutine TAcsf_calculate


  !> Calculates ACSF derivatives for given geometries.
  subroutine TAcsf_calculatePrime(this, geos, env, localAtToAtNum, extFeaturesInp)

    !> representation of ACSF mappings
    class(TAcsf), intent(inout) :: this

    !> system geometry container
    type(TGeometry), intent(in) :: geos(:)

    !> if compiled with mpi enabled, contains mpi communicator
    type(TEnv), intent(in) :: env

    !> index mapping local atom --> atomic number
    type(TIntArray1D), intent(in) :: localAtToAtNum(:)

    !> optional atom dependent scaling parameters for cutoff function
    type(TRealArray2D), intent(in), optional :: extFeaturesInp(:)

    !> atom dependent scaling parameters for cutoff function
    type(TRealArray2D), allocatable :: extFeatures(:)

    !> temporary storage of ACSF derivatives of each node
    type(TMultiAcsfPrimeVals) :: tmpValsPrime

    !> true, if current process is the lead
    logical :: tLead

    !> auxiliary variables
    integer :: iSys, iStart, iEnd

    ! prevent for accessing an unallocated array
    if (present(extFeaturesInp)) then
      extFeatures = extFeaturesInp
    else
      allocate(extFeatures(size(geos)))
    end if

  #:if WITH_MPI
    tLead = env%globalMpiComm%lead
    call getStartAndEndIndex(size(geos), env%globalMpiComm%size, env%globalMpiComm%rank, iStart,&
        & iEnd)
  #:else
    tLead = .true.
    iStart = 1
    iEnd = size(geos)
  #:endif

    call TMultiAcsfPrimeVals_init(tmpValsPrime, geos, size(this%gFunctions%func))
    call TMultiAcsfPrimeVals_init(this%valsPrime, geos, size(this%gFunctions%func))

    lpSystem: do iSys = iStart, iEnd
      call iGeoAcsfGrad(tmpValsPrime%vals(iSys), geos(iSys), this%gFunctions%func,&
          & localAtToAtNum(iSys)%array, extFeatures=extFeatures(iSys)%array)
    end do lpSystem

  #:if WITH_MPI
    ! sync ACSF mappings between MPI nodes
    do iSys = 1, size(geos)
      call mpifx_allreduce(env%globalMpiComm, tmpValsPrime%vals(iSys)%array,&
          & this%valsPrime%vals(iSys)%array, MPI_SUM)
    end do
  #:else
    this%valsPrime = tmpValsPrime
  #:endif

    if (tLead .and. this%tZscore) then
      call this%applyZscorePrime()
    end if

  #:if WITH_MPI
    do iSys = 1, size(this%valsPrime%vals)
      call mpifx_bcast(env%globalMpiComm, this%valsPrime%vals(iSys)%array)
    end do
  #:endif

  end subroutine TAcsf_calculatePrime


  !> Reduces a geometry to a subset of species as well as a given single atom.
  pure subroutine reduceGeometrySpecies(geo, iAtom, speciesIn, speciesOut, reducedGeo, iAtomOut,&
      & tKeep)

    !> geometry to reduce to a subset of species
    type(TGeometry), intent(in) :: geo

    !> index of an additional atom to keep
    integer, intent(in) :: iAtom

    !> species indices of current geometry
    integer, intent(in) :: speciesIn(:)

    !> desired species indices of output geometry
    integer, intent(in) :: speciesOut(:)

    !> reduced geometry
    type(TGeometry), intent(out) :: reducedGeo

    !> index of single atom after reduction
    integer, intent(out) :: iAtomOut

    !> mask to determine which atoms to keep (does also contain iAtom)
    integer, intent(out), allocatable, optional :: tKeep(:)

    !> temporary mask to determine which atoms to keep
    integer, allocatable :: tmp(:)

    !> auxiliary variables
    integer :: iAtom1, iAtom2, ind

    allocate(tmp(geo%nAtom))
    ind = 0

    outer: do iAtom1 = 1, geo%nAtom
      if (iAtom1 == iAtom) then
        ind = ind + 1
        tmp(ind) = iAtom1
        iAtomOut = ind
        cycle outer
      end if
      inner: do iAtom2 = 1, size(speciesOut)
        if (speciesIn(iAtom1) == speciesOut(iAtom2)) then
          ind = ind + 1
          tmp(ind) = iAtom1
          exit inner
        end if
      end do inner
    end do outer

    allocate(tKeep(ind))
    tKeep(:) = tmp(1:ind)

    reducedGeo = geo

    reducedGeo%coords = geo%coords(:, tKeep)
    reducedGeo%nAtom = size(reducedGeo%coords, dim=2)
    reducedGeo%species = geo%species(tKeep)
    call getNumberOfUniqueInt(reducedGeo%species, reducedGeo%nSpecies)

    ! remove atom to later build the neighborlist for
    ! ind = 0
    ! call move_alloc(tKeep, tmp)
    ! allocate(tKeep(size(tmp) - 1))
    ! do iAtom1 = 1, size(tmp)
    !   if (iAtom1 /= iAtomOut) then
    !     ind = ind + 1
    !     tKeep(ind) = tmp(iAtom)
    !   else
    !     cycle
    !   end if
    ! end do

  end subroutine reduceGeometrySpecies


  !> Calculates ACSF mappings for a single geometry.
  subroutine iGeoAcsf(this, geo, gFunctions, localAtToAtNum, extFeatures)

    !> symmetry function value instance of geometry
    type(TRealArray2D), intent(inout) :: this

    !> system geometry container
    type(TGeometry), intent(in) :: geo

    !> list of multi G-functions
    type(TGFunction), intent(in) :: gFunctions(:)

    !> index mapping local atom --> atomic number
    integer, intent(in) :: localAtToAtNum(:)

    !> atom dependent scaling parameters for cutoff function
    real(dp), intent(in) :: extFeatures(:,:)

    !> geometries reduced to atoms of a single species
    type(TGeometry) :: geo1, geo2

    !> atomic prefactor of central atom
    real(dp) :: atomId

    !> atomic prefactors of neighboring atoms
    real(dp), allocatable :: atomIds1(:), atomIds2(:)

    !> neighbor coordinates and squared distances
    real(dp), allocatable :: neighDists1(:), neighDists2(:), neighCoords1(:,:), neighCoords2(:,:)

    !> auxiliary variables
    integer :: iAtom, iAcsf, iAtomOut1

    do iAtom = 1, geo%nAtom
      do iAcsf = 1, size(gFunctions)

        call buildGFunctionNeighborlists(iAtom, geo, gFunctions(iAcsf), localAtToAtNum,&
            & extFeatures, geo1, geo2, iAtomOut1, atomIds1, atomIds2, neighDists1, neighDists2,&
            & neighCoords1, neighCoords2)

        if (gFunctions(iAcsf)%atomId > 0) then
          atomId = extFeatures(gFunctions(iAcsf)%atomId, iAtom)
        else
          atomId = 1.0_dp
        end if

        select case (gFunctions(iAcsf)%type)
        case ('g1')
          this%array(iAcsf, iAtom) = g1(neighDists1, atomId, atomIds1, gFunctions(iAcsf)%rCut)
        case ('g2')
          this%array(iAcsf, iAtom) = g2(neighDists1, atomId, atomIds1, gFunctions(iAcsf)%eta,&
              & gFunctions(iAcsf)%rs, gFunctions(iAcsf)%rCut)
        case ('g3')
          this%array(iAcsf, iAtom) = g3(neighDists1, atomId, atomIds1, gFunctions(iAcsf)%kappa,&
              & gFunctions(iAcsf)%rCut)
        case ('g4')
          this%array(iAcsf, iAtom) = g4(geo1%coords(:, iAtomOut1), neighCoords1, neighCoords2,&
              & neighDists1, neighDists2, atomId, atomIds1, atomIds2, gFunctions(iAcsf)%xi,&
              & gFunctions(iAcsf)%lambda, gFunctions(iAcsf)%eta, gFunctions(iAcsf)%rCut)
        case ('g5')
          this%array(iAcsf, iAtom) = g5(geo1%coords(:, iAtomOut1), neighCoords1, neighCoords2,&
              & neighDists1, neighDists2, atomId, atomIds1, atomIds2, gFunctions(iAcsf)%xi,&
              & gFunctions(iAcsf)%lambda, gFunctions(iAcsf)%eta, gFunctions(iAcsf)%rCut)
        case default
          call error('Invalid function type, aborting ACSF calculation.')
        end select

      end do
    end do

  end subroutine iGeoAcsf


  !> Calculates ACSF derivatives for a single geometry.
  subroutine iGeoAcsfGrad(this, geo, gFunctions, localAtToAtNum, extFeatures)

    !> symmetry function derivative instance of geometry
    type(TRealArray4D), intent(inout) :: this

    !> system geometry container
    type(TGeometry), intent(in) :: geo

    !> list of multiple G-functions
    type(TGFunction), intent(in) :: gFunctions(:)

    !> index mapping local atom --> atomic number
    integer, intent(in) :: localAtToAtNum(:)

    !> atom dependent scaling parameters for cutoff function
    real(dp), intent(in), optional :: extFeatures(:,:)

    !> geometries reduced to atoms of a single species
    type(TGeometry) :: geo1, geo2

    !> atomic prefactor of central atom
    real(dp) :: atomId

    !> atomic prefactors of neighboring atoms
    real(dp), allocatable :: atomIds1(:), atomIds2(:)

    !> neighbor coordinates and squared distances
    real(dp), allocatable :: neighDists1(:), neighDists2(:), neighCoords1(:,:), neighCoords2(:,:)

    !> global indices of neighboring atoms in the reduced geometry
    integer, allocatable  :: atomIndices1(:), atomIndices2(:)

    !> auxiliary variables
    integer :: iAtom, iAcsf, iAtomOut1

    do iAtom = 1, geo%nAtom
      do iAcsf = 1, size(gFunctions)

        call buildGFunctionNeighborlists(iAtom, geo, gFunctions(iAcsf), localAtToAtNum,&
            & extFeatures, geo1, geo2, iAtomOut1, atomIds1, atomIds2, neighDists1, neighDists2,&
            & neighCoords1, neighCoords2, atomIndices1_=atomIndices1, atomIndices2_=atomIndices2)

        if (gFunctions(iAcsf)%atomId > 0) then
          atomId = extFeatures(gFunctions(iAcsf)%atomId, iAtom)
        else
          atomId = 1.0_dp
        end if

        this%array(:, iAcsf, iAtom, atomIndices1) = gFuncGrad(gFunctions(iAcsf),&
            & geo%coords(:,iAtom), atomId, neighCoords1, neighDists1, atomIds1, neighCoords2,&
            & neighDists2, atomIds2)

        this%array(:, iAcsf, iAtom, iAtom) = - sum(this%array(:, iAcsf, iAtom, atomIndices1), dim=2)

        if (gFunctions(iAcsf)%tAngular) then
          if (gFunctions(iAcsf)%atomicNumbers(1) /= gFunctions(iAcsf)%atomicNumbers(2)) then
            ! Z1 <--> Z2
            this%array(:, iAcsf, iAtom, atomIndices2) = gFuncGrad(gFunctions(iAcsf),&
                & geo%coords(:,iAtom), atomId, neighCoords2, neighDists2, atomIds2, neighCoords1,&
                & neighDists1, atomIds1)

            this%array(:, iAcsf, iAtom, iAtom) = this%array(:, iAcsf, iAtom, iAtom)&
                & - sum(this%array(:, iAcsf, iAtom, atomIndices2), dim=2)
          end if
        end if

      end do
    end do

  end subroutine iGeoAcsfGrad


  !> Builds the G-function specific neighborlist for a given atom and geometry.
  subroutine buildGFunctionNeighborlists(iAtom, geo, gFunction, localAtToAtNum, extFeatures, geo1,&
      & geo2, iAtomOut1, atomIds1, atomIds2, neighDists1, neighDists2, neighCoords1, neighCoords2,&
      & atomIndices1_, atomIndices2_)

    !> atom index to calculate ACSF mapping force
    integer, intent(in) :: iAtom

    !> system geometry container
    type(TGeometry), intent(in) :: geo

    !> G-function parametrization
    type(TGFunction), intent(in) :: gFunction

    !> index mapping local atom --> atomic number
    integer, intent(in) :: localAtToAtNum(:)

    !> atom dependent scaling parameters for cutoff function
    real(dp), intent(in) :: extFeatures(:,:)

    !> geometries reduced to atoms of a single species
    type(TGeometry), intent(out) :: geo1, geo2

    !> index of iAtom after reduction of geo1
    integer, intent(out) :: iAtomOut1

    ! !> atomic prefactors of central atom
    ! real(dp), intent(out) :: atomId1, atomId2

    !> atomic prefactors of neighboring atoms
    real(dp), intent(out), allocatable :: atomIds1(:), atomIds2(:)

    !> squared distances of neighbors
    real(dp), intent(out), allocatable :: neighDists1(:), neighDists2(:)

    !> neighbor coordinates
    real(dp), intent(out), allocatable :: neighCoords1(:,:), neighCoords2(:,:)

    !> neighbour atom indices (dummy)
    integer,  intent(out), allocatable, optional :: atomIndices1_(:), atomIndices2_(:)

    !> neighbor atom indices
    integer,  allocatable :: atomIndices1(:), atomIndices2(:)

    !> masks to determine which atoms to keep while reducing the geometry
    integer, allocatable :: tKeep1(:), tKeep2(:)

    !> true, if a species-resolved neighborlist is desired
    logical :: tSpeciesResolved

    !> auxiliary variable
    integer :: iAtomOut2

    if (all(gFunction%atomicNumbers == 0)) then
      tSpeciesResolved = .false.
    else
      tSpeciesResolved = .true.
    end if

    if (.not. tSpeciesResolved) then
      geo1 = geo
      geo2 = geo1
      iAtomOut1 = iAtom
      iAtomOut2 = iAtomOut1
      if (gFunction%atomId > 0) then
        atomIds1 = extFeatures(gFunction%atomId, :)
        atomIds2 = atomIds1
      else
        if (allocated(atomIds1)) deallocate(atomIds1)
        allocate(atomIds1(geo1%nAtom))
        atomIds1(:) = 1.0_dp
        if (allocated(atomIds2)) deallocate(atomIds2)
        allocate(atomIds2(geo2%nAtom))
        atomIds2(:) = 1.0_dp
      end if
      call buildNeighborlist(geo1, gFunction%rCut, iAtom, neighDists1, neighCoords1, atomIndices1)
      neighDists2 = neighDists1
      neighCoords2 = neighCoords1
      atomIndices2 = atomIndices1
      atomIds1 = atomIds1(atomIndices1)
      atomIds2 = atomIds1

      if (present(atomIndices1_)) atomIndices1_ = atomIndices1
      if (present(atomIndices2_)) atomIndices2_ = atomIndices2

    else
      ! reduce geometry to a given species
      call reduceGeometrySpecies(geo, iAtom, localAtToAtNum, [gFunction%atomicNumbers(1)], geo1,&
          & iAtomOut1, tKeep=tKeep1)
      if (gFunction%atomId > 0) then
        atomIds1 = extFeatures(gFunction%atomId, tKeep1)
       else
        if (allocated(atomIds1)) deallocate(atomIds1)
        allocate(atomIds1(geo1%nAtom))
        atomIds1(:) = 1.0_dp
      end if
      call buildNeighborlist(geo1, gFunction%rCut, iAtomOut1, neighDists1, neighCoords1,&
          & atomIndices1)
      atomIds1 = atomIds1(atomIndices1)

      ! for angular functions species tupels, i.e. pairs of atomic numbers are required
      if (gFunction%tAngular) then
        call reduceGeometrySpecies(geo, iAtom, localAtToAtNum, [gFunction%atomicNumbers(2)], geo2,&
            & iAtomOut2, tKeep=tKeep2)
        if (gFunction%atomId > 0) then
          atomIds2 = extFeatures(gFunction%atomId, tKeep2)
         else
          if (allocated(atomIds2)) deallocate(atomIds2)
          allocate(atomIds2(geo2%nAtom))
          atomIds2(:) = 1.0_dp
        end if
        call buildNeighborlist(geo2, gFunction%rCut, iAtomOut2, neighDists2, neighCoords2,&
            & atomIndices2)
        atomIds2 = atomIds2(atomIndices2)
       else
        atomIds2 = atomIds1
        atomIndices2 = atomIndices1
        tKeep2 = tKeep1
      end if

      if (present(atomIndices1_)) atomIndices1_ = tKeep1(atomIndices1)
      if (present(atomIndices2_)) atomIndices2_ = tKeep2(atomIndices2)

    end if

  end subroutine buildGFunctionNeighborlists


  !> Builds the neighborlist of a given geometry and returns neighbor coordinates and distances.
  subroutine buildNeighborlist(geo, rCut, iAtom, allNeighDists, allNeighCoords, allAtomIndices)

    !> system geometry container
    type(TGeometry), intent(in) :: geo

    !> cutoff defining the neighborlist sphere
    real(dp), intent(in) :: rCut

    !> atom index to build neighbor list for
    integer, intent(in) :: iAtom

    !> neighbor coordinates and distances of all iterations
    real(dp), intent(out), allocatable :: allNeighDists(:), allNeighCoords(:,:)

    !> neighbor atom indices of all iterations
    integer, intent(out), allocatable :: allAtomIndices(:)

    !> instance of dynamic neighbour list
    type(TDynNeighList), target :: neighList

    !> pointer to dynamic neighbour list
    type(TDynNeighList), pointer :: pNeighList

    !> instance of dynamic neighbour list iterator
    type(TNeighIterator) :: neighIter

    !> obtained neighbor coordinates and distances
    real(dp) :: neighDists(iterChunkSize), neighCoords(3, iterChunkSize)

    !> obtained neighbor atom indices
    integer :: atomIndices(iterChunkSize)

    !> temporary coordinate storage
    real(dp), allocatable :: tmpCoords(:,:)

    !> auxiliary variables
    integer :: nNeigh, nTotNeigh

    call TDynNeighList_init(neighList, rCut, geo%nAtom, geo%tPeriodic)
    call neighList%updateCoords(geo%coords)

    if (geo%tPeriodic) then
      call neighList%updateLatVecs(geo%latVecs, geo%recVecs2p)
    end if
    pNeighList => neighList

    allocate(tmpCoords(3, 0))
    allocate(allNeighDists(0))
    allocate(allAtomIndices(0))

    call TNeighIterator_init(neighIter, pNeighList, iAtom, includeSelf=.false.)

    nNeigh = iterChunkSize

    do while (nNeigh == iterChunkSize)
      call neighIter%getNextNeighbours(nNeigh, coords=neighCoords, dists=neighDists,&
          & atomIndices=atomIndices)
      nTotNeigh = size(tmpCoords, dim=2) + nNeigh
      allocate(allNeighCoords(3, nTotNeigh))
      allNeighDists = [allNeighDists, neighDists(1:nNeigh)]
      allAtomIndices = [allAtomIndices, atomIndices(1:nNeigh)]
      allNeighCoords(:, 1:nTotNeigh - nNeigh) = tmpCoords
      allNeighCoords(:, nTotNeigh - nNeigh + 1:nTotNeigh) = neighCoords(:, 1:nNeigh)
      call move_alloc(allNeighCoords, tmpCoords)
    end do

    allNeighCoords = tmpCoords

  end subroutine buildNeighborlist


  !> Calculates the distance between two atoms of a geometry.
  pure function distance(geo, iAt1, iAt2)

    !> system geometry container
    type(TGeometry), intent(in) :: geo

    !> Indices of atoms to calculate distance for
    integer, intent(in) :: iAt1, iAt2

    !> Obtained distance betweeen the two atoms
    real(dp) :: distance

    distance = norm2(geo%coords(:, iAt2) - geo%coords(:, iAt1))

  end function distance


  !> Calculates the angle defined by three different atomic centers.
  pure function theta(atomCoords, neighCoords1, neighCoords2, neighDist1, neighDist2)

    !> coordinates of reference atom
    real(dp), intent(in) :: atomCoords(:)

    !> coordinates of two neighboring atoms
    real(dp), intent(in) :: neighCoords1(:), neighCoords2(:)

    !> distances of reference atom to first and second neighbor
    real(dp), intent(in) :: neighDist1, neighDist2

    !> obtained angle betweeen the three atoms
    real(dp) :: theta

    theta = acos(dot_product((neighCoords1 - atomCoords), (neighCoords2 - atomCoords))&
        & / (neighDist1 * neighDist2 + 1e-13_dp))

  end function theta


  !> Calculates the cutoff function for a given interatomic distance.
  !! Explicitly checks for rr < rcut.
  pure function cutoff(rr, atomId1, atomId2, rcut) result(res)

    !> atom distance
    real(dp), intent(in) :: rr

    !> atom ID of center atom
    real(dp), intent(in) :: atomId1

    !> atom ID of neighbor
    real(dp), intent(in) :: atomId2

    !> cutoff radius
    real(dp), intent(in) :: rcut

    !> resulting cutoff function value
    real(dp) :: res

    if (rr > rcut) then
      res = 0.0_dp
    else
      res = 0.5_dp * atomId1 * atomId2 * (cos(pi * rr / rcut) + 1.0_dp)
    end if

  end function cutoff


  !> Calculates the cutoff function for a given interatomic distance.
  !! Assumes that rr < rcut. For performance reasons, there is no explicit check!
  pure function cutoffWoCheck(rr, atomId1, atomId2, rcut) result(res)

    !> atom distance (in cutoff range)
    real(dp), intent(in) :: rr

    !> atom ID of center atom
    real(dp), intent(in) :: atomId1

    !> atom ID of neighbor
    real(dp), intent(in) :: atomId2

    !> cutoff radius
    real(dp), intent(in) :: rcut

    !> resulting cutoff function value
    real(dp) :: res

    res = 0.5_dp * atomId1 * atomId2 * (cos(pi * rr / rcut) + 1.0_dp)

  end function cutoffWoCheck


  !> Calculates the n-th derivative of the cutoff function.
  pure function dfCutoffWoCheck(rr, atomId1, atomId2, rcut, deriv) result(res)

    !> atom distance (in cutoff range)
    real(dp), intent(in) :: rr

    !> atom ID of center atom
    real(dp), intent(in) :: atomId1

    !> atom ID of neighbor
    real(dp), intent(in) :: atomId2

    !> cutoff radius
    real(dp), intent(in) :: rcut

    !> derivative order
    integer, intent(in) :: deriv

    !> resulting cutoff derivative value
    real(dp) :: res

    ! here we use the formula for the n-th derivative of cos(ax)
    res = 0.5_dp * atomId1 * atomId2 * (pi / rcut)**deriv&
        & * cos((pi * rr / rcut) + 0.5_dp * deriv * pi)

  end function dfCutoffWoCheck


  !> Calculates the cutoff function for an 1d array of interatomic distances.
  pure function cutoff1d(rr, atomId, atomIds, rcut) result(res)

    !> array of atom distances
    real(dp), intent(in) :: rr(:)

    !> atom ID of center atom
    real(dp), intent(in) :: atomId

    !> atom ID's of all neighbors
    real(dp), intent(in) :: atomIds(:)

    !> cutoff radius
    real(dp), intent(in) :: rcut

    !> corresponding activation function values
    real(dp) :: res(size(rr))

    where (rr > rcut)
      res(:) = 0.0_dp
    elsewhere
      res(:) = 0.5_dp * atomId * atomIds * (cos(pi * rr / rcut) + 1.0_dp)
    end where

  end function cutoff1d


  !> Calculates the G1 function.
  pure function g1(rr, atomId, atomIds, rcut)

    !> array of atom distances in cutoff range
    real(dp), intent(in) :: rr(:)

    !> atom ID of center atom
    real(dp), intent(in) :: atomId

    !> atom ID's of all neighbors
    real(dp), intent(in) :: atomIds(:)

    !> cutoff radius
    real(dp), intent(in) :: rcut

    !> corresponding symmetry function values
    real(dp) :: g1

    if (size(rr) == 0) then
      g1 = 0.0_dp
    else
      g1 = sum(cutoff1d(rr, atomId, atomIds, rcut))
    end if

  end function g1


  !> Calculates the G2 function.
  pure function g2(rr, atomId, atomIds, eta, rs, rcut)

    !> array of atom distances in cutoff range
    real(dp), intent(in) :: rr(:)

    !> atom ID of center atom
    real(dp), intent(in) :: atomId

    !> atom ID's of all neighbors
    real(dp), intent(in) :: atomIds(:)

    !> width of Gaussians
    real(dp), intent(in) :: eta

    !> center of Gaussians
    real(dp), intent(in) :: rs

    !> cutoff radius
    real(dp), intent(in) :: rcut

    !> corresponding symmetry function values
    real(dp) :: g2

    if (size(rr) == 0) then
      g2 = 0.0_dp
    else
      g2 = sum(exp(- eta * (rr - rs)**2) * cutoff1d(rr, atomId, atomIds, rcut))
    end if

  end function g2


  !> Calculates the G3 function.
  pure function g3(rr, atomId, atomIds, kappa, rcut)

    !> array of atom distances in cutoff range
    real(dp), intent(in) :: rr(:)

    !> atom ID of center atom
    real(dp), intent(in) :: atomId

    !> atom ID's of all neighbors
    real(dp), intent(in) :: atomIds(:)

    !> period length of damped cosine functions
    real(dp), intent(in) :: kappa

    !> cutoff radius
    real(dp), intent(in) :: rcut

    !> corresponding symmetry function values
    real(dp) :: g3

    if (size(rr) == 0) then
      g3 = 0.0_dp
    else
      g3 = sum(cos(kappa * rr) * cutoff1d(rr, atomId, atomIds, rcut))
    end if

  end function g3


  !> Calculates the G4 function.
  pure function g4(atomCoords, neighCoords1, neighCoords2, neighDists1, neighDists2, atomId,&
      & atomIds1, atomIds2, xi, lambda, eta, rcut)

    !> coordinates of reference atom
    real(dp), intent(in) :: atomCoords(:)

    !> coordinates of neighboring atoms of two species
    real(dp), intent(in) :: neighCoords1(:,:), neighCoords2(:,:)

    !> distances of reference atom to neighboring atoms
    real(dp), intent(in) :: neighDists1(:), neighDists2(:)

    !> atom ID of center atom
    real(dp), intent(in) :: atomId

    !> atom ID's of all neighbors
    real(dp), intent(in) :: atomIds1(:), atomIds2(:)

    !> parameter specifying angular resolution
    real(dp), intent(in) :: xi

    !> lambda parameter of angular part, recommended values: [+1, -1]
    real(dp), intent(in) :: lambda

    !> width of Gaussian part
    real(dp), intent(in) :: eta

    !> cutoff radius
    real(dp), intent(in) :: rcut

    ! distance between atom j of species 1 and atom k of species 2
    real(dp) :: distjk

    !> corresponding symmetry function values
    real(dp) :: g4

    !> auxiliary variables
    integer :: jj, kk

    g4 = 0.0_dp

    if ((size(neighDists1) == 0) .and. (size(neighDists2) == 0)) return

    do jj = 1, size(neighCoords1, dim=2)
      do kk = 1, size(neighCoords2, dim=2)
        ! calculate distance between atom j of species 1 and atom k of species 2
        distjk = norm2(neighCoords1(:, jj) - neighCoords2(:, kk))
        ! calculate G4 function contribution
        g4 = g4 + (1.0_dp + lambda * cos(theta(atomCoords, neighCoords1(:, jj),&
            & neighCoords2(:, kk), neighDists1(jj), neighDists2(kk))))**xi&
            & * exp(- eta * (neighDists1(jj)**2 + neighDists2(kk)**2 + distjk**2))&
            & * cutoffWoCheck(neighDists1(jj), atomId, atomIds1(jj), rcut)&
            & * cutoffWoCheck(neighDists2(kk), atomId, atomIds2(kk), rcut)&
            & * cutoff(distjk, atomIds1(jj), atomIds2(kk), rcut)
      end do
    end do

    g4 = g4 * 2.0_dp**(1.0_dp - xi)

  end function g4


  !> Calculates the G5 function.
  pure function g5(atomCoords, neighCoords1, neighCoords2, neighDists1, neighDists2, atomId,&
      & atomIds1, atomIds2, xi, lambda, eta, rcut)

    !> coordinates of reference atom
    real(dp), intent(in) :: atomCoords(:)

    !> coordinates of neighboring atoms of two species
    real(dp), intent(in) :: neighCoords1(:,:), neighCoords2(:,:)

    !> distances of reference atom to neighboring atoms
    real(dp), intent(in) :: neighDists1(:), neighDists2(:)

    !> atom ID of center atom
    real(dp), intent(in) :: atomId

    !> atom ID's of all neighbors
    real(dp), intent(in) :: atomIds1(:), atomIds2(:)

    !> parameter specifying angular resolution
    real(dp), intent(in) :: xi

    !> lambda parameter of angular part, recommended values: [+1, -1]
    real(dp), intent(in) :: lambda

    !> width of Gaussian part
    real(dp), intent(in) :: eta

    !> cutoff radius
    real(dp), intent(in) :: rcut

    !> corresponding symmetry function values
    real(dp) :: g5

    !> auxiliary variables
    integer :: jj, kk

    g5 = 0.0_dp

    if ((size(neighDists1) == 0) .and. (size(neighDists2) == 0)) return

    do jj = 1, size(neighCoords1, dim=2)
      do kk = 1, size(neighCoords2, dim=2)
        g5 = g5 + (1.0_dp + lambda * cos(theta(atomCoords, neighCoords1(:, jj),&
            & neighCoords2(:, kk), neighDists1(jj), neighDists2(kk))))**xi&
            & * exp(- eta * (neighDists1(jj)**2 + neighDists2(kk)**2))&
            & * cutoffWoCheck(neighDists1(jj), atomId, atomIds1(jj), rcut)&
            & * cutoffWoCheck(neighDists2(kk), atomId, atomIds2(kk), rcut)
      end do
    end do

    g5 = g5 * 2.0_dp**(1.0_dp - xi)

  end function g5


  !> Calculates G-function specific derivative for given neighborhoods.
  function gFuncGrad(gFunction, iCoords, atomId, neighCoords1, neighDists1, atomIds1, neighCoords2,&
      & neighDists2, atomIds2)

    !> list of multiple G-functions
    type(TGFunction), intent(in) :: gFunction

    !> cartesian coordinates of center atom
    real(dp), intent(in) :: iCoords(:)

    !> atom ID of center atom
    real(dp), intent(in) :: atomId

    !> atomic prefactors, neighbor coordinates and squared distances of first neighborlist
    real(dp), intent(in) :: atomIds1(:), neighDists1(:), neighCoords1(:,:)

    !> atomic prefactors, neighbor coordinates and squared distances of second neighborlist
    real(dp), intent(in) :: atomIds2(:), neighDists2(:), neighCoords2(:,:)

    !> xyz-component of G-function derivative for every neighboring atom
    real(dp) :: gFuncGrad(3, size(neighDists1))

    !> cutoff (derivative) values
    real(dp) :: dfc_ik, dfc_jk, fc_ik, fc_jk

    !> auxiliary variables
    real(dp) :: dist_jk, a_ijk, prefct_prod, T_ik(3, 3), uvec_jk(3), uvec_ik(3)
    real(dp), allocatable :: uvecs2(:,:), fc2(:), sqDists2(:)
    integer  :: jAtom, kAtom

    gFuncGrad(:,:) = 0.0_dp

    if (gFunction%type == 'g1') then
      do jAtom = 1, size(neighDists1)
        gFuncGrad(:, jAtom) = dfCutoffWoCheck(neighDists1(jAtom), atomId, atomIds1(jAtom),&
            & gFunction%rCut, 1) * (neighCoords1(:, jAtom) - iCoords) / neighDists1(jAtom)
      end do
      return
    end if

    if (gFunction%type == 'g2') then
      do jAtom = 1, size(neighDists1)
        gFuncGrad(:, jAtom) = (dfCutoffWoCheck(neighDists1(jAtom), atomId, atomIds1(jAtom),&
            & gFunction%rCut, 1) - cutoffWoCheck(neighDists1(jAtom), atomId, atomIds1(jAtom),&
            & gFunction%rCut) * 2.0_dp * gFunction%eta * (neighDists1(jAtom) - gFunction%rs))&
            & * exp(- gFunction%eta * (neighDists1(jAtom) - gFunction%rs)**2)&
            & * (neighCoords1(:, jAtom) - iCoords) / neighDists1(jAtom)
      end do
      return
    end if

    if (gFunction%type == 'g3') then
      associate (kappa => gFunction%kappa)
        do jAtom = 1, size(neighDists1)
          gFuncGrad(:, jAtom) = (cos(kappa * neighDists1(jAtom))&
              & * dfCutoffWoCheck(neighDists1(jAtom), atomId, atomIds1(jAtom), gFunction%rCut, 1)&
              & - sin(kappa * neighDists1(jAtom)) *  cutoffWoCheck(neighDists1(jAtom), atomId,&
              & atomIds1(jAtom), gFunction%rCut) * kappa) * (neighCoords1(:, jAtom) - iCoords)&
              & / neighDists1(jAtom)
        end do
      end associate
      return
    end if

    allocate(uvecs2, source=neighCoords2)
    allocate(sqDists2, source=neighDists2)
    allocate(fc2, mold=neighDists2)
    fc2(:) = 0.0_dp

    do jAtom = 1, size(neighDists2)
        uvecs2(:, jAtom) = (uvecs2(:, jAtom) - iCoords) / neighDists2(jAtom)
        sqDists2(jAtom) = neighDists2(jAtom)**2
        if (neighDists2(jAtom) > gFunction%rCut) cycle
        fc2(jAtom) = cutoffWoCheck(neighDists2(jAtom), atomId, atomIds2(jAtom), gFunction%rCut)
    end do

    ! first handle the case where both neighbourhoods are the same
    ! this happens when there's no species resolution or when Z1 == Z2
    if (gFunction%atomicNumbers(1) == gFunction%atomicNumbers(2)) then

      do kAtom = 1, size(neighDists1)
        T_ik = T_matr(uvecs2(:, kAtom), neighDists1(kAtom)) * gFunction%lambda * gFunction%xi
        dfc_ik = - 2.0_dp * gFunction%eta * neighDists1(kAtom) * fc2(kAtom)&
            & + dfCutoffWoCheck(neighDists1(kAtom), atomId, atomIds1(kAtom), gFunction%rCut, 1)

        ! case kAtom == jAtom
        prefct_prod = 1.0_dp
        if (gFunction%type == 'g4') prefct_prod = atomIds1(kAtom) * atomIds2(kAtom)
        gFuncGrad(:, kAtom) = gFuncGrad(:, kAtom) + (1.0_dp + gFunction%lambda)**gFunction%xi&
            & * fc2(kAtom) * dfc_ik * exp(- 2.0_dp * gFunction%eta * sqDists2(kAtom))&
            & * uvecs2(:, kAtom) * prefct_prod

        do jAtom = 1, size(neighDists2)
          ! skip jAtom == kAtom
          if (jAtom == kAtom) cycle

          a_ijk = 1.0_dp + gFunction%lambda * dot_product(uvecs2(:, jAtom), uvecs2(:, kAtom))

          if (gFunction%type == 'g5') then
            gFuncGrad(:, kAtom) = gFuncGrad(:, kAtom)&
                & + fc2(jAtom) * exp(- gFunction%eta * (sqDists2(kAtom) + sqDists2(jAtom)))&
                & * a_ijk**(gFunction%xi - 1.0_dp) * (fc2(kAtom) * matmul(T_ik, uvecs2(:, jAtom))&
                & + a_ijk * dfc_ik * uvecs2(:, kAtom))
            cycle
          end if

          uvec_jk = neighCoords1(:, kAtom) - neighCoords2(:, jAtom)
          dist_jk = norm2(uvec_jk)

          if (dist_jk > gFunction%rCut) cycle

          uvec_jk = uvec_jk / dist_jk
          fc_jk   = cutoffWoCheck(dist_jk, AtomIds1(kAtom), atomIds2(jAtom), gFunction%rCut)
          dfc_jk  = - 2.0_dp * gFunction%eta * dist_jk * fc_jk&
              & + dfCutoffWoCheck(dist_jk, AtomIds1(kAtom), atomIds2(jAtom), gFunction%rCut, 1)

          gFuncGrad(:, kAtom) = gFuncGrad(:, kAtom) + &
              & fc2(jAtom) * exp(- gFunction%eta * (sqDists2(kAtom) + sqDists2(jAtom)&
              & + dist_jk**2)) * a_ijk**(gFunction%xi - 1.0_dp) * (fc2(kAtom) * fc_jk&
              & * matmul(T_ik, uvecs2(:, jAtom)) + a_ijk * (fc_jk * dfc_ik * uvecs2(:, kAtom)&
              & + fc2(kAtom) * dfc_jk * uvec_jk))
        end do

      end do
      gFuncGrad = gFuncGrad * 2.0_dp**(2.0_dp - gFunction%xi)
      return
    end if

    ! finally, case Z1 /= Z2
    do kAtom = 1, size(neighDists1)

      fc_ik = cutoffWoCheck(neighDists1(kAtom), atomId, atomIds1(kAtom), gFunction%rCut)
      dfc_ik = - 2.0_dp * gFunction%eta * neighDists1(kAtom) * fc_ik&
          & + dfCutoffWoCheck(neighDists1(kAtom), atomId, atomIds1(kAtom), gFunction%rCut, 1)

      uvec_ik = (neighCoords1(:, kAtom) - iCoords) / neighDists1(kAtom)
      T_ik = T_matr(uvec_ik, neighDists1(kAtom)) * gFunction%lambda * gFunction%xi

      do jAtom = 1, size(neighDists2)
        a_ijk = 1.0_dp + gFunction%lambda * dot_product(uvecs2(:, jAtom), uvec_ik)

        if (gFunction%type == 'g5') then
          gFuncGrad(:, kAtom) = gFuncGrad(:, kAtom) + fc2(jAtom) * exp(- gFunction%eta&
              & * (neighDists1(kAtom)**2 + sqDists2(jAtom))) * a_ijk**(gFunction%xi - 1.0_dp)&
              & * (fc_ik * matmul(T_ik, uvecs2(:, jAtom)) + a_ijk * dfc_ik * uvec_ik)
          cycle
        end if

        uvec_jk = neighCoords1(:, kAtom) - neighCoords2(:, jAtom)
        dist_jk = norm2(uvec_jk)

        if (dist_jk > gFunction%rCut) cycle

        uvec_jk = uvec_jk / dist_jk
        fc_jk = cutoffWoCheck(dist_jk, AtomIds1(kAtom), atomIds2(jAtom), gFunction%rCut)
        dfc_jk = - 2.0_dp * gFunction%eta * dist_jk * fc_jk + dfCutoffWoCheck(dist_jk,&
            & AtomIds1(kAtom), atomIds2(jAtom), gFunction%rCut, 1)

        gFuncGrad(:, kAtom) = gFuncGrad(:, kAtom) + fc2(jAtom) * exp(- gFunction%eta&
            & * (neighDists1(kAtom)**2 + sqDists2(jAtom) + dist_jk**2))&
            & * a_ijk**(gFunction%xi - 1.0_dp) * (fc_ik * fc_jk * matmul(T_ik, uvecs2(:, jAtom))&
            & + a_ijk * (fc_jk * dfc_ik * uvec_ik + fc_ik * dfc_jk * uvec_jk))
      end do
    end do

    gFuncGrad(:,:) = gFuncGrad * 2.0_dp**(1.0_dp - gFunction%xi)

  end function gFuncGrad


  !>
  pure function T_matr(uvec, nrm)

    !>
    real(dp), intent(in)  :: uvec(3), nrm

    !>
    real(dp) :: T_matr(3, 3)

    !> auxiliary variables
    integer  :: ii, jj

    T_matr(:,:) = 0.0_dp

    do ii = 1, 3
      T_matr(ii, ii) = 1.0_dp - uvec(ii)**2
      do jj = ii + 1, 3
        T_matr(ii, jj) = - uvec(ii) * uvec(jj)
        T_matr(jj, ii) = T_matr(ii, jj)
      end do
    end do

    T_matr(:,:) = T_matr / nrm

  end function


#:if WITH_MPI

  !> Synchronizes the ACSF derived type, apart from raw values.
  subroutine TAcsf_syncConfig(this, comm)

    !> representation of ACSF mappings
    class(TAcsf), intent(inout) :: this

    !> mpi communicator with some additional information
    type(mpifx_comm), intent(in) :: comm

    !> auxiliary variables
    integer :: iFunc, dims0d

    call mpifx_bcast(comm, this%tZscore)

    if (comm%lead) then
      dims0d = size(this%gFunctions%func)
    end if
    call mpifx_bcast(comm, dims0d)

    if (.not. comm%lead) then
      if (allocated(this%gFunctions%func)) deallocate(this%gFunctions%func)
      allocate(this%gFunctions%func(dims0d))
    end if

    ! synchronize G-function parameters
    do iFunc = 1, size(this%gFunctions%func)
      call mpifx_bcast(comm, this%gFunctions%func(iFunc)%type)
      call mpifx_bcast(comm, this%gFunctions%func(iFunc)%tRadial)
      call mpifx_bcast(comm, this%gFunctions%func(iFunc)%tAngular)
      call mpifx_bcast(comm, this%gFunctions%func(iFunc)%rCut)
      call mpifx_bcast(comm, this%gFunctions%func(iFunc)%kappa)
      call mpifx_bcast(comm, this%gFunctions%func(iFunc)%rs)
      call mpifx_bcast(comm, this%gFunctions%func(iFunc)%eta)
      call mpifx_bcast(comm, this%gFunctions%func(iFunc)%lambda)
      call mpifx_bcast(comm, this%gFunctions%func(iFunc)%xi)
      call mpifx_bcast(comm, this%gFunctions%func(iFunc)%atomId)
      call mpifx_bcast(comm, this%gFunctions%func(iFunc)%atomicNumbers)
    end do

    if (this%tZscore) then
      if (comm%lead) then
        dims0d = size(this%zPrec, dim=1)
      end if
      call mpifx_bcast(comm, dims0d)
      if (.not. comm%lead) then
        if (allocated(this%zPrec)) deallocate(this%zPrec)
        allocate(this%zPrec(dims0d, 2))
      end if
      call mpifx_bcast(comm, this%zPrec)
    end if

  end subroutine TAcsf_syncConfig

#:endif


  !> Writes the ACSF configuration to a HDF5 netstat file.
  subroutine TAcsf_toFile(this, fname)

    !> representation of ACSF mappings
    class(TAcsf), intent(in) :: this

    !> filename or path to write to
    character(len=*), intent(in) :: fname

    !> various specifier flags
    integer(hid_t) :: file_id, netstat_id, mapping_id, func_id, precond_id

    !> name of current G-function
    character(len=:), allocatable :: funcname

    !> number of symmetry mappings
    integer :: nFunctions

    !> auxiliary variable
    integer(size_t) :: dim

    !> auxiliary variables
    integer :: iFunc, iErr, tmp, tExist

    ! open the hdf5 interface
    call h5open_f(iErr)

    ! open the netstat file
    call h5fopen_f(fname, H5F_ACC_RDWR_F, file_id, iErr)

    ! open the netstat group
    call h5gopen_f(file_id, 'netstat', netstat_id, iErr)

    ! check if a mapping entry is already present
    tExist = h5ltfind_dataset_f(netstat_id, 'mapping')

    if (tExist == 1) then
      call h5ldelete_f(netstat_id, 'mapping', iErr)
    end if

    ! create the mapping group
    call h5gcreate_f(netstat_id, 'mapping', mapping_id, iErr)

    ! set the type attribute to acsf
    call h5ltset_attribute_string_f(mapping_id, './', 'type', 'acsf', iErr)

    nFunctions = size(this%gFunctions%func)

    ! set the total number of mapping functions
    dim = 1
    call h5ltset_attribute_int_f(mapping_id, './', 'nfunctions', [nFunctions], dim, iErr)

    do iFunc = 1, nFunctions
      funcname = 'function' // i2c(iFunc)

      ! create the function group
      call h5gcreate_f(mapping_id, funcname, func_id, iErr)

      ! set the G-function type
      call h5ltset_attribute_string_f(func_id, './', 'type', this%gFunctions%func(iFunc)%type, iErr)

      ! set the dataset index for atomID
      dim = 1
      call h5ltset_attribute_int_f(func_id, './', 'atomid', [this%gFunctions%func(iFunc)%atomId],&
          & dim, iErr)

      ! set the cutoff radius
      dim = 1
      call h5ltset_attribute_double_f(func_id, './', 'cutoff', [this%gFunctions%func(iFunc)%rCut],&
          & dim, iErr)

      ! set the atomic numbers for the species-resolved scheme
      dim = 2
      call h5ltset_attribute_int_f(func_id, './', 'atomicnumbers',&
          & this%gFunctions%func(iFunc)%atomicNumbers, dim, iErr)

      ! specify whether this is a radial function
      dim = 1
      if (this%gFunctions%func(iFunc)%tRadial) then
        tmp = 1
      else
        tmp = 0
      end if
      call h5ltset_attribute_int_f(func_id, './', 'tradial', [tmp], dim, iErr)

      ! specify whether this is an angular function
      dim = 1
      if (this%gFunctions%func(iFunc)%tAngular) then
        tmp = 1
      else
        tmp = 0
      end if
      call h5ltset_attribute_int_f(func_id, './', 'tangular', [tmp], dim, iErr)

      ! write G-function dependent parameters
      select case (this%gFunctions%func(iFunc)%type)
      case ('g2')
        dim = 1
        call h5ltset_attribute_double_f(func_id, './', 'eta', [this%gFunctions%func(iFunc)%eta],&
            & dim, iErr)
        call h5ltset_attribute_double_f(func_id, './', 'rs', [this%gFunctions%func(iFunc)%rs],&
            & dim, iErr)
      case ('g3')
        dim = 1
        call h5ltset_attribute_double_f(func_id, './', 'kappa',&
            & [this%gFunctions%func(iFunc)%kappa], dim, iErr)
      case ('g4', 'g5')
        dim = 1
        call h5ltset_attribute_double_f(func_id, './', 'xi', [this%gFunctions%func(iFunc)%xi],&
            & dim, iErr)
        call h5ltset_attribute_double_f(func_id, './', 'lambda',&
            & [this%gFunctions%func(iFunc)%lambda], dim, iErr)
        call h5ltset_attribute_double_f(func_id, './', 'eta', [this%gFunctions%func(iFunc)%eta],&
            & dim, iErr)
      end select

      ! close the function group
      call h5gclose_f(func_id, iErr)

    end do

    if (this%tZscore) then

      ! create the preconditioning group
      call h5gcreate_f(mapping_id, 'preconditioning', precond_id, iErr)

      ! set the type attribute to z-score
      call h5ltset_attribute_string_f(precond_id, './', 'type', 'zscore', iErr)

      ! write means
      call h5ltfxmake_dataset_double_f(precond_id, 'means', this%zPrec(:, 1))

      ! write variances
      call h5ltfxmake_dataset_double_f(precond_id, 'variances', this%zPrec(:, 2))

      ! close the preconditioning group
      call h5gclose_f(precond_id, iErr)

    end if

    ! close the mapping group
    call h5gclose_f(mapping_id, iErr)

    ! close the netstat group
    call h5gclose_f(netstat_id, iErr)

    ! close the netstat file
    call h5fclose_f(file_id, iErr)

    ! close the hdf5 interface
    call h5close_f(iErr)

  end subroutine TAcsf_toFile


  !> Reads the ACSF configuration from a HDF5 netstat file.
  subroutine TAcsf_fromFile(this, fname, tReduce, tStandardize, nRadial, nAngular)

    !> representation of ACSF mappings
    class(TAcsf), intent(out) :: this

    !> filename or path to read from
    character(len=*), intent(in) :: fname

    !> true, if the species-resolved features should get summed up
    logical, intent(out), optional :: tReduce

    !> true, if standardization of features is desired
    logical, intent(out), optional :: tStandardize

    !> number of radial and angular mappings
    integer, intent(out), optional :: nRadial, nAngular

    !> various specifier flags
    integer(hid_t) :: file_id, netstat_id, mapping_id, func_id, precond_id

    !> name of current G-function
    character(len=:), allocatable :: funcname

    !> temporary storage container
    integer, allocatable :: acsfAtomicNumbers(:), tmpAtomicNumbers(:)

    !> auxiliary variables
    character(len=100) :: tmpStr
    real(dp) :: tmpReal(1)
    real(dp), allocatable :: tmpRealArray1d(:)
    integer :: iFunc, iErr, tExist, tmpInt(1)

    ! open the hdf5 interface
    call h5open_f(iErr)

    ! open the netstat file
    call h5fopen_f(fname, H5F_ACC_RDONLY_F, file_id, iErr)

    ! open the netstat group
    call h5gopen_f(file_id, 'netstat', netstat_id, iErr)

    ! check if a mapping entry is present
    tExist = h5ltfind_dataset_f(netstat_id, 'mapping')

    if (.not. (tExist == 1)) then
      call error('Error while reading ACSF from file. Mapping group absent.')
    end if

    ! open the mapping group
    call h5gopen_f(netstat_id, 'mapping', mapping_id, iErr)

    ! check the mapping type
    call h5ltget_attribute_string_f(mapping_id, './', 'type', tmpStr, iErr)
    if (.not. (tolower(trim(tmpStr)) == 'acsf')) then
      call error('Error while reading ACSF from netstat file. Invalid mapping type '&
          & // tolower(trim(tmpStr)) // '.')
    end if

    ! read the total number of mapping functions
    call h5ltget_attribute_int_f(mapping_id, './', 'nfunctions', tmpInt, iErr)
    if (tmpInt(1) <= 0) then
      call error('Error while reading ACSF from netstat file. Non-positive number of ACSF '&
          & // 'functions obtained.')
    else
      allocate(this%gFunctions%func(tmpInt(1)))
    end if

    do iFunc = 1, size(this%gFunctions%func)
      funcname = 'function' // i2c(iFunc)

      ! open the function group
      call h5gopen_f(mapping_id, funcname, func_id, iErr)

      ! read the G-function type
      call h5ltget_attribute_string_f(func_id, './', 'type', tmpStr, iErr)
      if (.not. ((tolower(trim(tmpStr)) == 'g1') .or. (tolower(trim(tmpStr)) == 'g2')&
          & .or. (tolower(trim(tmpStr)) == 'g3') .or. (tolower(trim(tmpStr)) == 'g4')&
          & .or. (tolower(trim(tmpStr)) == 'g5'))) then
        call error('Error while reading ACSF from netstat file. Invalid G-Function type '&
            & // tolower(trim(tmpStr)) // '.')
      end if
      this%gFunctions%func(iFunc)%type = tolower(trim(tmpStr))

      ! read the dataset index for atomID
      call h5ltget_attribute_int_f(func_id, './', 'atomid', tmpInt, iErr)
      if (tmpInt(1) < 0) then
        call error('Error while reading ACSF from netstat file. Negative atomId obtained.')
      end if
      this%gFunctions%func(iFunc)%atomId = tmpInt(1)

      ! read the cutoff radius
      call h5ltget_attribute_double_f(func_id, './', 'cutoff', tmpReal, iErr)
      if (tmpReal(1) <= 0.0_dp) then
        call error('Error while reading ACSF from netstat file. Invalid cutoff obtained.')
      end if
      this%gFunctions%func(iFunc)%rCut = tmpReal(1)

      ! read the atomic numbers for the species-resolved scheme
      call h5ltget_attribute_int_f(func_id, './', 'atomicnumbers',&
          & this%gFunctions%func(iFunc)%atomicNumbers, iErr)
      if ((this%gFunctions%func(iFunc)%atomicNumbers(1) < 0)&
          & .or. (this%gFunctions%func(iFunc)%atomicNumbers(2) < 0)) then
        call error('Error while reading ACSF from netstat file. Invalid atomic numbers obtained.')
      end if

      ! inquire whether this is a radial function
      call h5ltget_attribute_int_f(func_id, './', 'tradial', tmpInt, iErr)
      if (tmpInt(1) == 0) then
        this%gFunctions%func(iFunc)%tRadial = .false.
      elseif (tmpInt(1) == 1) then
        this%gFunctions%func(iFunc)%tRadial = .true.
      else
        call error('Error while reading ACSF from netstat file. Invalid tRadial entry.')
      end if

      ! inquire whether this is an angular function
      call h5ltget_attribute_int_f(func_id, './', 'tangular', tmpInt, iErr)
      if (tmpInt(1) == 0) then
        this%gFunctions%func(iFunc)%tAngular = .false.
      elseif (tmpInt(1) == 1) then
        this%gFunctions%func(iFunc)%tAngular = .true.
      else
        call error('Error while reading ACSF from netstat file. Invalid tAngular entry.')
      end if

      ! check consistency of tRadial and tAngular entries
      if (this%gFunctions%func(iFunc)%tRadial .eqv. this%gFunctions%func(iFunc)%tAngular) then
        call error('Error while reading ACSF from netstat file. Conflicting entries tRadial '&
            & // 'and tAngular.')
      end if

      ! provide dummy values for now
      this%gFunctions%func(iFunc)%kappa = 0.0_dp
      this%gFunctions%func(iFunc)%rs = 0.0_dp
      this%gFunctions%func(iFunc)%eta = 0.0_dp
      this%gFunctions%func(iFunc)%lambda = 0.0_dp
      this%gFunctions%func(iFunc)%xi = 0.0_dp

      ! read G-function dependent parameters
      select case (this%gFunctions%func(iFunc)%type)
      case ('g2')
        call h5ltget_attribute_double_f(func_id, './', 'eta', tmpReal, iErr)
        this%gFunctions%func(iFunc)%eta = tmpReal(1)
        call h5ltget_attribute_double_f(func_id, './', 'rs', tmpReal, iErr)
        this%gFunctions%func(iFunc)%rs = tmpReal(1)
      case ('g3')
        call h5ltget_attribute_double_f(func_id, './', 'kappa', tmpReal, iErr)
        this%gFunctions%func(iFunc)%kappa = tmpReal(1)
      case ('g4', 'g5')
        call h5ltget_attribute_double_f(func_id, './', 'xi', tmpReal, iErr)
        this%gFunctions%func(iFunc)%xi = tmpReal(1)
        call h5ltget_attribute_double_f(func_id, './', 'lambda', tmpReal, iErr)
        this%gFunctions%func(iFunc)%lambda = tmpReal(1)
        call h5ltget_attribute_double_f(func_id, './', 'eta', tmpReal, iErr)
        this%gFunctions%func(iFunc)%eta = tmpReal(1)
      end select

      ! close the function group
      call h5gclose_f(func_id, iErr)

    end do

    ! inquire preconditioning configuration
    tExist = h5ltfind_dataset_f(mapping_id, 'preconditioning')

    if (tExist == 0) then
      this%tZscore = .false.
    else
      this%tZscore = .true.
    end if

    if (this%tZscore) then

      ! open the preconditioning group
      call h5gopen_f(mapping_id, 'preconditioning', precond_id, iErr)

      ! read the preconditioning type
      call h5ltget_attribute_string_f(precond_id, './', 'type', tmpStr, iErr)
      if (.not. (tolower(trim(tmpStr)) == 'zscore')) then
        call error('Error while reading ACSF from netstat file. Invalid preconditioning type '&
            & // tolower(trim(tmpStr)) // '.')
      end if

      if (allocated(this%zPrec)) deallocate(this%zPrec)
      allocate(this%zPrec(size(this%gFunctions%func), 2))

      ! read means and variances
      call h5ltfx_read_dataset_double_f(precond_id, 'means', tmpRealArray1d)
      this%zPrec(:, 1) = tmpRealArray1d
      call h5ltfx_read_dataset_double_f(precond_id, 'variances', tmpRealArray1d)
      this%zPrec(:, 2) = tmpRealArray1d

      ! close the preconditioning group
      call h5gclose_f(precond_id, iErr)

    end if

    ! close the mapping group
    call h5gclose_f(mapping_id, iErr)

    ! close the netstat group
    call h5gclose_f(netstat_id, iErr)

    ! close the netstat file
    call h5fclose_f(file_id, iErr)

    ! close the hdf5 interface
    call h5close_f(iErr)

    nRadial = 0
    nAngular = 0
    allocate(tmpAtomicNumbers(0))
    ! determine if the ACSF configuration is species-resolved
    do iFunc = 1, size(this%gFunctions%func)
      if (this%gFunctions%func(iFunc)%tRadial) then
        nRadial = nRadial + 1
        tmpAtomicNumbers = [tmpAtomicNumbers, [this%gFunctions%func(iFunc)%atomicNumbers(1)]]
      else
        nAngular = nAngular + 1
        tmpAtomicNumbers = [tmpAtomicNumbers, this%gFunctions%func(iFunc)%atomicNumbers]
      end if
    end do

    ! reduce the ACSF atomic numbers to unique ones
    call getUniqueInt(tmpAtomicNumbers, acsfAtomicNumbers)

    if (all(acsfAtomicNumbers == 0)) then
      tReduce = .true.
    else
      tReduce = .false.
    end if
    tStandardize = this%tZscore

  end subroutine TAcsf_fromFile

end module fnet_acsf
