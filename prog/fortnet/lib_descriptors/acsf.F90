!--------------------------------------------------------------------------------------------------!
!  FORTNET: A Behler-Parrinello-Neural-Network Implementation                                      !
!  Copyright (C) 2020 - 2021  T. W. van der Heide                                                  !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

module fnet_acsf

  use dftbp_accuracy, only: dp
  use dftbp_message, only : error
  use dftbp_constants, only: pi
  use dftbp_globalenv, only : stdOut
  use dftbp_typegeometry, only : TGeometry
  use dftbp_dynneighlist, only : TDynNeighList, TDynNeighList_init
  use dftbp_dynneighlist, only : TNeighIterator, TNeighIterator_init

  use fnet_nestedtypes, only : TIntArray1D, TRealArray1D, TRealArray2D, TEnv
  use fnet_workarounds, only : myFindloc

#:if WITH_MPI
  use fnet_mpifx
  use fnet_parallel, only : getStartAndEndIndex
#:endif

  implicit none

  private

  public :: TAcsf, TAcsf_init, TAcsfParams_init, TMultiAcsfVals


  type :: TG2Params

    !> width of Gaussians
    real(dp) :: eta

    !> centers of Gaussians
    real(dp), allocatable :: rs(:)

  end type TG2Params


  type :: TG3Params

    !> kappa scaling for G3 function
    real(dp) :: kappa

  end type TG3Params


  type :: TG4Params

    !> width of Gaussians
    real(dp) :: eta

    !> lambda parameter of angular part
    real(dp), allocatable :: lambda(:)

    !> parameter specifying angular resolution
    real(dp), allocatable :: xi(:)

  end type TG4Params


  type :: TAcsfParams

    !> parameters of G2 functions
    type(TG2Params) :: g2

    !> parameters of G3 functions
    type(TG3Params) :: g3

    !> parameters of G4 functions
    type(TG4Params) :: g4

  end type TAcsfParams


  type :: TMultiAcsfVals

    !> array of ACSF value instances
    type(TRealArray2D), allocatable :: vals(:)

  end type TMultiAcsfVals


  type :: TAcsf

    !> number of radial symmetry functions
    integer :: nRadial

    !> number of angular symmetry functions
    integer :: nAngular

    !> cutoff radius
    real(dp) :: rcut

    !> species dependent scaling parameters for cutoff function
    real(dp), allocatable :: speciesIds(:)

    !> real valued atom identifier, expected shape: [nDatapoints]
    type(TRealArray1D), allocatable :: atomIds(:)

    !> container for all ACSF function parameters
    type(TAcsfParams) :: param

    !> wrapper around multiple ACSF value instances
    type(TMultiAcsfVals) :: vals

    !> storage container of means and variances to calculate z-score
    real(dp), allocatable :: zPrec(:,:)

    !> true, if z-score standardization should be applied
    logical :: tZscore

  contains

    procedure :: getMeansAndVariances => TAcsf_getMeansAndVariances
    procedure :: applyZscore => TAcsf_applyZscore
    procedure :: calculate => TAcsf_calculate
    procedure :: fromFile => TAcsf_fromFile
    procedure :: toFile => TAcsf_toFile

  end type TAcsf

  !> Chunk size to use when obtaining neighbours dynamically via an iterator
  !> (used to artificially restrict memory usage to a certain amount)
  integer, parameter :: iterChunkSize = 1000


contains


  subroutine TAcsf_init(this, nRadial, nAngular, rcut, speciesIds, kk, tZscore)

    !> representation of ACSF mappings
    type(TAcsf), intent(inout) :: this

    !> number of radial symmetry functions
    integer, intent(in) :: nRadial

    !> number of angular symmetry functions
    integer, intent(in) :: nAngular

    !> cutoff radius
    real(dp), intent(in) :: rcut

    !> species dependent scaling parameters for cutoff function
    real(dp), intent(in) :: speciesIds(:)

    !> optional kappa parameter for G3 function
    real(dp), intent(in), optional :: kk

    !> true, if z-score standardization should be applied
    logical, intent(in), optional :: tZscore

    !> if present, equals given kappa value; otherwise 1.0
    real(dp) :: kappa

    if (present(kk)) then
      kappa = kk
    else
      kappa = 1.0_dp
    end if

    if (present(tZscore)) then
      this%tZscore = tZscore
    else
      this%tZscore = .false.
    end if

    this%nRadial = nRadial
    this%nAngular = nAngular
    this%rcut = rcut
    this%speciesIds = speciesIds

  end subroutine TAcsf_init


  subroutine TAcsf_toFile(this, fname, globalSpNames)

    !> representation of ACSF mappings
    class(TAcsf), intent(in) :: this

    !> filename or path to save acsf parameters to
    character(len=*), intent(in) :: fname

    !> contains (unique) species of all dataset geometries
    character(len=*), intent(in) :: globalSpNames(:)

    !> unique fileunit
    integer :: fd

    !> auxiliary variables
    integer :: iAcsf, iSp

    open(newunit=fd, file=fname, form='formatted', status='replace', action='write')

    write(fd, '(ES26.16E3)') this%rCut

    write(fd, '(I26,ES26.16E3)') this%nRadial, this%param%g2%eta
    do iAcsf = 1, this%nRadial
      write(fd, '(ES26.16E3)') this%param%g2%rs(iAcsf)
    end do

    write(fd, '(I26,ES26.16E3)') this%nAngular, this%param%g4%eta
    do iAcsf = 1, this%nAngular
      write(fd, '(2ES26.16E3)') this%param%g4%xi(iAcsf), this%param%g4%lambda(iAcsf)
    end do

    write(fd, '(I26)') size(this%speciesIds)
    do iSp = 1, size(this%speciesIds)
      write(fd, '(A,ES26.16E3)') trim(globalSpNames(iSp)), this%speciesIds(iSp)
    end do

    if (this%tZscore) then
      write(fd, '(I26)') 1
      do iAcsf = 1, this%nRadial + this%nAngular
        write(fd, '(2ES26.16E3)') this%zPrec(iAcsf, 1), this%zPrec(iAcsf, 2)
      end do
    else
      write(fd, '(I26)') 0
    end if

    close(fd)

  end subroutine TAcsf_toFile


  subroutine TAcsf_fromFile(this, fname, globalSpNames)

    !> representation of ACSF mappings
    class(TAcsf), intent(out) :: this

    !> filename or path to save acsf parameters to
    character(len=*), intent(in) :: fname

    !> contains (unique) species of all dataset geometries
    character(len=*), intent(in) :: globalSpNames(:)

    !> unique fileunit
    integer :: fd

    !> temporary storage container
    character(len=50) :: speciesNameTmp
    real(dp) :: speciesIdTmp

    !> auxiliary variables
    integer :: isZscore, nSpecies

    !> auxiliary variables
    integer :: iAcsf, iSp

    open(newunit=fd, file=fname, form='formatted', status='old', action='read')

    read(fd, *) this%rCut

    read(fd, *) this%nRadial, this%param%g2%eta
    allocate(this%param%g2%rs(this%nRadial))
    do iAcsf = 1, this%nRadial
      read(fd, *) this%param%g2%rs(iAcsf)
    end do

    read(fd, *) this%nAngular, this%param%g4%eta
    allocate(this%param%g4%xi(this%nAngular))
    allocate(this%param%g4%lambda(this%nAngular))
    do iAcsf = 1, this%nAngular
      read(fd, *) this%param%g4%xi(iAcsf), this%param%g4%lambda(iAcsf)
    end do

    read(fd, *) nSpecies
    allocate(this%speciesIds(nSpecies))
    do iSp = 1, size(this%speciesIds)
      read(fd, *) speciesNameTmp, speciesIdTmp
      if (.not. (myFindloc(globalSpNames, speciesNameTmp) == 0)) then
        this%speciesIds(myFindloc(globalSpNames, speciesNameTmp)) = speciesIdTmp
      else
        call error('Species in acsf output file do not match with dataset.')
      end if
    end do

    read(fd, *) isZscore

    if (isZscore == 0) then
      this%tZscore = .false.
    else
      this%tZscore = .true.
      allocate(this%zPrec(this%nRadial + this%nAngular, 2))
      do iAcsf = 1, this%nRadial + this%nAngular
        read(fd, *) this%zPrec(iAcsf, 1), this%zPrec(iAcsf, 2)
      end do
    end if

    close(fd)

  end subroutine TAcsf_fromFile


  subroutine TAcsfParams_init(this, nRadial, nAngular, rcut, kappa)

    !> representation of ACSF parameters
    type(TAcsfParams), intent(inout) :: this

    !> number of radial symmetry functions
    integer, intent(in) :: nRadial

    !> number of angular symmetry functions
    integer, intent(in) :: nAngular

    !> cutoff radius
    real(dp), intent(in) :: rcut

    !> optional kappa parameter for G3 function
    real(dp), intent(in), optional :: kappa

    !> stepsize for center parameter
    real(dp) :: rsStep

    !> auxiliary variable
    real(dp) :: xi

    !> auxiliary variables
    integer :: ii, jj, ind

    rsStep = rcut / (nRadial - 1)
    this%g2%eta = 5.0_dp * log(10.0_dp) / (2.0_dp * rsStep)**2

    allocate(this%g2%rs(nRadial))
    this%g2%rs(1) = 0.0_dp

    if (nRadial > 1) then
      this%g2%rs(2:) = [(ii * rsStep, ii = 1, nRadial - 1)]
    end if

    if (present(kappa)) then
      this%g3%kappa = kappa
    else
      this%g3%kappa = 1.0_dp
    end if

    this%g4%eta = 2.0_dp * log(10.0_dp) / (rcut)**2

    allocate(this%g4%lambda(nAngular))
    allocate(this%g4%xi(nAngular))

    ind = 1

    do ii = 0, ceiling(nAngular / 2.0_dp - 1.0_dp)
      if (nAngular <= 2) then
        xi = 1.0_dp
      else
        xi = 1.0_dp + ii * 30.0_dp / (nAngular - 2.0_dp)
      end if
      do jj = 1, -1, -2
        this%g4%lambda(ind) = jj
        this%g4%xi(ind) = xi
        if (ind >= nAngular) exit
        ind = ind + 1
      end do
    end do

  end subroutine TAcsfParams_init


  subroutine TMultiAcsfVals_init(this, geos, nAcsfVals)

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
      allocate(this%zPrec(this%nRadial + this%nAngular, 2))
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


  subroutine TAcsf_applyZscore(this)

    !> representation of ACSF mappings
    class(TAcsf), intent(inout) :: this

    !> auxiliary variables
    integer :: iGeo, iAtom, iAcsf

    if (.not. allocated(this%zPrec)) then
      call error('Cannot apply z-score standardization without means and variances.')   
    end if

    ! apply z-score standardization to ACSF mappings
    do iGeo = 1, size(this%vals%vals)
      do iAtom = 1, size(this%vals%vals(iGeo)%array, dim=2)
        do iAcsf = 1, size(this%vals%vals(iGeo)%array, dim=1)
          if (this%zPrec(iAcsf, 2) < 1e-08_dp) then
            cycle
          end if
          this%vals%vals(iGeo)%array(iAcsf, iAtom) = (this%vals%vals(iGeo)%array(iAcsf, iAtom) -&
              & this%zPrec(iAcsf, 1)) / this%zPrec(iAcsf, 2)
        end do
      end do
    end do

  end subroutine TAcsf_applyZscore


  subroutine TAcsf_calculate(this, geos, env, localAtToGlobalSp, atomIds, weights, zPrec)

    !> representation of ACSF mappings
    class(TAcsf), intent(inout) :: this

    !> system geometry container
    type(TGeometry), intent(in) :: geos(:)

    !> if compiled with mpi enabled, contains mpi communicator
    type(TEnv), intent(in) :: env

    !> index mapping local atom --> global species index
    type(TIntArray1D), intent(in) :: localAtToGlobalSp(:)

    !> optional atom dependent scaling parameters for cutoff function
    type(TRealArray1D), intent(in), optional :: atomIds(:)

    !> optional weighting of each corresponding datapoint
    integer, intent(in), optional :: weights(:)

    !> storage container of means and variances to calculate z-score
    real(dp), intent(in), optional :: zPrec(:,:)

    !> instance of dynamic neighbour list
    type(TDynNeighList), target :: neighList

    !> pointer to dynamic neighbour list
    type(TDynNeighList), pointer :: pNeighList

    !> temporary storage of ACSF values of each node
    type(TMultiAcsfVals) :: tmpVals

    !> atom dependent scaling parameters for cutoff function
    type(TRealArray1D), allocatable :: atomIdentifier(:)

    !> weighting of each corresponding datapoint
    integer, allocatable :: weighting(:)

    !> true, if current process is the lead
    logical :: tLead

    !> auxiliary variables
    integer :: iSys, iStart, iEnd

    allocate(atomIdentifier(size(geos)))
    allocate(weighting(size(geos)))

    if (present(atomIds)) then
      atomIdentifier(:) = atomIds
    else
      do iSys = 1, size(geos)
        allocate(atomIdentifier(iSys)%array(geos(iSys)%nAtom))
        atomIdentifier(iSys)%array(:) = 1.0_dp
      end do
    end if

    if (present(weights)) then
      weighting(:) = weights
    else
      weighting(:) = 1
    end if

    if (present(zPrec)) then
      this%zPrec = zPrec
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

    call TMultiAcsfVals_init(tmpVals, geos, this%nRadial + this%nAngular)
    call TMultiAcsfVals_init(this%vals, geos, this%nRadial + this%nAngular)

    lpSystem: do iSys = iStart, iEnd

      call TDynNeighList_init(neighList, this%rcut, geos(iSys)%nAtom, geos(iSys)%tPeriodic)
      call neighList%updateCoords(geos(iSys)%coords, localAtToGlobalSp(iSys)%array)
      if (geos(iSys)%tPeriodic) then
        call neighList%updateLatVecs(geos(iSys)%latVecs, geos(iSys)%recVecs2p)
      end if
      pNeighList => neighList
      call iAtomAcsf(tmpVals%vals(iSys), geos(iSys), pNeighList, this%param, this%nRadial,&
          & this%nAngular, this%rcut, this%speciesIds, atomIdentifier(iSys)%array)

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


  subroutine iAtomAcsf(this, geo, pNeighList, param, nRadial, nAngular, rcut, speciesIds, atomIds)

    !> symmetry function value instance of geometry
    type(TRealArray2D), intent(inout) :: this

    !> system geometry container
    type(TGeometry), intent(in) :: geo

    !> Pointer to dynamic neighbour list
    type(TDynNeighList), intent(in), pointer :: pNeighList

    !> container for all ACSF function parameters
    type(TAcsfParams), intent(in) :: param

    !> number of radial symmetry functions
    integer, intent(in) :: nRadial

    !> number of angular symmetry functions
    integer, intent(in) :: nAngular

    !> cutoff radius
    real(dp), intent(in) :: rcut

    !> species dependent scaling parameters for G1 function
    real(dp), intent(in) :: speciesIds(:)

    !> atom dependent scaling parameters for cutoff function
    real(dp), intent(in) :: atomIds(:)

    !> instance of dynamic neighbour list iterator
    type(TNeighIterator) :: neighIter

    !> obtained neighbor coordinates and distances
    real(dp) :: neighDists(iterChunkSize), neighCoords(3, iterChunkSize)

    !> obtained neighbor species indices
    integer :: speciesIndices(iterChunkSize)

    !> obtained neighbor atom indices
    integer :: atomIndices(iterChunkSize)

    !> neighbor coordinates, distances and species/atom ID's of all iterator stepsizes
    real(dp), allocatable :: allNeighDists(:), allNeighCoords(:,:), tmpCoords(:,:)
    real(dp), allocatable :: allSpeciesIds(:), allAtomIds(:)

    !> neighbor species/atom indices of all iterator stepsizes
    integer, allocatable :: allSpeciesIndices(:), allAtomIndices(:)

    !> auxiliary variables
    integer :: iAtom, iAcsf, iSp, iAtomInd, nNeigh, nTotNeigh

    do iAtom = 1, geo%nAtom

      allocate(tmpCoords(3, 0))
      allocate(allNeighDists(0))
      allocate(allSpeciesIndices(0))
      allocate(allAtomIndices(0))

      call TNeighIterator_init(neighIter, pNeighList, iAtom, includeSelf=.false.)

      nNeigh = iterChunkSize

      do while (nNeigh == iterChunkSize)
        call neighIter%getNextNeighbours(nNeigh, coords=neighCoords, dists=neighDists,&
            & speciesIndices=speciesIndices, atomIndices=atomIndices)
        nTotNeigh = size(tmpCoords, dim=2) + nNeigh
        allocate(allNeighCoords(3, nTotNeigh))
        allNeighDists = [allNeighDists, neighDists(1:nNeigh)]
        allSpeciesIndices = [allSpeciesIndices, speciesIndices(1:nNeigh)]
        allAtomIndices = [allAtomIndices, atomIndices(1:nNeigh)]
        allNeighCoords(:, 1:nTotNeigh - nNeigh) = tmpCoords
        allNeighCoords(:, nTotNeigh - nNeigh + 1:nTotNeigh) = neighCoords(:, 1:nNeigh)
        call move_alloc(allNeighCoords, tmpCoords)
      end do

      allNeighCoords = tmpCoords

      allocate(allSpeciesIds(size(allSpeciesIndices)))
      allocate(allAtomIds(size(allAtomIndices)))

      do iSp = 1, size(allSpeciesIds)
        allSpeciesIds(iSp) = speciesIds(allSpeciesIndices(iSp))
      end do

      do iAtomInd = 1, size(allAtomIds)
        allAtomIds(iAtomInd) = atomIds(allAtomIndices(iAtomInd))
      end do

      do iAcsf = 1, nRadial
        this%array(iAcsf, iAtom) = g2(allNeighDists, allSpeciesIds, allAtomIds, param%g2%eta,&
            & param%g2%rs(iAcsf), rcut)
      end do

      do iAcsf = 1 + nRadial, nRadial + nAngular
        this%array(iAcsf, iAtom) = g4(geo%coords(:, iAtom), allNeighCoords, allNeighDists,&
            & allSpeciesIds, allAtomIds, param%g4%xi(iAcsf - nRadial),&
            & param%g4%lambda(iAcsf - nRadial), param%g4%eta, rcut)
      end do

      deallocate(tmpCoords)
      deallocate(allNeighDists)
      deallocate(allNeighCoords)
      deallocate(allSpeciesIds)
      deallocate(allAtomIds)
      deallocate(allSpeciesIndices)
      deallocate(allAtomIndices)

    end do

  end subroutine iAtomAcsf


  pure function distance(geo, iAt1, iAt2)

    !> system geometry container
    type(TGeometry), intent(in) :: geo

    !> Indices of atoms to calculate distance for
    integer, intent(in) :: iAt1, iAt2

    !> Obtained distance betweeen the two atoms
    real(dp) :: distance

    distance = norm2(geo%coords(:, iAt2) - geo%coords(:, iAt1))

  end function distance


  pure function theta(atomCoords, neighCoords1, neighCoords2, neighDist1, neighDist2)

    !> coordinates of reference atom
    real(dp), intent(in) :: atomCoords(:)

    !> coordinates of two neighboring atoms
    real(dp), intent(in) :: neighCoords1(:), neighCoords2(:)

    !> distances of reference atom to first and second neighbor
    real(dp), intent(in) :: neighDist1, neighDist2

    !> obtained angle betweeen the three atoms
    real(dp) :: theta

    theta = acos(dot_product((neighCoords1(:) - atomCoords(:)), (neighCoords2(:) - atomCoords(:)))&
        & / (neighDist1 * neighDist2 + 1e-13_dp))

  end function theta


  pure function cutoff(rr, speciesId, atomId, rcut) result(res)

    !> atom distance (in cutoff range)
    real(dp), intent(in) :: rr

    !> species ID of neighbor
    real(dp), intent(in) :: speciesId

    !> atom ID of neighbor
    real(dp), intent(in) :: atomId

    !> cutoff radius
    real(dp), intent(in) :: rcut

    !> corresponding activation function values
    real(dp) :: res

    if (rr > rcut) then
      res = 0.0_dp
    else
      res = 0.5_dp * speciesId * atomId * (cos(pi * rr / rcut) + 1.0_dp)
    end if

  end function cutoff


  pure function cutoff_vec(rr, speciesIds, atomIds, rcut) result(res)

    !> array of atom distances (in cutoff range)
    real(dp), intent(in) :: rr(:)

    !> species ID's of all neighbors
    real(dp), intent(in) :: speciesIds(:)

    !> atom ID's of all neighbors
    real(dp), intent(in) :: atomIds(:)

    !> cutoff radius
    real(dp), intent(in) :: rcut

    !> corresponding activation function values
    real(dp) :: res(size(rr))

    where (rr > rcut)
      res = 0.0_dp
    elsewhere
      res = 0.5_dp * speciesIds * atomIds * (cos(pi * rr / rcut) + 1.0_dp)
    end where

  end function cutoff_vec


  pure function g1(rr, speciesIds, atomIds, rcut)

    !> array of atom distances in cutoff range
    real(dp), intent(in) :: rr(:)

    !> species ID's of all neighbors
    real(dp), intent(in) :: speciesIds(:)

    !> atom ID's of all neighbors
    real(dp), intent(in) :: atomIds(:)

    !> cutoff radius
    real(dp), intent(in) :: rcut

    !> corresponding symmetry function values
    real(dp) :: g1

    g1 = sum(cutoff_vec(rr, speciesIds, atomIds, rcut))

  end function g1


  pure function g2(rr, speciesIds, atomIds, eta, rs, rcut)

    !> array of atom distances in cutoff range
    real(dp), intent(in) :: rr(:)

    !> species ID's of all neighbors
    real(dp), intent(in) :: speciesIds(:)

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

    g2 = sum(exp(- eta * (rr - rs)**2) * cutoff_vec(rr, speciesIds, atomIds, rcut))

  end function g2


  pure function g3(rr, speciesIds, atomIds, kappa, rcut)

    !> array of atom distances in cutoff range
    real(dp), intent(in) :: rr(:)

    !> species ID's of all neighbors
    real(dp), intent(in) :: speciesIds(:)

    !> atom ID's of all neighbors
    real(dp), intent(in) :: atomIds(:)

    !> period length of damped cosine functions
    real(dp), intent(in) :: kappa

    !> cutoff radius
    real(dp), intent(in) :: rcut

    !> corresponding symmetry function values
    real(dp) :: g3

    g3 = sum(cos(kappa * rr) * cutoff_vec(rr, speciesIds, atomIds, rcut))

  end function g3


  pure function g4(atomCoords, neighCoords, neighDists, speciesIds, atomIds, xi, lambda, eta, rcut)

    !> coordinates of reference atom
    real(dp), intent(in) :: atomCoords(:)

    !> coordinates of neighboring atoms
    real(dp), intent(in) :: neighCoords(:,:)

    !> distances of reference atom to neighboring atoms
    real(dp), intent(in) :: neighDists(:)

    !> species ID's of all neighbors
    real(dp), intent(in) :: speciesIds(:)

    !> atom ID's of all neighbors
    real(dp), intent(in) :: atomIds(:)

    !> parameter specifying angular resolution
    real(dp), intent(in) :: xi

    !> lambda parameter of angular part, recommended values: [+1, -1]
    real(dp), intent(in) :: lambda

    !> width of Gaussian part
    real(dp), intent(in) :: eta

    !> cutoff radius
    real(dp), intent(in) :: rcut

    !> corresponding symmetry function values
    real(dp) :: g4

    !> auxiliary variables
    integer :: jj, kk

    g4 = 0.0_dp

    do jj = 1, size(neighCoords, dim=2)
      do kk = 1, size(neighCoords, dim=2)
        g4 = g4 + (1.0_dp + lambda * cos(theta(atomCoords, neighCoords(:, jj), neighCoords(:, kk),&
            & neighDists(jj), neighDists(kk))))**xi * exp(- eta * (neighDists(jj)**2 +&
            & neighDists(kk)**2)) * cutoff(neighDists(jj), speciesIds(jj), atomIds(jj), rcut) *&
            & cutoff(neighDists(kk), speciesIds(kk), atomIds(kk), rcut)
      end do
    end do

    g4 = g4 * 2.0_dp**(1.0_dp - xi)

  end function g4

end module fnet_acsf
