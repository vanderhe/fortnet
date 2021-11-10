!--------------------------------------------------------------------------------------------------!
!  FORTNET: A Behler-Parrinello-Neural-Network Implementation                                      !
!  Copyright (C) 2020 - 2021  T. W. van der Heide                                                  !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Provides the functionality to calculate atomic forces.
module fnet_forces

  use dftbp_assert
  use dftbp_accuracy, only: dp
  use dftbp_typegeometry, only : TGeometry

  use fnet_nestedtypes, only : TEnv, TPredicts, TIntArray1D, TRealArray2D, TJacobians
  use fnet_bpnn, only : TBpnn
  use fnet_acsf, only : TAcsf

#:if WITH_MPI
  use fnet_mpifx
  use fnet_parallel, only : getStartAndEndIndex
#:endif

  implicit none

  private

  public :: TForces, forceAnalysis
  public :: TGeometriesForFiniteDiff, TGeometriesForFiniteDiff_init
  public :: TGeometriesForFiniteDiff_deserialize


  type :: TForces

    !> contains atomic forces of geometries, exp. shape: [3 * nGlobalTargets, nAtom]
    type(TRealArray2D), allocatable :: geos(:)

  end type TForces


  type :: T6Geometries

    !> contains a geometry for every coordinate shift
    type(TGeometry) :: coord(6)

  end type T6Geometries


  type :: TAtomGeometries

    !> contains a geometry for every coordinate
    type(T6Geometries), allocatable :: atom(:)

  end type TAtomGeometries


  type :: TGeometriesForFiniteDiff

    !> contains 3N geometries per original geometry with shifted atom coordinates
    type(TAtomGeometries), allocatable :: geos(:)

    !> index mapping local atom --> global species index
    type(TIntArray1D), allocatable :: localAtToGlobalSp(:)

    !> index mapping local atom --> atomic number
    type(TIntArray1D), allocatable :: localAtToAtNum(:)

  contains

    procedure :: serialize => TGeometriesForFiniteDiff_serialize

  end type TGeometriesForFiniteDiff


  interface forceAnalysis
    module procedure :: forceAnalysis_finiteDiff, forceAnalysis_analytical
  end interface forceAnalysis


contains

  !> Initialises a structure of geometries with shifted atomic coordinates.
  pure subroutine TGeometriesForFiniteDiff_init(this, geos, delta)

    !> representation of geometries with shifted atomic coordinates
    type(TGeometriesForFiniteDiff), intent(out) :: this

    !> original reference structure to rattle
    type(TGeometry), intent(in) :: geos(:)

    !> absolute shift of each atomic coordinate
    real(dp), intent(in) :: delta

    !> auxiliary variables
    integer :: iGeo, iAtom

    allocate(this%geos(size(geos)))

    do iGeo = 1, size(this%geos)
      allocate(this%geos(iGeo)%atom(geos(iGeo)%nAtom))
      do iAtom = 1, size(this%geos(iGeo)%atom)
        this%geos(iGeo)%atom(iAtom)%coord(:) = geos(iGeo)

        ! shift entries manually, according to the central difference
        this%geos(iGeo)%atom(iAtom)%coord(1)%coords(1, iAtom)&
            & = this%geos(iGeo)%atom(iAtom)%coord(1)%coords(1, iAtom) - delta
        this%geos(iGeo)%atom(iAtom)%coord(2)%coords(1, iAtom)&
            & = this%geos(iGeo)%atom(iAtom)%coord(2)%coords(1, iAtom) + delta

        this%geos(iGeo)%atom(iAtom)%coord(3)%coords(2, iAtom)&
            & = this%geos(iGeo)%atom(iAtom)%coord(3)%coords(2, iAtom) - delta
        this%geos(iGeo)%atom(iAtom)%coord(4)%coords(2, iAtom)&
            & = this%geos(iGeo)%atom(iAtom)%coord(4)%coords(2, iAtom) + delta

        this%geos(iGeo)%atom(iAtom)%coord(5)%coords(3, iAtom)&
            & = this%geos(iGeo)%atom(iAtom)%coord(5)%coords(3, iAtom) - delta
        this%geos(iGeo)%atom(iAtom)%coord(6)%coords(3, iAtom)&
            & = this%geos(iGeo)%atom(iAtom)%coord(6)%coords(3, iAtom) + delta
      end do
    end do    

  end subroutine TGeometriesForFiniteDiff_init


  !> Serializes a nested structure of geometries with shifted atomic coordinates.
  pure subroutine TGeometriesForFiniteDiff_serialize(this, geos)

    !> representation of geometries with shifted atomic coordinates
    class(TGeometriesForFiniteDiff), intent(in) :: this

    !> serialized geometries
    type(TGeometry), intent(out), allocatable :: geos(:)

    !> auxiliary variables
    integer :: iGeo, iAtom, iCoord, nTotGeometries, ind

    nTotGeometries = 0
    ind = 1

    do iGeo = 1, size(this%geos)
      nTotGeometries = nTotGeometries + 6 * size(this%geos(iGeo)%atom)
    end do

    allocate(geos(nTotGeometries))

    do iGeo = 1, size(this%geos)
      do iAtom = 1, size(this%geos(iGeo)%atom)
        do iCoord = 1, 6
          geos(ind) = this%geos(iGeo)%atom(iAtom)%coord(iCoord)
          ind = ind + 1
        end do
      end do
    end do

  end subroutine TGeometriesForFiniteDiff_serialize


  !> Deserializes a flat structure of geometries with shifted atomic coordinates.
  pure subroutine TGeometriesForFiniteDiff_deserialize(this, geos)

    !> representation of geometries with shifted atomic coordinates
    type(TGeometriesForFiniteDiff), intent(out) :: this

    !> flat reference structure to nest
    type(TGeometry), intent(in) :: geos(:)

    !> auxiliary variables
    integer :: iGeo, iAtom, iCoord, ind

    ind = 1
    allocate(this%geos(size(geos)))

    do iGeo = 1, size(this%geos)
      allocate(this%geos(iGeo)%atom(geos(iGeo)%nAtom))
      do iAtom = 1, size(this%geos(iGeo)%atom)
        do iCoord = 1, 6
          this%geos(iGeo)%atom(iAtom)%coord(iCoord) = geos(ind)
          ind = ind + 1
        end do
      end do
    end do

  end subroutine TGeometriesForFiniteDiff_deserialize


  !> Predicts points of energy hypersurface and calculates forces based on the central difference.
  function forceAnalysis_finiteDiff(bpnn, forcesGeos, forcesAcsf, env, localAtToGlobalSp, delta)&
      & result(resForces)

    !> representation of a Behler-Parrinello neural network
    type(TBpnn), intent(in) :: bpnn

    !> representation of geometries with shifted atomic coordinates
    type(TGeometriesForFiniteDiff), intent(in) :: forcesGeos

    !> representation of ACSF mappings
    type(TAcsf), intent(in) :: forcesAcsf

    !> if compiled with mpi enabled, contains mpi communicator
    type(TEnv), intent(in) :: env

    !> index mapping local atom --> global species index
    type(TIntArray1D), intent(in) :: localAtToGlobalSp(:)

    !> absolute shift of each atomic coordinate
    real(dp), intent(in) :: delta

    !> extended index mapping local atom --> global species index
    type(TIntArray1D), allocatable :: localAtToGlobalSpExtended(:)

    !> network predictions for rattled geometries with one perturbed atom
    type(TPredicts) :: predicts, resPredicts

    !> obtained atomic forces
    type(TForces) :: forces, resForces

    !> true, if current process is the lead
    logical :: tLead

    !> auxiliary variables
    integer :: iSys, iGeo, iAtom, iCoord, iStart, iEnd, ii, ind

  #:if WITH_MPI
    tLead = env%globalMpiComm%lead
    call getStartAndEndIndex(size(forcesAcsf%vals%vals), env%globalMpiComm%size,&
        & env%globalMpiComm%rank, iStart, iEnd)
  #:else
    tLead = .true.
    iStart = 1
    iEnd = size(forcesAcsf%vals%vals)
  #:endif

    ! prepare structure that stores the predicted hypersurface of rattled geometries
    allocate(predicts%sys(size(forcesAcsf%vals%vals)))
    do iSys = 1, size(forcesAcsf%vals%vals)
      allocate(predicts%sys(iSys)%array(bpnn%dims(size(bpnn%dims)), 1))
      predicts%sys(iSys)%array(:,:) = 0.0_dp
    end do
    resPredicts = predicts

    ! expand index mapping local atom --> global species index to shifted geometries
    ind = 1
    allocate(localAtToGlobalSpExtended(size(forcesAcsf%vals%vals)))
    do iGeo = 1, size(forcesGeos%geos)
      do ii = 1, 6 * size(forcesGeos%geos(iGeo)%atom)
        localAtToGlobalSpExtended(ind+ii-1)%array = localAtToGlobalSp(iGeo)%array
      end do
      ind = ind + 6 * size(forcesGeos%geos(iGeo)%atom)
    end do

    do iSys = iStart, iEnd
      predicts%sys(iSys)%array(:,:) = bpnn%iPredict(forcesAcsf%vals%vals(iSys)%array,&
          & localAtToGlobalSpExtended(iSys)%array, .false.)
    end do

  #:if WITH_MPI
    do iSys = 1, size(forcesAcsf%vals%vals)
      call mpifx_allreduce(env%globalMpiComm, predicts%sys(iSys)%array,&
          & resPredicts%sys(iSys)%array, MPI_SUM)
    end do
    ! wait for all the predictions to finish
    call mpifx_barrier(env%globalMpiComm)
  #:else
    resPredicts = predicts
  #:endif

    ! prepare structure that stores the atomic forces
    allocate(forces%geos(size(forcesGeos%geos)))
    do iGeo = 1, size(forcesGeos%geos)
      allocate(forces%geos(iGeo)%array(3 * bpnn%dims(size(bpnn%dims)),&
          & forcesGeos%geos(iGeo)%atom(1)%coord(1)%nAtom))
      forces%geos(iGeo)%array(:,:) = 0.0_dp
    end do
    resForces = forces

  #:if WITH_MPI
    call getStartAndEndIndex(size(forcesGeos%geos), env%globalMpiComm%size,&
        & env%globalMpiComm%rank, iStart, iEnd)
  #:else
    iStart = 1
    iEnd = size(forcesGeos%geos)
  #:endif

    ! skip entries of previous MPI processes
    ind = 1
    do iGeo = 1, iStart - 1
      do iAtom = 1, size(forcesGeos%geos(iGeo)%atom)
        do iCoord = 1, 3
          ind = ind + 2
        end do
      end do
    end do
    do iGeo = iStart, iEnd
      do iAtom = 1, size(forcesGeos%geos(iGeo)%atom)
        do iCoord = 1, 3
          ! calculate negative gradients, i.e. forces based on the central finite difference
          forces%geos(iGeo)%array(iCoord::3, iAtom) =&
              & (predicts%sys(ind)%array(:, 1) - predicts%sys(ind+1)%array(:, 1)) / (2.0_dp * delta)
          ind = ind + 2
        end do
      end do
    end do

  #:if WITH_MPI
    do iGeo = 1, size(forcesGeos%geos)
      call mpifx_allreduce(env%globalMpiComm, forces%geos(iGeo)%array, resForces%geos(iGeo)%array,&
          & MPI_SUM)
    end do
  #:else
    resForces = forces
  #:endif

  end function forceAnalysis_finiteDiff


  !> Calculates forces based on the BPNN's Jacobian matrix and analytical ACSF derivatives.
  function forceAnalysis_analytical(bpnn, features, jacobian, forcesAcsf, env, localAtToGlobalSp)&
      & result(resForces)

    !> representation of a Behler-Parrinello neural network
    type(TBpnn), intent(in) :: bpnn

    !> atomic features as network input
    type(TRealArray2D), intent(in) :: features(:)

    !> contains Jacobians of multiple systems
    type(TJacobians), intent(in) :: jacobian

    !> contains ACSF derivatives, exp. shape: [forcesAcsf%vals%vals(iSys)%array(3, nAcsf, nAtom)]
    type(TAcsf), intent(in) :: forcesAcsf

    !> if compiled with mpi enabled, contains mpi communicator
    type(TEnv), intent(in) :: env

    !> index mapping local atom --> global species index
    type(TIntArray1D), intent(in) :: localAtToGlobalSp(:)

    !> network predictions for rattled geometries with one perturbed atom
    type(TPredicts) :: predicts, resPredicts

    !> obtained atomic forces
    type(TForces) :: forces, resForces

    !> true, if current process is the lead
    logical :: tLead

    !> auxiliary variables
    integer :: iSys, iAtom, iAcsf, iCoord, iStart, iEnd

  #:if WITH_MPI
    tLead = env%globalMpiComm%lead
    call getStartAndEndIndex(size(features), env%globalMpiComm%size, env%globalMpiComm%rank,&
        & iStart, iEnd)
  #:else
    tLead = .true.
    iStart = 1
    iEnd = size(features)
  #:endif

    ! prepare structure that stores the predicted hypersurface of rattled geometries
    allocate(predicts%sys(size(features)))
    do iSys = 1, size(features)
      allocate(predicts%sys(iSys)%array(bpnn%dims(size(bpnn%dims)), 1))
      predicts%sys(iSys)%array(:,:) = 0.0_dp
    end do
    resPredicts = predicts

    do iSys = iStart, iEnd
      predicts%sys(iSys)%array(:,:) = bpnn%iPredict(features(iSys)%array,&
          & localAtToGlobalSp(iSys)%array, .false.)
    end do

  #:if WITH_MPI
    do iSys = 1, size(features)
      call mpifx_allreduce(env%globalMpiComm, predicts%sys(iSys)%array,&
          & resPredicts%sys(iSys)%array, MPI_SUM)
    end do
    ! wait for all the predictions to finish
    call mpifx_barrier(env%globalMpiComm)
  #:else
    resPredicts = predicts
  #:endif

    ! prepare structure that stores the atomic forces
    allocate(forces%geos(size(features)))
    do iSys = 1, size(features)
      allocate(forces%geos(iSys)%array(3 * bpnn%dims(size(bpnn%dims)),&
          & size(features(iSys)%array, dim=2)))
      forces%geos(iSys)%array(:,:) = 0.0_dp
    end do
    resForces = forces

  #:if WITH_MPI
    call getStartAndEndIndex(size(features), env%globalMpiComm%size, env%globalMpiComm%rank,&
        & iStart, iEnd)
  #:else
    iStart = 1
    iEnd = size(features)
  #:endif

    do iSys = iStart, iEnd
      do iAtom = 1, size(features(iSys)%array, dim=2)
        do iAcsf = 1, size(features(iSys)%array, dim=1)
          do iCoord = 1, 3
            forces%geos(iSys)%array(iCoord::3, iAtom) = forces%geos(iSys)%array(iCoord::3, iAtom)&
                & - jacobian%sys(iSys)%atom(iAtom)%array(:, iAcsf)&
                & * forcesAcsf%valsPrime%vals(iSys)%array(iCoord, iAcsf, iAtom)
          end do
        end do
      end do
    end do

  #:if WITH_MPI
    do iSys = 1, size(features)
      call mpifx_allreduce(env%globalMpiComm, forces%geos(iSys)%array, resForces%geos(iSys)%array,&
          & MPI_SUM)
    end do
  #:else
    resForces = forces
  #:endif

  end function forceAnalysis_analytical

end module fnet_forces
