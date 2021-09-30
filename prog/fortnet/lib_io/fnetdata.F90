!--------------------------------------------------------------------------------------------------!
!  FORTNET: A Behler-Parrinello-Neural-Network Implementation                                      !
!  Copyright (C) 2020 - 2021  T. W. van der Heide                                                  !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Implements the generic dataset type as well as corresponding IO routines.
module fnet_fnetdata

  use h5lt
  use hdf5

  use dftbp_assert
  use dftbp_accuracy, only: dp, mc
  use dftbp_message, only : error
  use dftbp_typegeometry, only : TGeometry, normalize
  use dftbp_simplealgebra, only : invert33, determinant33
  use dftbp_charmanip, only : i2c
  use dftbp_constants, only : elementSymbol
  use dftbp_sorting, only : heap_sort
  use dftbp_globalenv, only : synchronizeAll

  use fnet_acsf, only : TAcsf
  use fnet_nestedtypes, only : TIntArray1D, TRealArray2D
  use fnet_hdf5fx, only : h5ltfx_read_dataset_int_f, h5ltfx_read_dataset_double_f
  use fnet_intmanip, only : getUniqueInt, getNumberOfUniqueInt

#:if WITH_MPI
  use fnet_mpifx
#:endif

  implicit none

  private

  public :: TDataset
  public :: readHdfDataset, checkDatasetCompatibility, checkBpnnDatasetCompatibility
  public :: checkAcsfDatasetCompatibility, checkExtFeaturesDatasetCompatibility
  public :: inquireStructures, inquireTargets, inquireExtFeatures

#:if WITH_MPI
  public :: syncDataset
#:endif


  !> Representation of a Fortnet compatible dataset.
  type TDataset

    !> number of datapoints in dataset
    integer :: nDatapoints

    !> number of dataset species
    integer :: nSpecies

    !> total number of atoms in the dataset structures
    integer :: nTotalAtoms

    !> true, if targets are provided
    logical :: tTargets

    !> true, if targets are atomic properties
    logical :: tAtomicTargets

    !> number of target values per atom (if tAtomicTargets = .true.) or system
    integer :: nTargets

    !> global atomic numbers present in the dataset
    integer, allocatable :: atomicNumbers(:)

    !> number of targets per network parameter
    real(dp) :: nTargetsPerParam

    !> target values for training or validation
    type(TRealArray2D), allocatable :: targets(:)

    !> contains obtained datapoint weights
    integer, allocatable :: weights(:)

    !> true, if geometries are provided
    logical :: tStructures

    !> contains obtained dataset geometries
    type(TGeometry), allocatable :: geos(:)

    !> index mapping local atom --> global species index
    type(TIntArray1D), allocatable :: localAtToGlobalSp(:)

    !> index mapping local atom --> atomic number
    type(TIntArray1D), allocatable :: localAtToAtNum(:)

    !> true, if external features are provided
    logical :: tExtFeatures

    !> total number of available external, atomic features in dataset
    integer :: nExtFeatures

    !> additional external, atomic features of the dataset
    type(TRealArray2D), allocatable :: extFeatures(:)

  contains

    procedure :: checkConsistency => TDataset_checkConsistency

  end type TDataset


contains

#:if WITH_MPI
  !> Synchronizes a dataset between the MPI nodes.
  subroutine syncDataset(dataset, comm)

    !> representation of a dataset
    type(TDataset), intent(inout) :: dataset

    !> mpi communicator with some additional information
    type(mpifx_comm), intent(in) :: comm

    !> auxiliary variables
    integer :: iData, dims0d
    integer, allocatable :: dims1d(:)

    ! synchronize scalar integers
    call mpifx_bcast(comm, dataset%nDatapoints)
    call mpifx_bcast(comm, dataset%nSpecies)
    call mpifx_bcast(comm, dataset%nTotalAtoms)
    call mpifx_bcast(comm, dataset%nTargets)
    call mpifx_bcast(comm, dataset%nTargetsPerParam)
    call mpifx_bcast(comm, dataset%nExtFeatures)

    ! synchronize scalar logicals
    call mpifx_bcast(comm, dataset%tTargets)
    call mpifx_bcast(comm, dataset%tAtomicTargets)
    call mpifx_bcast(comm, dataset%tStructures)
    call mpifx_bcast(comm, dataset%tExtFeatures)

    if (.not. comm%lead) then

      if (allocated(dataset%atomicNumbers)) deallocate(dataset%atomicNumbers)
      allocate(dataset%atomicNumbers(dataset%nSpecies))

      if (allocated(dataset%weights)) deallocate(dataset%weights)
      allocate(dataset%weights(dataset%nDatapoints))

    end if

    ! synchronize 1d integer arrays
    call mpifx_bcast(comm, dataset%atomicNumbers)
    call mpifx_bcast(comm, dataset%weights)

    if (.not. comm%lead) then

      if (dataset%tTargets) then
        if (allocated(dataset%targets)) deallocate(dataset%targets)
        allocate(dataset%targets(dataset%nDatapoints))
      end if

      if (dataset%tStructures) then
        if (allocated(dataset%geos)) deallocate(dataset%geos)
        allocate(dataset%geos(dataset%nDatapoints))
        if (allocated(dataset%localAtToGlobalSp)) deallocate(dataset%localAtToGlobalSp)
        allocate(dataset%localAtToGlobalSp(dataset%nDatapoints))
        if (allocated(dataset%localAtToAtNum)) deallocate(dataset%localAtToAtNum)
        allocate(dataset%localAtToAtNum(dataset%nDatapoints))
      end if

      if (dataset%tExtFeatures) then
        if (allocated(dataset%extFeatures)) deallocate(dataset%extFeatures)
        allocate(dataset%extFeatures(dataset%nDatapoints))
      end if

    end if

    ! synchronize derived types
    do iData = 1, dataset%nDatapoints
      if (dataset%tStructures) then
        call syncGeometry(comm, dataset%geos(iData))
        if (comm%lead) then
          dims0d = size(dataset%localAtToGlobalSp(iData)%array)
        end if
        call mpifx_bcast(comm, dims0d)
        if (.not. comm%lead) then
          if (allocated(dataset%localAtToGlobalSp(iData)%array)) then
            deallocate(dataset%localAtToGlobalSp(iData)%array)
          end if
          if (allocated(dataset%localAtToAtNum(iData)%array)) then
            deallocate(dataset%localAtToAtNum(iData)%array)
          end if
          allocate(dataset%localAtToGlobalSp(iData)%array(dims0d))
          allocate(dataset%localAtToAtNum(iData)%array(dims0d))
        end if
        call mpifx_bcast(comm, dataset%localAtToGlobalSp(iData)%array)
        call mpifx_bcast(comm, dataset%localAtToAtNum(iData)%array)
      end if
      if (dataset%tTargets) then
        if (comm%lead) then
          dims1d = shape(dataset%targets(iData)%array)
        else
          if (allocated(dims1d)) deallocate(dims1d)
          allocate(dims1d(2))
        end if
        call mpifx_bcast(comm, dims1d)
        if (.not. comm%lead) then
          if (allocated(dataset%targets(iData)%array)) then
            deallocate(dataset%targets(iData)%array)
          end if
          allocate(dataset%targets(iData)%array(dims1d(1), dims1d(2)))
        end if
        call mpifx_bcast(comm, dataset%targets(iData)%array)
      end if
      if (dataset%tExtFeatures) then
        if (comm%lead) then
          dims1d = shape(dataset%extFeatures(iData)%array)
        else
          if (allocated(dims1d)) deallocate(dims1d)
          allocate(dims1d(2))
        end if
        call mpifx_bcast(comm, dims1d)
        if (.not. comm%lead) then
          if (allocated(dataset%extFeatures(iData)%array)) then
            deallocate(dataset%extFeatures(iData)%array)
          end if
          allocate(dataset%extFeatures(iData)%array(dims1d(1), dims1d(2)))
        end if
        call mpifx_bcast(comm, dataset%extFeatures(iData)%array)
      end if
      ! ensure that the current loop iteration is completed
      call synchronizeAll()
    end do

  end subroutine syncDataset


  !> Synchronizes a geometry between the MPI nodes.
  subroutine syncGeometry(comm, geo)

    !> mpi communicator with some additional information
    type(mpifx_comm), intent(in) :: comm

    !> representation of a geometry
    type(TGeometry), intent(inout) :: geo

    ! synchronize scalar integers
    call mpifx_bcast(comm, geo%nAtom)
    call mpifx_bcast(comm, geo%nSpecies)

    ! synchronize scalar logicals
    call mpifx_bcast(comm, geo%tPeriodic)
    call mpifx_bcast(comm, geo%tFracCoord)
    call mpifx_bcast(comm, geo%tHelical)

    if (.not. comm%lead) then

      if (allocated(geo%species)) deallocate(geo%species)
      allocate(geo%species(geo%nAtom))

      if (allocated(geo%origin)) deallocate(geo%origin)
      allocate(geo%origin(3))

      if (allocated(geo%speciesNames)) deallocate(geo%speciesNames)
      allocate(geo%speciesNames(geo%nSpecies))

    end if

    ! synchronize 1d arrays
    call mpifx_bcast(comm, geo%species)
    if (geo%tPeriodic) then
      call mpifx_bcast(comm, geo%origin)
    end if
    call mpifx_bcast(comm, geo%speciesNames)

    if (.not. comm%lead) then

      if (allocated(geo%coords)) deallocate(geo%coords)
      allocate(geo%coords(3, geo%nAtom))

      if (geo%tPeriodic) then

        if (allocated(geo%latVecs)) deallocate(geo%latVecs)
        allocate(geo%latVecs(3, 3))

        if (allocated(geo%recVecs2p)) deallocate(geo%recVecs2p)
        allocate(geo%recVecs2p(3, 3))

      end if

    end if

    ! synchronize 2d arrays
    call mpifx_bcast(comm, geo%coords)
    if (geo%tPeriodic) then
      call mpifx_bcast(comm, geo%latVecs)
      call mpifx_bcast(comm, geo%recVecs2p)
    end if

  end subroutine syncGeometry
#:endif


  !> Inquires whether geometries are provided by the dataset file.
  subroutine inquireStructures(fname, tStructures)

    !> filename of the dataset
    character(len=*), intent(in) :: fname

    !> true, if the dataset provides structural information
    logical, intent(out) :: tStructures

    !> file identification
    integer(hid_t) :: file_id

    !> temporary storage container
    integer :: tmp(1)

    !> auxiliary variable
    integer :: iErr

    ! open the hdf5 interface
    call h5open_f(iErr)

    ! open the dataset file
    call h5fopen_f(fname, H5F_ACC_RDONLY_F, file_id, iErr)

    call h5ltget_attribute_int_f(file_id, 'fnetdata/dataset', 'withstructures', tmp, iErr)
    if (tmp(1) == 1) then
      tStructures = .true.
    else
      tStructures = .false.
    end if

    ! close the dataset file
    call h5fclose_f(file_id, iErr)

    ! close the hdf5 interface
    call h5close_f(iErr)

  end subroutine inquireStructures


  !> Inquires whether targets are provided by the dataset file.
  subroutine inquireTargets(fname, tTargets, tAtomic, nTargets)

    !> filename of the dataset
    character(len=*), intent(in) :: fname

    !> true, if the dataset provides target information
    logical, intent(out) :: tTargets

    !> optional information regarding the target kind
    logical, intent(out), optional :: tAtomic

    !> optional output for the number of targets
    integer, intent(out), optional :: nTargets

    !> file identification
    integer(hid_t) :: file_id

    !> temporary storage container
    integer :: tmp(1)

    !> auxiliary variables
    integer :: iErr

    ! open the hdf5 interface
    call h5open_f(iErr)

    ! open the dataset file
    call h5fopen_f(fname, H5F_ACC_RDONLY_F, file_id, iErr)

    call h5ltget_attribute_int_f(file_id, 'fnetdata/dataset/training', 'ntargets', tmp, iErr)

    if (tmp(1) > 0) then
      tTargets = .true.
    else
      tTargets = .false.
    end if

    if (present(nTargets)) then
      nTargets = tmp(1)
    end if

    call h5ltget_attribute_int_f(file_id, 'fnetdata/dataset/training', 'atomic', tmp, iErr)
    if (present(tAtomic)) then
      if (tmp(1) == 1) then
        tAtomic = .true.
      else
        tAtomic = .false.
      end if
    end if

    ! close the dataset file
    call h5fclose_f(file_id, iErr)

    ! close the hdf5 interface
    call h5close_f(iErr)

  end subroutine inquireTargets


  !> Inquires the number of external atomic features in fnetdata tree.
  subroutine inquireExtFeatures(fname, tFeatures, nFeatures)

    !> filename of the dataset
    character(len=*), intent(in) :: fname

    !> true, if the dataset provides external atomic input features
    logical, intent(out) :: tFeatures

    !> optional number of features per atom, provided by the dataset
    integer, intent(out), optional :: nFeatures

    !> file identification
    integer(hid_t) :: file_id

    !> temporary storage container
    integer :: tmp(1)

    !> auxiliary variables
    integer :: iErr

    ! open the hdf5 interface
    call h5open_f(iErr)

    ! open the dataset file
    call h5fopen_f(fname, H5F_ACC_RDONLY_F, file_id, iErr)

    ! read nextfeatures attribute
    call h5ltget_attribute_int_f(file_id, 'fnetdata/dataset', 'nextfeatures', tmp, iErr)
    if (tmp(1) > 0) then
      tFeatures = .true.
    else
      tFeatures = .false.
    end if

    if (present(nFeatures)) then
      nFeatures = tmp(1)
    end if

    ! close the dataset file
    call h5fclose_f(file_id, iErr)

    ! close the hdf5 interface
    call h5close_f(iErr)

  end subroutine inquireExtFeatures


  !> Checks the compatibility of two datasets.
  subroutine checkDatasetCompatibility(ref, comp)

    !> representation of a dataset
    type(TDataset), intent(in) :: ref, comp

    !> logicals to determine if species exist in reference setchildvalue
    logical, allocatable :: tFound(:)

    !> auxiliary variable
    integer :: iSp1, iSp2

    if (.not. (ref%tTargets .eqv. comp%tTargets)) then
      call error('Incompatible training and validation dataset. Check targets.')
    end if

    if (.not. (ref%tAtomicTargets .eqv. comp%tAtomicTargets)) then
      call error('Incompatible training and validation dataset. Check target kind.')
    end if

    if (.not. (ref%nTargets == comp%nTargets)) then
      call error('Incompatible training and validation dataset. Check number of targets.')
    end if

    if (.not. (ref%tStructures .eqv. comp%tStructures)) then
      call error('Incompatible training and validation dataset. Check structures.')
    end if

    if (.not. (ref%tExtFeatures .eqv. comp%tExtFeatures)) then
      call error('Incompatible training and validation dataset. Check external features.')
    end if

    if (.not. (ref%nExtFeatures == comp%nExtFeatures)) then
      call error('Incompatible training and validation dataset. Check number of external features.')
    end if

    if (ref%tStructures) then

      allocate(tFound(comp%nSpecies))
      tFound(:) = .false.

      ! check for global atomic numbers as species identifier
      outer: do iSp1 = 1, comp%nSpecies
        ! check if current element is in reference dataset
        inner: do iSp2 = 1, ref%nSpecies
          if (ref%atomicNumbers(iSp2) == comp%atomicNumbers(iSp1)) then
            tFound(iSp1) = .true.
            exit inner
          end if
        end do inner
      end do outer

      ! evaluate results of global species names
      if (any(.not. tFound)) then
        call error('Incompatible training and validation dataset. Check species.')
      end if

    end if

  end subroutine checkDatasetCompatibility


  !> Check dataset consistency by evaluating several assertions.
  subroutine TDataset_checkConsistency(this)

    !> representation of a dataset
    class(TDataset), intent(in) :: this

    @:ASSERT(this%nDatapoints >= 1)
    @:ASSERT(this%nSpecies >= 1)
    @:ASSERT(this%nExtFeatures >= 0)

    @:ASSERT(minval(this%weights) >= 1)

    if (this%tTargets) then
      @:ASSERT(this%nTargets > 0)
    else
      @:ASSERT(this%nTargets == 0)
    end if

    if (this%tExtFeatures) then
      @:ASSERT(this%nExtFeatures > 0)
    else
      @:ASSERT(this%nExtFeatures == 0)
    end if

    if (this%tStructures) then
      @:ASSERT(size(this%geos) >= 1)
      @:ASSERT(size(this%weights) == size(this%geos))
      @:ASSERT(size(this%atomicNumbers) == this%nSpecies)
      @:ASSERT(size(this%localAtToGlobalSp) == size(this%geos))
      @:ASSERT(size(this%localAtToAtNum) == size(this%geos))
    else
      @:ASSERT(this%nTargets == 0)
    end if

  end subroutine TDataset_checkConsistency


  !> Checks the compatibility of a dataset with an external features configuration.
  subroutine checkAcsfDatasetCompatibility(dataset, acsf)

    !> representation of a dataset
    type(TDataset), intent(in) :: dataset

    !> representation of acsf mapping information
    type(TAcsf), intent(in) :: acsf

    !> true, if the ACSF configuration is fully species-unresolved
    logical :: tReduce

    !> temporary storage container
    integer, allocatable :: datasetAtomicNumbers(:), acsfAtomicNumbers(:), tmpAtomicNumbers(:)

    !> auxiliary variable
    integer :: iFunc

    if (.not. dataset%tStructures) then
      call error('Error while parsing ACSF from netstat file. The dataset does not provide '&
          & // ' any structural information.')
    end if

    ! determine if the ACSF configuration is species-resolved
    allocate(tmpAtomicNumbers(0))
    do iFunc = 1, size(acsf%gFunctions%func)
      if (acsf%gFunctions%func(iFunc)%tRadial) then
        tmpAtomicNumbers = [tmpAtomicNumbers, [acsf%gFunctions%func(iFunc)%atomicNumbers(1)]]
      else
        tmpAtomicNumbers = [tmpAtomicNumbers, acsf%gFunctions%func(iFunc)%atomicNumbers]
      end if
    end do

    ! reduce the ACSF atomic numbers to unique ones
    call getUniqueInt(tmpAtomicNumbers, acsfAtomicNumbers)

    if (all(acsfAtomicNumbers == 0)) then
      tReduce = .true.
    else
      tReduce = .false.
    end if

    if (.not. tReduce) then
      ! compare their sorted equivalents
      datasetAtomicNumbers = dataset%atomicNumbers
      call heap_sort(datasetAtomicNumbers)
      call heap_sort(acsfAtomicNumbers)
      if (.not. all(shape(datasetAtomicNumbers) == shape(acsfAtomicNumbers))) then
        call error('Incompatibility in atomic number shape of dataset and ACSF.')
      end if
      if (.not. all(datasetAtomicNumbers == acsfAtomicNumbers)) then
        call error('Incompatibility in atomic numbers of dataset and ACSF.')
      end if
    end if

  end subroutine checkAcsfDatasetCompatibility


  !> Checks the compatibility of a dataset with an external features configuration.
  subroutine checkExtFeaturesDatasetCompatibility(dataset, indices)

    !> representation of a dataset
    type(TDataset), intent(in) :: dataset

    !> dataset indices of requested external features
    integer, intent(in) :: indices(:)

    if (any(indices <= 0)) then
      call error('Error while parsing external feature indices from netstat file. Zero or '&
          & // 'negative values obtained.')
    end if

    if (any(indices > dataset%nExtFeatures)) then
      call error('Error while parsing external feature indices from netstat file. One or more '&
          & // 'indices exceed the range of provided features of the dataset.')
    end if

  end subroutine checkExtFeaturesDatasetCompatibility


  !> Checks the compatibility of a dataset with the species of a BPNN.
  subroutine checkBpnnDatasetCompatibility(dataset, atomicNumbers, tAtomicTargets, allowSpSubset)

    !> representation of a dataset
    type(TDataset), intent(in) :: dataset

    !> atomic numbers of a BPNN representation
    integer, intent(in) :: atomicNumbers(:)

    !> true, if the BPNN was created for atomic targets
    logical, intent(in) :: tAtomicTargets

    !> true, if the dataset is allowed to only hold a subset of BPNN species
    logical, intent(in), optional :: allowSpSubset

    !> logicals to determine if species exist in reference setchildvalue
    logical, allocatable :: tFound(:)

    !> if present, equals optional dummy argument, otherwise false
    logical :: tAllowSpSubset

    !> temporary storage container
    integer, allocatable :: datasetAtomicNumbers(:), bpnnAtomicNumbers(:)

    !> auxiliary variables
    integer :: iAtNum1, iAtNum2

    if (present(allowSpSubset)) then
      tAllowSpSubset = allowSpSubset
    else
      tAllowSpSubset = .false.
    end if

    if (.not. (dataset%tAtomicTargets .eqv. tAtomicTargets)) then
      call error('Incompatibility in target kind of dataset and BPNN.')
    end if

    ! compare atomic numbers
    if (tAllowSpSubset) then
      allocate(tFound(size(dataset%atomicNumbers)))
      tFound(:) = .false.
      outer: do iAtNum1 = 1, size(dataset%atomicNumbers)
        ! check if current species name is in reference dataset
        inner: do iAtNum2 = 1, size(atomicNumbers)
          if (atomicNumbers(iAtNum2) == dataset%atomicNumbers(iAtNum1)) then
            tFound(iAtNum1) = .true.
            exit inner
          end if
        end do inner
      end do outer
      ! evaluate results of atomic numbers
      if (any(.not. tFound)) then
        call error('Incompatibility in atomic numbers of dataset and BPNN.')
      end if
    else
      datasetAtomicNumbers = dataset%atomicNumbers
      call heap_sort(datasetAtomicNumbers)
      bpnnAtomicNumbers = atomicNumbers
      call heap_sort(bpnnAtomicNumbers)
      if (.not. all(shape(datasetAtomicNumbers) == shape(bpnnAtomicNumbers))) then
        call error('Incompatibility in atomic number shape of dataset and BPNN.')
      end if
      if (.not. all(datasetAtomicNumbers == bpnnAtomicNumbers)) then
        call error('Incompatibility in atomic numbers of dataset and BPNN.')
      end if
    end if

  end subroutine checkBpnnDatasetCompatibility


  !> Reads a binary HDF5 dataset file.
  subroutine readHdfDataset(fname, dataset)

    !> filename of the dataset
    character(len=*), intent(in) :: fname

    !> representation of a dataset
    type(TDataset), intent(out) :: dataset

    !> file identification
    integer(hid_t) :: file_id, dataset_grp, training_grp, datapoint_grp, geometry_grp

    !> temporary storage container
    integer :: tmp(1)

    !> determinant to calculate inverse lattice vectors
    real(dp) :: det

    !> name of current datapoint
    character(len=:), allocatable :: dname

    !> auxiliary variable
    integer :: iErr, iDatapoint

    ! open the hdf5 interface
    call h5open_f(iErr)

    ! open the dataset file
    call h5fopen_f(fname, H5F_ACC_RDONLY_F, file_id, iErr)

    ! open the dataset group
    call h5gopen_f(file_id, 'fnetdata/dataset', dataset_grp, iErr)

    ! open the training group
    call h5gopen_f(dataset_grp, 'training', training_grp, iErr)

    call h5ltget_attribute_int_f(dataset_grp, './', 'ndatapoints', tmp, iErr)
    dataset%nDatapoints = tmp(1)

    call h5ltget_attribute_int_f(dataset_grp, './', 'nextfeatures', tmp, iErr)
    dataset%nExtFeatures = tmp(1)
    if (dataset%nExtFeatures > 0) then
      dataset%tExtFeatures = .true.
      allocate(dataset%extFeatures(dataset%nDatapoints))
    else
      dataset%tExtFeatures = .false.
    end if

    call h5ltget_attribute_int_f(dataset_grp, './', 'ntotatoms', tmp, iErr)
    dataset%nTotalAtoms = tmp(1)

    call h5ltget_attribute_int_f(dataset_grp, './', 'withstructures', tmp, iErr)
    if (tmp(1) == 1) then
      dataset%tStructures = .true.
      allocate(dataset%geos(dataset%nDatapoints))
      allocate(dataset%localAtToAtNum(dataset%nDatapoints))
      allocate(dataset%localAtToGlobalSp(dataset%nDatapoints))
    else
      dataset%tStructures = .false.
    end if

    call h5ltget_attribute_int_f(training_grp, './', 'ntargets', tmp, iErr)
    dataset%nTargets = tmp(1)
    if (dataset%nTargets > 0) then
      dataset%tTargets = .true.
      allocate(dataset%targets(dataset%nDatapoints))
    else
      dataset%tTargets = .false.
    end if

    call h5ltget_attribute_int_f(training_grp, './', 'atomic', tmp, iErr)
    if (tmp(1) == 1) then
      dataset%tAtomicTargets = .true.
    else
      dataset%tAtomicTargets = .false.
    end if

    call h5ltfx_read_dataset_int_f(dataset_grp, 'atomicnumbers', dataset%atomicNumbers)
    dataset%nSpecies = size(dataset%atomicNumbers)

    allocate(dataset%weights(dataset%nDatapoints))

    do iDatapoint = 1, dataset%nDatapoints
      dname = 'datapoint' // i2c(iDatapoint)

      ! open the datapoint group
      call h5gopen_f(dataset_grp, dname, datapoint_grp, iErr)

      ! get weight of datapoint
      call h5ltget_attribute_int_f(datapoint_grp, './', 'weight', tmp, iErr)
      dataset%weights(iDatapoint) = tmp(1)

      if (dataset%tStructures) then
        ! open the geometry group
        call h5gopen_f(datapoint_grp, 'geometry', geometry_grp, iErr)

        ! check for periodic boundary conditions
        call h5ltget_attribute_int_f(geometry_grp, './', 'periodic', tmp, iErr)
        if (tmp(1) == 1) then
          dataset%geos(iDatapoint)%tPeriodic = .true.
        else
          dataset%geos(iDatapoint)%tPeriodic = .false.
        end if

        ! no support for helical boundary conditions
        dataset%geos(iDatapoint)%tHelical = .false.

        ! check for fractional coordinates
        call h5ltget_attribute_int_f(geometry_grp, './', 'fractional', tmp, iErr)
        if (tmp(1) == 1) then
          dataset%geos(iDatapoint)%tFracCoord = .true.
        else
          dataset%geos(iDatapoint)%tFracCoord = .false.
        end if

        ! read mappings from local atom index to atomic number
        call h5ltfx_read_dataset_int_f(geometry_grp, 'localattoatnum',&
            & dataset%localAtToAtNum(iDatapoint)%array)
        call getNumberOfUniqueInt(dataset%localAtToAtNum(iDatapoint)%array,&
            & dataset%geos(iDatapoint)%nSpecies)

        ! read mappings from local atom index to local species index
        call h5ltfx_read_dataset_int_f(geometry_grp, 'localattolocalsp',&
            & dataset%geos(iDatapoint)%species)

        ! read mappings from local atom index to global species index
        call h5ltfx_read_dataset_int_f(geometry_grp, 'localattoglobalsp',&
            & dataset%localAtToGlobalSp(iDatapoint)%array)

        ! extract local species names from atomic numbers
        call getReducedSpeciesNames(dataset%geos(iDatapoint)%species,&
            & dataset%localAtToAtNum(iDatapoint)%array, dataset%geos(iDatapoint)%speciesNames)

        ! read coordinates
        call h5ltfx_read_dataset_double_f(geometry_grp, 'coordinates',&
            & dataset%geos(iDatapoint)%coords)
        dataset%geos(iDatapoint)%nAtom = size(dataset%geos(iDatapoint)%coords, dim=2)

        ! check validity of species
        if (any(dataset%geos(iDatapoint)%species < 1&
            & .or. dataset%geos(iDatapoint)%species > dataset%geos(iDatapoint)%nSpecies)) then
          call error('Error while parsing dataset species. Type index must be between 1 and '&
              & // i2c(dataset%geos(iDatapoint)%nSpecies) // '.')
        end if

        if (dataset%geos(iDatapoint)%tPeriodic) then
          call h5ltfx_read_dataset_double_f(geometry_grp, 'basis',&
              & dataset%geos(iDatapoint)%latVecs)
          ! set dummy origin
          dataset%geos(iDatapoint)%origin = [0.0_dp, 0.0_dp, 0.0_dp]
          dataset%geos(iDatapoint)%coords(:,:) = dataset%geos(iDatapoint)%coords&
              & - spread(dataset%geos(iDatapoint)%origin, 2, dataset%geos(iDatapoint)%nAtom)
          if (dataset%geos(iDatapoint)%tFracCoord) then
            dataset%geos(iDatapoint)%coords(:,:) = matmul(dataset%geos(iDatapoint)%latVecs,&
                & dataset%geos(iDatapoint)%coords)
            dataset%geos(iDatapoint)%origin(:) = matmul(dataset%geos(iDatapoint)%latVecs,&
                & dataset%geos(iDatapoint)%origin)
          end if
          allocate(dataset%geos(iDatapoint)%recVecs2p(3, 3))
          det = determinant33(dataset%geos(iDatapoint)%latVecs)
          if (abs(det) < 1e-12_dp) then
            call error('Dependent lattice vectors in ' // dname // '.')
          end if
          call invert33(dataset%geos(iDatapoint)%recVecs2p, dataset%geos(iDatapoint)%latVecs, det)
          dataset%geos(iDatapoint)%recVecs2p(:,:) =&
              & reshape(dataset%geos(iDatapoint)%recVecs2p, [3, 3], order=[2, 1])
        end if

        call normalize(dataset%geos(iDatapoint))

        ! close the geometry group
        call h5gclose_f(geometry_grp, iErr)
      end if

      if (dataset%tTargets) then
        call h5ltfx_read_dataset_double_f(datapoint_grp, 'targets',&
            & dataset%targets(iDatapoint)%array)
        dataset%targets(iDatapoint)%array = transpose(dataset%targets(iDatapoint)%array)
      end if

      if (dataset%tExtFeatures) then
        call h5ltfx_read_dataset_double_f(datapoint_grp, 'extfeatures',&
            & dataset%extFeatures(iDatapoint)%array)
      end if

      ! close the datapoint group
      call h5gclose_f(datapoint_grp, iErr)

    end do

    ! close the training group
    call h5gclose_f(training_grp, iErr)

    ! close the dataset group
    call h5gclose_f(dataset_grp, iErr)

    ! close the dataset file
    call h5fclose_f(file_id, iErr)

    ! close the hdf5 interface
    call h5close_f(iErr)

    ! perform some basic consistency checks
    call dataset%checkConsistency()

  end subroutine readHdfDataset


  !> Returns the local species names of a geometry in the local order.
  subroutine getReducedSpeciesNames(localAtToLocalSp, localAtToAtNum, speciesNames)

    !> index mapping local atom --> local species index
    integer, intent(in) :: localAtToLocalSp(:)

    !> index mapping local atom --> atomic number
    integer, intent(in) :: localAtToAtNum(:)

    !> unique species name in order of local species indices
    character(mc), intent(out), allocatable :: speciesNames(:)

    !> auxiliary variables
    integer :: ii, jj, kk
    integer :: tmp1(size(localAtToAtNum)), tmp2(size(localAtToLocalSp))
    integer, allocatable :: uniqueAtNum(:), uniqueLocalSp(:)

    kk = 1
    tmp1(1) = localAtToAtNum(1)
    tmp2(1) = localAtToLocalSp(1)

    outer: do ii = 2, size(localAtToAtNum)
      do jj = 1, kk
        if (tmp1(jj) == localAtToAtNum(ii)) then
          cycle outer
        end if
      end do
      kk = kk + 1
      tmp1(kk) = localAtToAtNum(ii)
      tmp2(kk) = localAtToLocalSp(ii)
    end do outer

    uniqueAtNum = tmp1(1:kk)
    uniqueLocalSp = tmp2(1:kk)

    allocate(speciesNames(size(uniqueAtNum)))

    ! map unique indices to species names
    do ii = 1, size(uniqueAtNum)
      speciesNames(uniqueLocalSp(ii)) = trim(elementSymbol(uniqueAtNum(ii)))
    end do

  end subroutine getReducedSpeciesNames

end module fnet_fnetdata
