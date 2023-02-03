!--------------------------------------------------------------------------------------------------!
!  FORTNET: A Behler-Parrinello-Neural-Network Implementation                                      !
!  Copyright (C) 2020 - 2023  T. W. van der Heide                                                  !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Defines a derived type as well as routines to handle atomic input features.
module fnet_features

  use dftbp_accuracy, only: dp
  use dftbp_message, only : error

#:if WITH_SOCKETS
  use dftbp_typegeometry, only : TGeometry
#:endif

  use fnet_acsf, only : TAcsf, TGFunctions
  use fnet_fnetdata, only : TDataset
  use fnet_nestedtypes, only : TRealArray2D

#:if WITH_MPI
  use fnet_mpifx
#:endif

  implicit none

  private

  public :: TMappingBlock, TExternalBlock, TFeaturesBlock
  public :: TFeatures, TFeatures_init, TFeatures_collect

#:if WITH_SOCKETS
  public :: TFeatures_initForSocketComm, TFeatures_collectForSocketComm
#:endif


  !> Data type containing variables of the Mapping block.
  type TMappingBlock

    !> true, if the species-resolved features should get summed up
    logical :: tReduce

    !> true, if standardization of features is desired
    logical :: tStandardize

    !> number of radial and angular mappings
    integer :: nRadial, nAngular

    !> wrapper around multiple G-functions
    type(TGFunctions) :: functions

    !> structural mapping type (currently only ACSF's are available)
    character(len=:), allocatable :: type

  end type TMappingBlock


  !> Data type containing variables of the External block.
  type TExternalBlock

    !> total number of requested external, atomic features
    integer :: nExtFeatures

    !> indices of external features of dataset to incorporate
    integer, allocatable :: indices(:)

#:if WITH_MPI
  contains
    procedure :: syncConfig => TExternalBlock_syncConfig
#:endif

  end type TExternalBlock


  !> Data type containing variables of the Features block.
  type TFeaturesBlock

    !> true, if mapping features are requested
    logical :: tMappingFeatures

    !> data type containing variables of the Mapping block
    type(TMappingBlock) :: mapping

    !> true, if external features are requested
    logical :: tExtFeatures

    !> data type containing variables of the External block
    type(TExternalBlock) :: ext

    !> total number of features, considering all inputs
    integer :: nFeatures

  end type TFeaturesBlock


  !> Data type containing the collected features.
  type TFeatures

    !> collection of input features for training
    type(TRealArray2D), allocatable :: trainFeatures(:)

    !> collection of input features for validation
    type(TRealArray2D), allocatable :: validFeatures(:)

#:if WITH_MPI
  contains
    procedure :: sync => TFeatures_sync
#:endif

  end type TFeatures


contains

  !> Initialises a feature instance.
  subroutine TFeatures_init(features, inpFeatures, trainDataset, validDataset, trainAcsf,&
      & tMonitorValid)

    !> collected features
    type(TFeatures), intent(inout) :: features

    !> input information of the features HSD block
    type(TFeaturesBlock), intent(in) :: inpFeatures

    !> representation of a training and validation dataset
    type(TDataset), intent(in) :: trainDataset, validDataset

    !> representation of acsf mapping informations
    type(TAcsf), intent(in) :: trainAcsf

    !> true, if validation monitoring is desired
    logical, intent(in) :: tMonitorValid

    !! total number of input features
    integer :: nFeatures

    !! auxiliary variable
    integer :: iData

    nFeatures = 0

    if (inpFeatures%tMappingFeatures) then
      nFeatures = nFeatures + size(trainAcsf%gFunctions%func)
    end if

    if (inpFeatures%tExtFeatures) then
      nFeatures = nFeatures + inpFeatures%ext%nExtFeatures
    end if

    if (nFeatures == 0) then
      call error('No features present to collect. Aborting.')
    elseif (nFeatures /= inpFeatures%nFeatures) then
      call error('Mismatch in number of features found. Aborting.')
    end if

    allocate(features%trainFeatures(trainDataset%nDatapoints))

    if (trainDataset%tStructures) then
      do iData = 1, trainDataset%nDatapoints
        allocate(features%trainFeatures(iData)%array(nFeatures, trainDataset%geos(iData)%nAtom))
        features%trainFeatures(iData)%array(:,:) = 0.0_dp
      end do
    elseif (trainDataset%tExtFeatures) then
      do iData = 1, trainDataset%nDatapoints
        allocate(features%trainFeatures(iData)%array(nFeatures,&
            & size(trainDataset%extFeatures(iData)%array, dim=2)))
        features%trainFeatures(iData)%array(:,:) = 0.0_dp
      end do
    else
      call error('Error while initializing features structure. Neither structures nor '&
          & // 'external features present.')
    end if

    if (tMonitorValid) then
      allocate(features%validFeatures(validDataset%nDatapoints))
      if (validDataset%tStructures) then
        do iData = 1, validDataset%nDatapoints
          allocate(features%validFeatures(iData)%array(nFeatures, validDataset%geos(iData)%nAtom))
          features%validFeatures(iData)%array(:,:) = 0.0_dp
        end do
      elseif (validDataset%tExtFeatures) then
        do iData = 1, validDataset%nDatapoints
          allocate(features%validFeatures(iData)%array(nFeatures,&
              & size(validDataset%extFeatures(iData)%array, dim=2)))
          features%validFeatures(iData)%array(:,:) = 0.0_dp
        end do
      else
        call error('Error while initializing features structure. Neither structures nor '&
            & // 'external features present.')
      end if
    end if

  end subroutine TFeatures_init


  !> Collects features from ACSF calculations and external sources.
  subroutine TFeatures_collect(features, inpFeatures, trainDataset, validDataset, trainAcsf,&
      & validAcsf, tMonitorValid)

    !> collected features of data and mapping block
    type(TFeatures), intent(inout) :: features

    !> input information of the features HSD block
    type(TFeaturesBlock), intent(in) :: inpFeatures

    !> representation of a training and validation dataset
    type(TDataset), intent(in) :: trainDataset, validDataset

    !> representation of acsf mapping informations
    type(TAcsf), intent(in) :: trainAcsf, validAcsf

    !> true, if validation monitoring is desired
    logical, intent(in) :: tMonitorValid

    !! auxiliary variables
    integer :: iData, nAcsf

    if (inpFeatures%tMappingFeatures) then
      nAcsf = size(trainAcsf%gFunctions%func)
    else
      nAcsf = 0
    end if

    if (inpFeatures%tMappingFeatures .and. inpFeatures%tExtFeatures) then
      do iData = 1, trainDataset%nDatapoints
        features%trainFeatures(iData)%array(1:nAcsf, :) = trainAcsf%vals%vals(iData)%array
        features%trainFeatures(iData)%array(nAcsf+1:inpFeatures%nFeatures, :)&
            & = trainDataset%extFeatures(iData)%array(inpFeatures%ext%indices, :)
      end do
    elseif (inpFeatures%tMappingFeatures .and. (.not. inpFeatures%tExtFeatures)) then
      do iData = 1, trainDataset%nDatapoints
        features%trainFeatures(iData)%array = trainAcsf%vals%vals(iData)%array
      end do
    elseif ((.not. inpFeatures%tMappingFeatures) .and. inpFeatures%tExtFeatures) then
      do iData = 1, trainDataset%nDatapoints
        features%trainFeatures(iData)%array&
            & = trainDataset%extFeatures(iData)%array(inpFeatures%ext%indices, :)
      end do
    else
      call error('No features present to collect. Aborting.')
    end if

    if (tMonitorValid) then
      if (inpFeatures%tMappingFeatures .and. inpFeatures%tExtFeatures) then
        do iData = 1, validDataset%nDatapoints
          features%validFeatures(iData)%array(1:nAcsf, :) = validAcsf%vals%vals(iData)%array
          features%validFeatures(iData)%array(nAcsf+1:inpFeatures%nFeatures, :)&
              & = validDataset%extFeatures(iData)%array(inpFeatures%ext%indices, :)
        end do
      elseif (inpFeatures%tMappingFeatures .and. (.not. inpFeatures%tExtFeatures)) then
        do iData = 1, validDataset%nDatapoints
          features%validFeatures(iData)%array = validAcsf%vals%vals(iData)%array
        end do
      elseif ((.not. inpFeatures%tMappingFeatures) .and. inpFeatures%tExtFeatures) then
        do iData = 1, validDataset%nDatapoints
          features%validFeatures(iData)%array&
              & = validDataset%extFeatures(iData)%array(inpFeatures%ext%indices, :)
        end do
      end if
    end if

  end subroutine TFeatures_collect


#:if WITH_MPI
  !> Synchronizes the collected features between the MPI nodes.
  subroutine TFeatures_sync(this, comm)

    !> collected features of data and mapping block
    class(TFeatures), intent(in) :: this

    !> mpi communicator with some additional information
    type(mpifx_comm), intent(in) :: comm

    !! auxiliary variable
    integer :: iData

    do iData = 1, size(this%trainFeatures)
      call mpifx_bcast(comm, this%trainFeatures(iData)%array)
    end do

    if (allocated(this%validFeatures)) then
      do iData = 1, size(this%validFeatures)
        call mpifx_bcast(comm, this%validFeatures(iData)%array)
      end do
    end if

  end subroutine TFeatures_sync


  !> Synchronizes the collected features between the MPI nodes.
  subroutine TExternalBlock_syncConfig(this, comm)

    !> external feature configuration to synchronize
    class(TExternalBlock), intent(inout) :: this

    !> mpi communicator with some additional information
    type(mpifx_comm), intent(in) :: comm

    call mpifx_bcast(comm, this%nExtFeatures)

    if (.not. comm%lead) then
      if (allocated(this%indices)) deallocate(this%indices)
      allocate(this%indices(this%nExtFeatures))
    end if
    call mpifx_bcast(comm, this%indices)

  end subroutine TExternalBlock_syncConfig
#:endif


#:if WITH_SOCKETS
  !> Initialises a feature instance.
  subroutine TFeatures_initForSocketComm(features, acsf, geo)

    !> collected features
    type(TFeatures), intent(inout) :: features

    !> representation of acsf mapping informations
    type(TAcsf), intent(in) :: acsf

    !> geometry instance
    type(TGeometry), intent(in) :: geo

    !! total number of input features
    integer :: nFeatures

    ! potentially dangerous, assumes that acsf are present/allocated
    nFeatures = size(acsf%gFunctions%func)

    if (nFeatures == 0) then
      call error('No features present to collect. Aborting.')
    end if

    ! we only want to calculate a single structure
    if (allocated(features%trainFeatures)) deallocate(features%trainFeatures)
    allocate(features%trainFeatures(1))
    allocate(features%trainFeatures(1)%array(nFeatures, geo%nAtom))
    features%trainFeatures(1)%array(:,:) = 0.0_dp

  end subroutine TFeatures_initForSocketComm


  !> Collects features from ACSF calculations and external sources.
  subroutine TFeatures_collectForSocketComm(features, acsf)

    !> collected features of data and mapping block
    type(TFeatures), intent(inout) :: features

    !> representation of acsf mapping informations
    type(TAcsf), intent(in) :: acsf

    !! number of ACSF mappings
    integer :: nAcsf

    nAcsf = size(acsf%gFunctions%func)

    if (nAcsf > 0) then
      features%trainFeatures(1)%array = acsf%vals%vals(1)%array
    else
      call error('No features present to collect. Aborting.')
    end if

  end subroutine TFeatures_collectForSocketComm
#:endif

end module fnet_features
