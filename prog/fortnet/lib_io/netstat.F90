!--------------------------------------------------------------------------------------------------!
!  FORTNET: A Behler-Parrinello-Neural-Network Implementation                                      !
!  Copyright (C) 2020 - 2024  T. W. van der Heide                                                  !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Implements the IO routines for processing the netstat files as Fortnet's configuration file.
module fnet_netstat

  use h5lt
  use hdf5

  use dftbp_message, only : error
  use dftbp_charmanip, only : i2c, tolower
  use dftbp_constants, only : elementSymbol

  use fnet_acsf, only : TAcsf
  use fnet_bpnn, only : TBpnn
  use fnet_hdf5fx, only : h5ltfx_read_dataset_int_f, h5ltfxmake_dataset_int_f
  use fnet_features, only : TExternalBlock

  implicit none

  private

  public :: readExtFeaturesConfig
  public :: createNetstat, readSubnetArchitecture
  public :: writeBpnnHeader, inquireAcsf
  public :: inquireExtFeatures, writeExtFeaturesConfig


contains

  !> Creates an empty netstat file with base group.
  subroutine createNetstat(fname)

    !> filename of netstat file
    character(len=*), intent(in) :: fname

    !! file and group identification
    integer(hid_t) :: file_id, netstat_id

    !! auxiliary variable
    integer :: iErr

    ! open the hdf5 interface
    call h5open_f(iErr)

    ! create the netstat file
    call h5fcreate_f(fname, H5F_ACC_TRUNC_F, file_id, iErr)

    ! create netstat group
    call h5gcreate_f(file_id, 'netstat', netstat_id, iErr)

    ! close netstat group
    call h5gclose_f(netstat_id, iErr)

    ! close the netstat file
    call h5fclose_f(file_id, iErr)

    ! close the hdf5 interface
    call h5close_f(iErr)

  end subroutine createNetstat


  !> Inquires whether netstat provides external feature information.
  subroutine inquireExtFeatures(fname, tExtFeatures)

    !> filename of netstat file
    character(len=*), intent(in) :: fname

    !> true, if ACSF mapping features are present
    logical, intent(out) :: tExtFeatures

    !! file and group identifier
    integer(hid_t) :: file_id, netstat_id, external_id

    !! auxiliary variables
    integer :: iErr, tExist, tmp(1)

    ! open the hdf5 interface
    call h5open_f(iErr)

    ! open the netstat file
    call h5fopen_f(fname, H5F_ACC_RDONLY_F, file_id, iErr)

    ! open the netstat group
    call h5gopen_f(file_id, 'netstat', netstat_id, iErr)

    ! check if the external group is present
    tExist = h5ltfind_dataset_f(netstat_id, 'external')

    if (tExist == 1) then
      ! open the external group
      call h5gopen_f(netstat_id, 'external', external_id, iErr)
      ! check the number of external features
      call h5ltget_attribute_int_f(external_id, './', 'nextfeatures', tmp, iErr)
      if (tmp(1) > 0) then
        tExtFeatures = .true.
      else
        tExtFeatures = .false.
      end if
      ! close the external group
      call h5gclose_f(external_id, iErr)
    else
      tExtFeatures = .false.
    end if

    ! close the netstat group
    call h5gclose_f(netstat_id, iErr)

    ! close the netstat file
    call h5fclose_f(file_id, iErr)

    ! close the hdf5 interface
    call h5close_f(iErr)

  end subroutine inquireExtFeatures


  !> Inquires whether netstat provides ACSF mapping information.
  subroutine inquireAcsf(fname, tAcsf)

    !> filename of netstat file
    character(len=*), intent(in) :: fname

    !> true, if ACSF mapping features are present
    logical, intent(out) :: tAcsf

    !! file and group identifier
    integer(hid_t) :: file_id, netstat_id, mapping_id

    !! auxiliary variables
    character(len=100) :: tmpStr
    integer :: iErr, tExist

    ! open the hdf5 interface
    call h5open_f(iErr)

    ! open the netstat file
    call h5fopen_f(fname, H5F_ACC_RDONLY_F, file_id, iErr)

    ! open the netstat group
    call h5gopen_f(file_id, 'netstat', netstat_id, iErr)

    ! check if a mapping entry is present
    tExist = h5ltfind_dataset_f(netstat_id, 'mapping')

    if (tExist == 1) then
      ! open the mapping group
      call h5gopen_f(netstat_id, 'mapping', mapping_id, iErr)
      ! check the mapping type
      call h5ltget_attribute_string_f(mapping_id, './', 'type', tmpStr, iErr)
      if (tolower(trim(tmpStr)) == 'acsf') then
        tAcsf = .true.
      else
        tAcsf = .false.
      end if
      ! close the mapping group
      call h5gclose_f(mapping_id, iErr)
    else
      tAcsf = .false.
    end if

    ! close the netstat group
    call h5gclose_f(netstat_id, iErr)

    ! close the netstat file
    call h5fclose_f(file_id, iErr)

    ! close the hdf5 interface
    call h5close_f(iErr)

  end subroutine inquireAcsf


  !> Reads external feature configuration from netstat file.
  subroutine readExtFeaturesConfig(fname, ext)

    !> filename of netstat file
    character(len=*), intent(in) :: fname

    !> variables of the External block read from a netstat file
    type(TExternalBlock), intent(out) :: ext

    !! file and group identifier
    integer(hid_t) :: file_id, netstat_id, external_id

    !! auxiliary variables
    integer :: iErr, tExist, tmp(1)

    ! open the hdf5 interface
    call h5open_f(iErr)

    ! open the netstat file
    call h5fopen_f(fname, H5F_ACC_RDONLY_F, file_id, iErr)

    ! open the netstat group
    call h5gopen_f(file_id, 'netstat', netstat_id, iErr)

    ! check if the external group is present
    tExist = h5ltfind_dataset_f(netstat_id, 'external')

    if (tExist == 1) then
      ! open the external group
      call h5gopen_f(netstat_id, 'external', external_id, iErr)
      ! check the number of external features
      call h5ltget_attribute_int_f(external_id, './', 'nextfeatures', tmp, iErr)
      if (tmp(1) >= 0) then
        ext%nExtFeatures = tmp(1)
      else
        call error('Error while reading external feature configuration from netstat file. '&
            & // 'Negative number of external features obtained.')
      end if
      call h5ltfx_read_dataset_int_f(external_id, 'indices', ext%indices)
      ! close the external group
      call h5gclose_f(external_id, iErr)
    else
      call error('Error while reading external feature configuration from netstat file. '&
          & // 'Corresponding group is not present.')
    end if

    ! close the netstat group
    call h5gclose_f(netstat_id, iErr)

    ! close the netstat file
    call h5fclose_f(file_id, iErr)

    ! close the hdf5 interface
    call h5close_f(iErr)

  end subroutine readExtFeaturesConfig


  !> Reads architecture of the BPNN header from netstat file.
  subroutine readSubnetArchitecture(fname, type, activation, topology)

    !> filename of netstat file
    character(len=*), intent(in) :: fname

    !> architecture type (currently, only the BPNN topology is available)
    character(len=:), intent(out), allocatable :: type

    !> type of activation functions
    character(len=:), intent(out), allocatable :: activation

    !> topology of the sub-nn's
    integer, intent(out), allocatable :: topology(:)

    !! various specifier flags
    integer(hid_t) :: file_id, netstat_id, bpnn_id, subnet_id

    !! atomic number of species
    integer, allocatable :: atomicNumbers(:)

    !! temporary network topology of current sub-nn
    integer, allocatable :: tmpTopology(:)

    !! temporary type of activation functions of current sub-nn
    character(len=:), allocatable :: tmpActivation

    !! name of current subnetwork
    character(len=:), allocatable :: netname

    !! auxiliary variables
    character(len=100) :: tmp
    integer :: iErr, iNet

    ! currently only the BPNN topology is supported
    type = 'bpnn'

    ! open the hdf5 interface
    call h5open_f(iErr)

    ! open the netstat file
    call h5fopen_f(fname, H5F_ACC_RDONLY_F, file_id, iErr)

    ! open the netstat group
    call h5gopen_f(file_id, 'netstat', netstat_id, iErr)

    ! open the bpnn group
    call h5gopen_f(netstat_id, 'bpnn', bpnn_id, iErr)

    ! read atomic numbers of sub-nn
    call h5ltfx_read_dataset_int_f(bpnn_id, 'atomicnumbers', atomicNumbers)

    do iNet = 1, size(atomicNumbers)
      netname = trim(elementSymbol(atomicNumbers(iNet))) // '-subnetwork'

      ! open the subnetwork group
      call h5gopen_f(bpnn_id, netname, subnet_id, iErr)

      ! read layer topology
      call move_alloc(topology, tmpTopology)
      call h5ltfx_read_dataset_int_f(subnet_id, 'topology', topology)
      if (allocated(tmpTopology)) then
        if (.not. all(shape(tmpTopology) == shape(topology))) then
          call error('Error while reading sub-network topology. Currently only equally shaped'&
              & //' networks are supported.')
        end if
        if (.not. all(tmpTopology == topology)) then
          call error('Error while reading sub-network topology. Currently only equal topologies'&
              & //' are supported.')
        end if
      end if

      ! read transfer function type
      call move_alloc(activation, tmpActivation)
      call h5ltget_attribute_string_f(subnet_id, './', 'activation', tmp, iErr)
      activation = trim(tmp)
      if (allocated(tmpActivation)) then
        if (.not. (trim(tmpActivation) == trim(activation))) then
          call error('Error while reading sub-network transfer functions. Currently only networks'&
              & //' with equal transfer are supported.')
        end if
      end if

      ! close the subnetwork group
      call h5gclose_f(subnet_id, iErr)

    end do

    ! close the bpnn group
    call h5gclose_f(bpnn_id, iErr)

    ! close the netstat group
    call h5gclose_f(netstat_id, iErr)

    ! close the netstat file
    call h5fclose_f(file_id, iErr)

    ! close the hdf5 interface
    call h5close_f(iErr)

  end subroutine readSubnetArchitecture


  !> Writes a BPNN header to netstat file that establishes rudimentary groups.
  subroutine writeBpnnHeader(fname, bpnn, nGlobalTargets, nAtomicTargets)

    !> filename of netstat file
    character(len=*), intent(in) :: fname

    !> Behler-Parrinello-Neural-Network instance
    type(TBpnn), intent(in) :: bpnn

    !> number of system-wide targets of BPNN
    integer, intent(in) :: nGlobalTargets

    !> number of atomic targets of BPNN
    integer, intent(in) :: nAtomicTargets

    !! various specifier flags
    integer(hid_t) :: file_id, netstat_id, bpnn_id, subnet_id, layer_id

    !! name of current subnetwork and layer
    character(len=:), allocatable :: netname, layername

    !! auxiliary variables
    integer(hsize_t) :: dims(1)
    integer :: iErr, iNet, iLayer

    ! open the hdf5 interface
    call h5open_f(iErr)

    ! open the netstat file
    call h5fopen_f(fname, H5F_ACC_RDWR_F, file_id, iErr)

    ! open the netstat group
    call h5gopen_f(file_id, 'netstat', netstat_id, iErr)

    ! create bpnn group
    call h5gcreate_f(netstat_id, 'bpnn', bpnn_id, iErr)

    ! write atomic numbers of sub-nn
    call h5ltfxmake_dataset_int_f(bpnn_id, 'atomicnumbers', bpnn%atomicNumbers)

    ! define the output type, i.e. number of atomic and global targets during training
    dims(1) = 1
    call h5ltset_attribute_int_f(bpnn_id, './', 'nglobaltargets', [nGlobalTargets], dims(1), iErr)
    call h5ltset_attribute_int_f(bpnn_id, './', 'natomictargets', [nAtomicTargets], dims(1), iErr)

    do iNet = 1, bpnn%nSpecies
      netname = trim(elementSymbol(bpnn%atomicNumbers(iNet))) // '-subnetwork'

      ! create subnetwork group
      call h5gcreate_f(bpnn_id, netname, subnet_id, iErr)

      ! write atomic number of sub-nn element
      dims(1) = 1
      call h5ltset_attribute_int_f(subnet_id, './', 'element', bpnn%atomicNumbers(iNet), dims(1),&
          & iErr)

      ! write number of layers of sub-nn
      dims(1) = 1
      call h5ltset_attribute_int_f(subnet_id, './', 'nlayer', [size(bpnn%nets(iNet)%dims)],&
          & dims(1), iErr)

      ! write layer topology
      call h5ltfxmake_dataset_int_f(subnet_id, 'topology', bpnn%dims)

      ! write transfer function type
      call h5ltset_attribute_string_f(subnet_id, './', 'activation', bpnn%nets(iNet)%transferType,&
          & iErr)

      ! write dummy layer structure
      do iLayer = 1, size(bpnn%nets(iNet)%dims) - 1
        layername = 'layer' // i2c(iLayer)

        ! create the layer group
        call h5gcreate_f(subnet_id, layername, layer_id, iErr)

        ! close layer group
        call h5gclose_f(layer_id, iErr)

      end do

      ! close subnetwork group
      call h5gclose_f(subnet_id, iErr)

    end do

    ! close bpnn group
    call h5gclose_f(bpnn_id, iErr)

    ! close netstat group
    call h5gclose_f(netstat_id, iErr)

    ! close the netstat file
    call h5fclose_f(file_id, iErr)

    ! close the hdf5 interface
    call h5close_f(iErr)

  end subroutine writeBpnnHeader


  !> Writes configuration of external atomic input features to netstat file.
  subroutine writeExtFeaturesConfig(fname, ext)

    !> filename of netstat file
    character(len=*), intent(in) :: fname

    !> data type containing variables of the External block
    type(TExternalBlock), intent(in) :: ext

    !! various specifier flags
    integer(hid_t) :: file_id, netstat_id, external_id

    !! auxiliary variables
    integer(hsize_t) :: dim
    integer :: iErr, tExist

    ! open the hdf5 interface
    call h5open_f(iErr)

    ! open the netstat file
    call h5fopen_f(fname, H5F_ACC_RDWR_F, file_id, iErr)

    ! open the netstat group
    call h5gopen_f(file_id, 'netstat', netstat_id, iErr)

    ! check if an external entry is already present
    tExist = h5ltfind_dataset_f(netstat_id, 'external')

    if (tExist == 1) then
      call h5ldelete_f(netstat_id, 'external', iErr)
    end if

    ! create the external group
    call h5gcreate_f(netstat_id, 'external', external_id, iErr)

    ! set the total number of external features
    dim = 1
    call h5ltset_attribute_int_f(external_id, './', 'nextfeatures', [ext%nExtFeatures], dim, iErr)

    ! set the dataset indices of external features
    call h5ltfxmake_dataset_int_f(external_id, 'indices', ext%indices)

    ! close the external group
    call h5gclose_f(external_id, iErr)

    ! close the netstat group
    call h5gclose_f(netstat_id, iErr)

    ! close the netstat file
    call h5fclose_f(file_id, iErr)

    ! close the hdf5 interface
    call h5close_f(iErr)

  end subroutine writeExtFeaturesConfig

end module fnet_netstat
