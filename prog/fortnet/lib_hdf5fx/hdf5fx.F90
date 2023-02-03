!--------------------------------------------------------------------------------------------------!
!  FORTNET: A Behler-Parrinello-Neural-Network Implementation                                      !
!  Copyright (C) 2020 - 2023  T. W. van der Heide                                                  !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Wraps around the high-level HDF5 Fortran API to increase the ease of use.
module fnet_hdf5fx

  use h5lt
  use hdf5

  use dftbp_accuracy, only: dp

  implicit none

  private

  public :: h5ltfx_read_dataset_int_f
  public :: h5ltfx_read_dataset_double_f

  public :: h5ltfxmake_dataset_int_f
  public :: h5ltfxmake_dataset_double_f

  interface h5ltfx_read_dataset_int_f
    module procedure h5ltfx_read_dataset_int1d_f
    module procedure h5ltfx_read_dataset_int2d_f
  end interface h5ltfx_read_dataset_int_f

  interface h5ltfx_read_dataset_double_f
    module procedure h5ltfx_read_dataset_double1d_f
    module procedure h5ltfx_read_dataset_double2d_f
  end interface h5ltfx_read_dataset_double_f

  interface h5ltfxmake_dataset_int_f
    module procedure h5ltfxmake_dataset_int1d_f
    module procedure h5ltfxmake_dataset_int2d_f
  end interface h5ltfxmake_dataset_int_f

  interface h5ltfxmake_dataset_double_f
    module procedure h5ltfxmake_dataset_double1d_f
    module procedure h5ltfxmake_dataset_double2d_f
  end interface h5ltfxmake_dataset_double_f


contains

  !> Wraps the h5ltread_dataset_int_f subroutine for a 1d integer array.
  subroutine h5ltfx_read_dataset_int1d_f(group, dname, array)

    !> group identifier
    integer(hid_t), intent(in) :: group

    !> name of the dataset to read
    character(len=*), intent(in) :: dname

    !> representation of a dataset
    integer, intent(out), allocatable :: array(:)

    !! auxiliary variables
    integer(hsize_t) :: dims(1)
    integer(size_t) :: type_size
    integer :: iErr, type_class

    call h5ltget_dataset_info_f(group, dname, dims, type_class, type_size, iErr)
    allocate(array(dims(1)))
    call h5ltread_dataset_int_f(group, dname, array, dims, iErr)

  end subroutine h5ltfx_read_dataset_int1d_f


  !> Wraps the h5ltread_dataset_int_f subroutine for a 2d integer array.
  subroutine h5ltfx_read_dataset_int2d_f(group, dname, array)

    !> group identifier
    integer(hid_t), intent(in) :: group

    !> name of the dataset to read
    character(len=*), intent(in) :: dname

    !> representation of a dataset
    integer, intent(out), allocatable :: array(:,:)

    !! auxiliary variables
    integer(hsize_t) :: dims(2)
    integer(size_t) :: type_size
    integer :: iErr, type_class

    call h5ltget_dataset_info_f(group, dname, dims, type_class, type_size, iErr)
    allocate(array(dims(1), dims(2)))
    call h5ltread_dataset_int_f(group, dname, array, dims, iErr)

  end subroutine h5ltfx_read_dataset_int2d_f


  !> Wraps the h5lt_read_dataset_double_f subroutine for a real-valued 1d array.
  subroutine h5ltfx_read_dataset_double1d_f(group, dname, array)

    !> group identifier
    integer(hid_t), intent(in) :: group

    !> name of the dataset to read
    character(len=*), intent(in) :: dname

    !> representation of a dataset
    real(dp), intent(out), allocatable :: array(:)

    !! auxiliary variables
    integer(hsize_t) :: dims(1)
    integer(size_t) :: type_size
    integer :: iErr, type_class

    call h5ltget_dataset_info_f(group, dname, dims, type_class, type_size, iErr)
    allocate(array(dims(1)))
    call h5ltread_dataset_double_f(group, dname, array, dims, iErr)

  end subroutine h5ltfx_read_dataset_double1d_f


  !> Wraps the h5lt_read_dataset_double2d_f subroutine for a real-valued 2d array.
  subroutine h5ltfx_read_dataset_double2d_f(group, dname, array)

    !> group identifier
    integer(hid_t), intent(in) :: group

    !> name of the dataset to read
    character(len=*), intent(in) :: dname

    !> representation of a dataset
    real(dp), intent(out), allocatable :: array(:,:)

    !! auxiliary variables
    integer(hsize_t) :: dims(2)
    integer(size_t) :: type_size
    integer :: iErr, type_class

    call h5ltget_dataset_info_f(group, dname, dims, type_class, type_size, iErr)
    allocate(array(dims(1), dims(2)))
    call h5ltread_dataset_double_f(group, dname, array, dims, iErr)

  end subroutine h5ltfx_read_dataset_double2d_f


  !> Wraps the h5ltmake_dataset_int_f subroutine for a 1d integer array.
  subroutine h5ltfxmake_dataset_int1d_f(group, dname, array)

    !> group identifier
    integer(hid_t), intent(in) :: group

    !> name of the dataset to write
    character(len=*), intent(in) :: dname

    !> representation of a dataset
    integer, intent(in) :: array(:)

    !! auxiliary variables
    integer(hsize_t) :: dims(1)
    integer :: iErr

    dims(1) = size(array)
    call h5ltmake_dataset_int_f(group, dname, 1, dims, array, iErr)

  end subroutine h5ltfxmake_dataset_int1d_f


  !> Wraps the h5ltmake_dataset_int_f subroutine for a 2d integer array.
  subroutine h5ltfxmake_dataset_int2d_f(group, dname, array)

    !> group identifier
    integer(hid_t), intent(in) :: group

    !> name of the dataset to write
    character(len=*), intent(in) :: dname

    !> representation of a dataset
    integer, intent(in) :: array(:,:)

    !! auxiliary variables
    integer(hsize_t) :: dims(2)
    integer :: iErr, iDim

    do iDim = 1, 2
      dims(iDim) = size(array, dim=iDim)
    end do

    call h5ltmake_dataset_int_f(group, dname, 2, dims, array, iErr)

  end subroutine h5ltfxmake_dataset_int2d_f


  !> Wraps the h5ltmake_dataset_double_f subroutine for a real-valued 1d array.
  subroutine h5ltfxmake_dataset_double1d_f(group, dname, array)

    !> group identifier
    integer(hid_t), intent(in) :: group

    !> name of the dataset to read
    character(len=*), intent(in) :: dname

    !> representation of a dataset
    real(dp), intent(in) :: array(:)

    !! auxiliary variables
    integer(hsize_t) :: dims(1)
    integer :: iErr

    dims(1) = size(array, dim=1)
    call h5ltmake_dataset_double_f(group, dname, 1, dims, array, iErr)

  end subroutine h5ltfxmake_dataset_double1d_f


  !> Wraps the h5ltmake_dataset_double_f subroutine for a real-valued 2d array.
  subroutine h5ltfxmake_dataset_double2d_f(group, dname, array)

    !> group identifier
    integer(hid_t), intent(in) :: group

    !> name of the dataset to read
    character(len=*), intent(in) :: dname

    !> representation of a dataset
    real(dp), intent(in) :: array(:,:)

    !! auxiliary variables
    integer(hsize_t) :: dims(2)
    integer :: iErr, iDim

    do iDim = 1, 2
      dims(iDim) = size(array, dim=iDim)
    end do

    call h5ltmake_dataset_double_f(group, dname, 2, dims, array, iErr)

  end subroutine h5ltfxmake_dataset_double2d_f

end module fnet_hdf5fx
