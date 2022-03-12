!--------------------------------------------------------------------------------------------------!
!  FORTNET: A Behler-Parrinello-Neural-Network Implementation                                      !
!  Copyright (C) 2020 - 2022  T. W. van der Heide                                                  !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Implements the Fortnet output file format.
module fnet_fnetout

  use h5lt
  use hdf5

  use dftbp_message, only : error
  use dftbp_accuracy, only: dp
  use dftbp_charmanip, only : i2c

  use fnet_forces, only : TForces
  use fnet_nestedtypes, only : TRealArray1D, TRealArray2D, TPredicts
  use fnet_hdf5fx, only : h5ltfxmake_dataset_double_f

  implicit none

  private

  public :: writeFnetout


contains

  !> Writes the obtained predictions (+ targets) to an Fnetout file.
  subroutine writeFnetout(fname, mode, globalTargets, atomicTargets, output, forces,&
      & tForcesSupplied)

    !> filename (will be fnetout.hdf5)
    character(len=*), intent(in) :: fname

    !> mode of calculation (train, validate, run)
    character(len=*), intent(in) :: mode

    !> system-wide target values for training
    type(TRealArray1D), intent(in) :: globalTargets(:)

    !> atomic target values for training
    type(TRealArray2D), intent(in) :: atomicTargets(:)

    !> obtained output values of network
    type(TPredicts), intent(in) :: output

    !> obtained atomic forces
    type(TForces), intent(in) :: forces

    !> true, if writing out atomic forces is desired
    logical, intent(in) :: tForcesSupplied

    !> file and group identification
    integer(hid_t) :: file_id, fnetout_id, output_id, datapoint_id

    !> temporary storage container
    real(dp), allocatable :: tmpOutput(:,:)

    !> total number of training targets
    integer :: nTargets

    !> auxiliary variables
    integer(hsize_t) :: dims(1)
    integer :: iSys, iErr

    if (mode /= 'validate' .and. mode /= 'predict') then
      call error('Invalid program running mode selected.')
    end if

    nTargets = output%nGlobalTargets + output%nAtomicTargets

    if ((mode == 'validate') .and. (nTargets == 0)) then
      call error('Validation mode only valid in combination with target data.')
    end if

    ! open the hdf5 interface
    call h5open_f(iErr)

    ! create the fnetout file
    call h5fcreate_f(fname, H5F_ACC_TRUNC_F, file_id, iErr)

    ! create the fnetout group
    call h5gcreate_f(file_id, 'fnetout', fnetout_id, iErr)

    ! write program running mode
    call h5ltset_attribute_string_f(fnetout_id, './', 'mode', mode, iErr)

    ! create the output group
    call h5gcreate_f(fnetout_id, 'output', output_id, iErr)

    ! write number of datapoints
    dims(1) = 1
    call h5ltset_attribute_int_f(output_id, './', 'ndatapoints', [output%nDatapoints], dims(1),&
        & iErr)

    ! indicate whether atomic forces are supplied
    dims(1) = 1
    if (tForcesSupplied) then
      call h5ltset_attribute_int_f(output_id, './', 'tforces', [1], dims(1), iErr)
    else
      call h5ltset_attribute_int_f(output_id, './', 'tforces', [0], dims(1), iErr)
    end if

    ! write number of predictions (system-wide/atomic)
    dims(1) = 1
    call h5ltset_attribute_int_f(output_id, './', 'nglobaltargets', [output%nGlobalTargets],&
        & dims(1), iErr)
    call h5ltset_attribute_int_f(output_id, './', 'natomictargets', [output%nAtomicTargets],&
        & dims(1), iErr)

    select case(mode)
    case('validate')
      do iSys = 1, output%nDatapoints
        ! create the datapoint group
        call h5gcreate_f(output_id, 'datapoint' // trim(i2c(iSys)), datapoint_id, iErr)
        ! write targets and predictions to file
        if (output%nGlobalTargets > 0) then
          call h5ltfxmake_dataset_double_f(datapoint_id, 'globalpredictions',&
              & sum(output%sys(iSys)%array(1:output%nGlobalTargets, :), dim=2))
        end if
        call h5ltfxmake_dataset_double_f(datapoint_id, 'rawpredictions', output%sys(iSys)%array)
        call h5ltfxmake_dataset_double_f(datapoint_id, 'globaltargets', globalTargets(iSys)%array)
        call h5ltfxmake_dataset_double_f(datapoint_id, 'atomictargets', atomicTargets(iSys)%array)
        if (tForcesSupplied) then
          ! write atomic forces to file
          call h5ltfxmake_dataset_double_f(datapoint_id, 'forces', forces%geos(iSys)%array)
        end if
        ! close the datapoint group
        call h5gclose_f(datapoint_id, iErr)
      end do
    case('predict')
      do iSys = 1, output%nDatapoints
        ! create the datapoint group
        call h5gcreate_f(output_id, 'datapoint' // trim(i2c(iSys)), datapoint_id, iErr)
        if (output%nGlobalTargets > 0) then
          call h5ltfxmake_dataset_double_f(datapoint_id, 'globalpredictions',&
              & sum(output%sys(iSys)%array(1:output%nGlobalTargets, :), dim=2))
        end if
        call h5ltfxmake_dataset_double_f(datapoint_id, 'rawpredictions', output%sys(iSys)%array)
        if (tForcesSupplied) then
          ! write atomic forces to file
          call h5ltfxmake_dataset_double_f(datapoint_id, 'forces', forces%geos(iSys)%array)
        end if
        ! close the datapoint group
        call h5gclose_f(datapoint_id, iErr)
      end do
    end select

    ! close the output group
    call h5gclose_f(output_id, iErr)

    ! close the fnetout group
    call h5gclose_f(fnetout_id, iErr)

    ! close the fnetout file
    call h5fclose_f(file_id, iErr)

    ! close the hdf5 interface
    call h5close_f(iErr)

  end subroutine writeFnetout

end module fnet_fnetout
