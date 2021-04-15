!--------------------------------------------------------------------------------------------------!
!  FORTNET: A Behler-Parrinello-Neural-Network Implementation                                      !
!  Copyright (C) 2020 - 2021  T. W. van der Heide                                                  !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

module fnet_loss

  use dftbp_accuracy, only : dp
  use fnet_nestedtypes, only : TPredicts, TRealArray2D

  implicit none

  private

  public :: lossFunc
  public :: minError, maxError
  public :: deviation, maLoss, msLoss, mslLoss, rmsLoss


  interface

    pure function lossFunc(predicts, targets) result(loss)

      use dftbp_accuracy, only: dp
      use fnet_nestedtypes, only : TPredicts, TRealArray2D

      implicit none

      !> neural network predictions
      type(TPredicts), intent(in) :: predicts

      !> target reference data
      type(TRealArray2D), intent(in) :: targets(:)

      !> summed root mean square loss of predictions, in comparison to targets
      real(dp) :: loss

    end function lossFunc

  end interface


contains


  !> simple deviation between predictions and targets
  pure function deviation(predicts, targets) result(dev)

    !> neural network predictions
    real(dp), intent(in) :: predicts(:)

    !> target reference data
    real(dp), intent(in) :: targets(:)

    !> deviation of predictions, in comparison to targets
    real(dp), allocatable :: dev(:)

    dev = predicts - targets

  end function deviation


  !> mean absolute loss function
  pure function simpleMaLoss(predicts, targets) result(loss)

    !> neural network predictions
    real(dp), intent(in) :: predicts(:)

    !> target reference data
    real(dp), intent(in) :: targets(:)

    !> summed mean absolute loss of predictions, in comparison to targets
    real(dp) :: loss

    loss = sum(abs(targets - predicts)) / real(size(predicts), dp)

  end function simpleMaLoss


  !> mean squared loss function
  pure function simpleMsLoss(predicts, targets) result(loss)

    !> neural network predictions
    real(dp), intent(in) :: predicts(:)

    !> target reference data
    real(dp), intent(in) :: targets(:)

    !> summed mean squared loss of predictions, in comparison to targets
    real(dp) :: loss

    loss = sum((targets - predicts)**2) / real(size(predicts), dp)

  end function simpleMsLoss


  !> mean squared logarithmic loss function
  pure function simpleMslLoss(predicts, targets) result(loss)

    !> neural network predictions
    real(dp), intent(in) :: predicts(:)

    !> target reference data
    real(dp), intent(in) :: targets(:)

    !> summed mean squared logarithmic loss of predictions, in comparison to targets
    real(dp) :: loss

    loss = sum((log(targets + 1.0_dp) - log(predicts + 1.0_dp))**2) / real(size(predicts), dp)

  end function simpleMslLoss


  !> root mean square loss function
  pure function simpleRmsLoss(predicts, targets) result(loss)

    !> neural network predictions
    real(dp), intent(in) :: predicts(:)

    !> target reference data
    real(dp), intent(in) :: targets(:)

    !> summed root mean square loss of predictions, in comparison to targets
    real(dp) :: loss

    loss = sqrt(sum((targets - predicts)**2) / real(size(predicts), dp))

  end function simpleRmsLoss


  !> mean absolute loss function
  pure function maLoss(predicts, targets) result(loss)

    !> neural network predictions
    type(TPredicts), intent(in) :: predicts

    !> target reference data
    type(TRealArray2D), intent(in) :: targets(:)

    !> summed mean absolute loss of predictions, in comparison to targets
    real(dp) :: loss

    !> temporary loss storage
    real(dp) :: tmpLoss

    !> auxiliary variables
    integer :: iSys, iAtom, nValues

    loss = 0.0_dp
    nValues = 0

    do iSys = 1, size(predicts%sys)
      do iAtom = 1, size(predicts%sys(iSys)%array, dim=2)
        tmpLoss = simpleMaLoss(predicts%sys(iSys)%array(:, iAtom), targets(iSys)%array(:, iAtom))
        loss = loss + tmpLoss
        nValues = nValues + 1
      end do
    end do

    loss = loss / real(nValues, dp)

  end function maLoss


  !> mean squared loss function
  pure function msLoss(predicts, targets) result(loss)

    !> neural network predictions
    type(TPredicts), intent(in) :: predicts

    !> target reference data
    type(TRealArray2D), intent(in) :: targets(:)

    !> summed mean squared loss of predictions, in comparison to targets
    real(dp) :: loss

    !> temporary loss storage
    real(dp) :: tmpLoss

    !> auxiliary variables
    integer :: iSys, iAtom, nValues

    loss = 0.0_dp
    nValues = 0

    do iSys = 1, size(predicts%sys)
      do iAtom = 1, size(predicts%sys(iSys)%array, dim=2)
        tmpLoss = simpleMsLoss(predicts%sys(iSys)%array(:, iAtom), targets(iSys)%array(:, iAtom))
        loss = loss + tmpLoss
        nValues = nValues + 1
      end do
    end do

    loss = loss / real(nValues, dp)

  end function msLoss


  !> mean squared logarithmic loss function
  pure function mslLoss(predicts, targets) result(loss)

    !> neural network predictions
    type(TPredicts), intent(in) :: predicts

    !> target reference data
    type(TRealArray2D), intent(in) :: targets(:)

    !> summed mean square logarithmic loss of predictions, in comparison to targets
    real(dp) :: loss

    !> temporary loss storage
    real(dp) :: tmpLoss

    !> auxiliary variables
    integer :: iSys, iAtom, nValues

    loss = 0.0_dp
    nValues = 0

    do iSys = 1, size(predicts%sys)
      do iAtom = 1, size(predicts%sys(iSys)%array, dim=2)
        tmpLoss = simpleMslLoss(predicts%sys(iSys)%array(:, iAtom), targets(iSys)%array(:, iAtom))
        loss = loss + tmpLoss
        nValues = nValues + 1
      end do
    end do

    loss = loss / real(nValues, dp)

  end function mslLoss


  !> root mean square loss function
  pure function rmsLoss(predicts, targets) result(loss)

    !> neural network predictions
    type(TPredicts), intent(in) :: predicts

    !> target reference data
    type(TRealArray2D), intent(in) :: targets(:)

    !> summed root mean square loss of predictions, in comparison to targets
    real(dp) :: loss

    !> temporary loss storage
    real(dp) :: tmpLoss

    !> auxiliary variables
    integer :: iSys, iAtom, nValues

    loss = 0.0_dp
    nValues = 0

    do iSys = 1, size(predicts%sys)
      do iAtom = 1, size(predicts%sys(iSys)%array, dim=2)
        tmpLoss = simpleRmsLoss(predicts%sys(iSys)%array(:, iAtom), targets(iSys)%array(:, iAtom))
        loss = loss + tmpLoss
        nValues = nValues + 1
      end do
    end do

    loss = loss / real(nValues, dp)

  end function rmsLoss


  !> calculates the absolute minium deviation between predictions and targets
  pure function minError(predicts, targets) result(err)

    !> neural network predictions
    type(TPredicts), intent(in) :: predicts

    !> target reference data
    type(TRealArray2D), intent(in) :: targets(:)

    !> absolute minium deviation between predictions and targets
    real(dp) :: err

    !> temporary error storage
    real(dp) :: tmpErr

    !> auxiliary variable
    integer :: iSys

    err = minval(abs(predicts%sys(1)%array - targets(1)%array))

    do iSys = 2, size(predicts%sys)
      tmpErr = minval(abs(predicts%sys(iSys)%array - targets(iSys)%array))
      if (tmpErr < err) then
        err = tmpErr
      end if
    end do

  end function minError


  !> calculates the absolute maximum deviation between predictions and targets
  pure function maxError(predicts, targets) result(err)

    !> neural network predictions
    type(TPredicts), intent(in) :: predicts

    !> target reference data
    type(TRealArray2D), intent(in) :: targets(:)

    !> absolute maximum deviation between predictions and targets
    real(dp) :: err

    !> temporary error storage
    real(dp) :: tmpErr

    !> auxiliary variable
    integer :: iSys

    err = maxval(abs(predicts%sys(1)%array - targets(1)%array))

    do iSys = 2, size(predicts%sys)
      tmpErr = maxval(abs(predicts%sys(iSys)%array - targets(iSys)%array))
      if (tmpErr > err) then
        err = tmpErr
      end if
    end do

  end function maxError

end module fnet_loss
