!--------------------------------------------------------------------------------------------------!
!  FORTNET: A Behler-Parrinello-Neural-Network Implementation                                      !
!  Copyright (C) 2020 - 2025  T. W. van der Heide                                                  !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Provides several loss functions as well as their gradients.
module fnet_loss

  use dftbp_accuracy, only : dp
  use fnet_nestedtypes, only : TPredicts, TRealArray1D, TRealArray2D

  implicit none

  private

  public :: lossFunc, lossGradientFunc
  public :: minError, maxError
  public :: deviation, maLoss, mapLoss, msLoss, rmsLoss
  public :: maGradients, mapGradients, msGradients, rmsGradients

  public :: TRegularizationBlock
  public :: reguFunc, nullReguLoss, lassoReguLoss, ridgeReguLoss, elasticNetReguLoss


  abstract interface

    pure function lossFunc(predicts, globalTargets, atomicTargets, atomicWeights, weights)&
        & result(loss)

      use dftbp_accuracy, only: dp
      use fnet_nestedtypes, only : TPredicts, TRealArray1D, TRealArray2D

      implicit none

      !> neural network predictions
      type(TPredicts), intent(in) :: predicts

      !> system-wide target reference data
      type(TRealArray1D), intent(in) :: globalTargets(:)

      !> atomic target reference data
      type(TRealArray2D), intent(in) :: atomicTargets(:)

      !> contains atomic gradient weights
      type(TRealArray1D), intent(in) :: atomicWeights(:)

      !> optional weighting of individual datapoints
      integer, intent(in), optional :: weights(:)

      !! total loss of predictions, in comparison to targets
      real(dp) :: loss

    end function lossFunc


    pure function lossGradientFunc(predicts, targets) result(grads)

      use dftbp_accuracy, only: dp

      implicit none

      !> neural network predictions
      real(dp), intent(in) :: predicts(:,:)

      !> target reference data
      real(dp), intent(in) :: targets(:,:)

      !! loss gradients w.r.t predictions and targets
      real(dp) :: grads(size(targets, dim=1), size(targets, dim=2))

    end function lossGradientFunc


    pure function reguFunc(weights, lambda, alpha) result(loss)

      use dftbp_accuracy, only: dp

      implicit none

      !> serialized weights of the network
      real(dp), intent(in) :: weights(:)

      !> strength of the regularization
      real(dp), intent(in) :: lambda

      !> elastic net weighting between lasso and ridge (0 = ridge, 1 = lasso)
      real(dp), intent(in) :: alpha

      !! total regularization loss
      real(dp) :: loss

    end function reguFunc

  end interface


  !> Data type containing variables of the Regularization block.
  type TRegularizationBlock

    !> regularization strength
    real(dp) :: strength

    !> elastic net weighting between lasso and ridge (0 = ridge, 1 = lasso)
    real(dp) :: alpha

    !> regularization type (lasso, ridge, elasticnet)
    character(len=:), allocatable :: type

  end type TRegularizationBlock


contains

  !> Calculates a dummy, zero-valued regularization loss contribution.
  pure function nullReguLoss(weights, lambda, alpha) result(loss)

    !> serialized weights of the network
    real(dp), intent(in) :: weights(:)

    !> strength of the regularization
    real(dp), intent(in) :: lambda

    !> elastic net weighting between lasso and ridge (0 = ridge, 1 = lasso)
    real(dp), intent(in) :: alpha

    !! total regularization loss
    real(dp) :: loss

    loss = 0.0_dp

  end function nullReguLoss


  !> Calculates the L1 regularization (lasso regression) loss contribution.
  pure function lassoReguLoss(weights, lambda, alpha) result(loss)

    !> serialized weights of the network
    real(dp), intent(in) :: weights(:)

    !> strength of the regularization
    real(dp), intent(in) :: lambda

    !> elastic net weighting between lasso and ridge (0 = ridge, 1 = lasso)
    real(dp), intent(in) :: alpha

    !! total regularization loss
    real(dp) :: loss

    loss = lambda / real(size(weights), dp) * sum(abs(weights))

  end function lassoReguLoss


  !> Calculates the L2 regularization (ridge regression) loss contribution.
  pure function ridgeReguLoss(weights, lambda, alpha) result(loss)

    !> serialized weights of the network
    real(dp), intent(in) :: weights(:)

    !> strength of the regularization
    real(dp), intent(in) :: lambda

    !> elastic net weighting between lasso and ridge (0 = ridge, 1 = lasso)
    real(dp), intent(in) :: alpha

    !! total regularization loss
    real(dp) :: loss

    loss = lambda / (2.0_dp * real(size(weights), dp)) * sum(weights**2)

  end function ridgeReguLoss


  !> Calculates an elastic net regularization loss contribution.
  pure function elasticNetReguLoss(weights, lambda, alpha) result(loss)

    !> serialized weights of the network
    real(dp), intent(in) :: weights(:)

    !> strength of the regularization
    real(dp), intent(in) :: lambda

    !> elastic net weighting between lasso and ridge (0 = ridge, 1 = lasso)
    real(dp), intent(in) :: alpha

    !! total regularization loss
    real(dp) :: loss

    loss = lambda / real(size(weights), dp) * ((1 - alpha) / 2.0_dp * sum(weights**2)&
        & + alpha * sum(abs(weights)))

  end function elasticNetReguLoss


  !> Calculates the simple deviation between predictions and targets.
  pure function deviation(predicts, targets) result(dev)

    !> neural network predictions
    real(dp), intent(in) :: predicts(:)

    !> target reference data
    real(dp), intent(in) :: targets(:)

    !! deviation of predictions, in comparison to targets
    real(dp) :: dev(size(targets))

    dev(:) = predicts - targets

  end function deviation


  !> Calculates the mean absolute error derivation w.r.t predictions and targets.
  pure function maGradients(predicts, targets) result(grads)

    !> neural network predictions
    real(dp), intent(in) :: predicts(:,:)

    !> target reference data
    real(dp), intent(in) :: targets(:,:)

    !! loss gradients w.r.t predictions and targets
    real(dp) :: grads(size(targets, dim=1), size(targets, dim=2))

    grads(:,:) = (predicts - targets) / abs(predicts - targets)

  end function maGradients


  !> Calculates the mean absolute percentage error derivation w.r.t predictions and targets.
  pure function mapGradients(predicts, targets) result(grads)

    !> neural network predictions
    real(dp), intent(in) :: predicts(:,:)

    !> target reference data
    real(dp), intent(in) :: targets(:,:)

    !! loss gradients w.r.t predictions and targets
    real(dp) :: grads(size(targets, dim=1), size(targets, dim=2))

    grads(:,:) = 100.0_dp * (predicts - targets) / (targets**2 * abs(predicts / targets - 1.0_dp))

  end function mapGradients


  !> Calculates the mean squared error derivation w.r.t predictions and targets.
  pure function msGradients(predicts, targets) result(grads)

    !> neural network predictions
    real(dp), intent(in) :: predicts(:,:)

    !> target reference data
    real(dp), intent(in) :: targets(:,:)

    !! loss gradients w.r.t predictions and targets
    real(dp) :: grads(size(targets, dim=1), size(targets, dim=2))

    grads(:,:) = 2.0_dp * (predicts - targets)

  end function msGradients


  !> Calculates the root mean squared error derivation w.r.t predictions and targets.
  pure function rmsGradients(predicts, targets) result(grads)

    !> neural network predictions
    real(dp), intent(in) :: predicts(:,:)

    !> target reference data
    real(dp), intent(in) :: targets(:,:)

    !! loss gradients w.r.t predictions and targets
    real(dp) :: grads(size(targets, dim=1), size(targets, dim=2))

    grads(:,:) = (predicts - targets) / (sqrt((predicts - targets)**2))

  end function rmsGradients


  !> Calculates the mean absolute loss function.
  pure function simpleMaLoss(predicts, targets) result(loss)

    !> neural network predictions
    real(dp), intent(in) :: predicts(:)

    !> target reference data
    real(dp), intent(in) :: targets(:)

    !! summed mean absolute loss of predictions, in comparison to targets
    real(dp) :: loss

    loss = sum(abs(targets - predicts)) / real(size(predicts), dp)

  end function simpleMaLoss


  !> Calculates the mean absolute percentage loss function.
  pure function simpleMapLoss(predicts, targets) result(loss)

    !> neural network predictions
    real(dp), intent(in) :: predicts(:)

    !> target reference data
    real(dp), intent(in) :: targets(:)

    !! summed mean absolute percentage loss of predictions, in comparison to targets
    real(dp) :: loss

    loss = 100.0_dp * sum(abs((targets - predicts) / targets)) / real(size(predicts), dp)

  end function simpleMapLoss


  !> Calculates the mean squared loss function.
  pure function simpleMsLoss(predicts, targets) result(loss)

    !> neural network predictions
    real(dp), intent(in) :: predicts(:)

    !> target reference data
    real(dp), intent(in) :: targets(:)

    !! summed mean squared loss of predictions, in comparison to targets
    real(dp) :: loss

    loss = sum((targets - predicts)**2) / real(size(predicts), dp)

  end function simpleMsLoss


  !> Calculates the mean squared logarithmic loss function.
  pure function simpleMslLoss(predicts, targets) result(loss)

    !> neural network predictions
    real(dp), intent(in) :: predicts(:)

    !> target reference data
    real(dp), intent(in) :: targets(:)

    !! summed mean squared logarithmic loss of predictions, in comparison to targets
    real(dp) :: loss

    loss = sum((log(targets + 1.0_dp) - log(predicts + 1.0_dp))**2) / real(size(predicts), dp)

  end function simpleMslLoss


  !> Calculates the root mean square loss function.
  pure function simpleRmsLoss(predicts, targets) result(loss)

    !> neural network predictions
    real(dp), intent(in) :: predicts(:)

    !> target reference data
    real(dp), intent(in) :: targets(:)

    !! summed root mean square loss of predictions, in comparison to targets
    real(dp) :: loss

    loss = sqrt(sum((targets - predicts)**2) / real(size(predicts), dp))

  end function simpleRmsLoss


  !> Calculates the mean absolute loss function.
  pure function maLoss(predicts, globalTargets, atomicTargets, atomicWeights, weights) result(loss)

    !> neural network predictions
    type(TPredicts), intent(in) :: predicts

    !> system-wide target reference data
    type(TRealArray1D), intent(in) :: globalTargets(:)

    !> atomic target reference data
    type(TRealArray2D), intent(in) :: atomicTargets(:)

    !> contains atomic gradient weights
    type(TRealArray1D), intent(in) :: atomicWeights(:)

    !> optional weighting of individual datapoints
    integer, intent(in), optional :: weights(:)

    !! temporary real valued storage for summed up system-wide predictions
    real(dp), allocatable :: globalPredicts(:)

    !! summed mean absolute loss of predictions, in comparison to targets
    real(dp) :: loss

    !! temporary loss storage
    real(dp) :: tmpLoss

    !! weighting of individual datapoints
    integer, allocatable :: weighting(:)

    !! auxiliary variables
    integer :: iSys, iAtom
    real(dp) :: nValues

    allocate(weighting(size(predicts%sys)))

    if (present(weights)) then
      weighting(:) = weights
    else
      weighting(:) = 1
    end if

    loss = 0.0_dp
    nValues = 0.0_dp

    if (predicts%nGlobalTargets > 0) then
      do iSys = 1, predicts%nDatapoints
        globalPredicts = sum(predicts%sys(iSys)%array(1:predicts%nGlobalTargets, :), dim=2)
        tmpLoss = simpleMaLoss(globalPredicts, globalTargets(iSys)%array)
        ! multiply loss with datapoint weighting and number of atoms contributing to the prediction
        loss = loss + real(weighting(iSys), dp) * sum(atomicWeights(iSys)%array) * tmpLoss
        nValues = nValues + real(weighting(iSys), dp) * sum(atomicWeights(iSys)%array)
      end do
    end if

    if (predicts%nAtomicTargets > 0) then
      do iSys = 1, predicts%nDatapoints
        do iAtom = 1, size(predicts%sys(iSys)%array, dim=2)
          tmpLoss = simpleMaLoss(predicts%sys(iSys)%array(predicts%nGlobalTargets + 1:, iAtom),&
              & atomicTargets(iSys)%array(:, iAtom))
          loss = loss + real(weighting(iSys), dp) * atomicWeights(iSys)%array(iAtom) * tmpLoss
          nValues = nValues + real(weighting(iSys), dp) * atomicWeights(iSys)%array(iAtom)
        end do
      end do
    end if

    loss = loss / nValues

  end function maLoss


  !> Calculates the mean absolute percentage loss function.
  pure function mapLoss(predicts, globalTargets, atomicTargets, atomicWeights, weights) result(loss)

    !> neural network predictions
    type(TPredicts), intent(in) :: predicts

    !> system-wide target reference data
    type(TRealArray1D), intent(in) :: globalTargets(:)

    !> atomic target reference data
    type(TRealArray2D), intent(in) :: atomicTargets(:)

    !> contains atomic gradient weights
    type(TRealArray1D), intent(in) :: atomicWeights(:)

    !> optional weighting of individual datapoints
    integer, intent(in), optional :: weights(:)

    !! temporary real valued storage for summed up system-wide predictions
    real(dp), allocatable :: globalPredicts(:)

    !! summed mean absolute loss of predictions, in comparison to targets
    real(dp) :: loss

    !! temporary loss storage
    real(dp) :: tmpLoss

    !! weighting of individual datapoints
    integer, allocatable :: weighting(:)

    !! auxiliary variables
    integer :: iSys, iAtom
    real(dp) :: nValues

    allocate(weighting(size(predicts%sys)))

    if (present(weights)) then
      weighting(:) = weights
    else
      weighting(:) = 1
    end if

    loss = 0.0_dp
    nValues = 0.0_dp

    if (predicts%nGlobalTargets > 0) then
      do iSys = 1, predicts%nDatapoints
        globalPredicts = sum(predicts%sys(iSys)%array(1:predicts%nGlobalTargets, :), dim=2)
        tmpLoss = simpleMapLoss(globalPredicts, globalTargets(iSys)%array)
        ! multiply loss with datapoint weighting and number of atoms contributing to the prediction
        loss = loss + real(weighting(iSys), dp) * sum(atomicWeights(iSys)%array) * tmpLoss
        nValues = nValues + real(weighting(iSys), dp) * sum(atomicWeights(iSys)%array)
      end do
    end if

    if (predicts%nAtomicTargets > 0) then
      do iSys = 1, predicts%nDatapoints
        do iAtom = 1, size(predicts%sys(iSys)%array, dim=2)
          tmpLoss = simpleMapLoss(predicts%sys(iSys)%array(predicts%nGlobalTargets + 1:, iAtom),&
              & atomicTargets(iSys)%array(:, iAtom))
          loss = loss + real(weighting(iSys), dp) * atomicWeights(iSys)%array(iAtom) * tmpLoss
          nValues = nValues + real(weighting(iSys), dp) * atomicWeights(iSys)%array(iAtom)
        end do
      end do
    end if

    loss = loss / nValues

  end function mapLoss


  !> Calculates the mean squared loss function.
  pure function msLoss(predicts, globalTargets, atomicTargets, atomicWeights, weights) result(loss)

    !> neural network predictions
    type(TPredicts), intent(in) :: predicts

    !> system-wide target reference data
    type(TRealArray1D), intent(in) :: globalTargets(:)

    !> atomic target reference data
    type(TRealArray2D), intent(in) :: atomicTargets(:)

    !> contains atomic gradient weights
    type(TRealArray1D), intent(in) :: atomicWeights(:)

    !> optional weighting of individual datapoints
    integer, intent(in), optional :: weights(:)

    !! temporary real valued storage for summed up system-wide predictions
    real(dp), allocatable :: globalPredicts(:)

    !! summed mean absolute loss of predictions, in comparison to targets
    real(dp) :: loss

    !! temporary loss storage
    real(dp) :: tmpLoss

    !! weighting of individual datapoints
    integer, allocatable :: weighting(:)

    !! auxiliary variables
    integer :: iSys, iAtom
    real(dp) :: nValues

    allocate(weighting(size(predicts%sys)))

    if (present(weights)) then
      weighting(:) = weights
    else
      weighting(:) = 1
    end if

    loss = 0.0_dp
    nValues = 0.0_dp

    if (predicts%nGlobalTargets > 0) then
      do iSys = 1, predicts%nDatapoints
        globalPredicts = sum(predicts%sys(iSys)%array(1:predicts%nGlobalTargets, :), dim=2)
        tmpLoss = simpleMsLoss(globalPredicts, globalTargets(iSys)%array)
        ! multiply loss with datapoint weighting and number of atoms contributing to the prediction
        loss = loss + real(weighting(iSys), dp) * sum(atomicWeights(iSys)%array) * tmpLoss
        nValues = nValues + real(weighting(iSys), dp) * sum(atomicWeights(iSys)%array)
      end do
    end if

    if (predicts%nAtomicTargets > 0) then
      do iSys = 1, predicts%nDatapoints
        do iAtom = 1, size(predicts%sys(iSys)%array, dim=2)
          tmpLoss = simpleMsLoss(predicts%sys(iSys)%array(predicts%nGlobalTargets + 1:, iAtom),&
              & atomicTargets(iSys)%array(:, iAtom))
          loss = loss + real(weighting(iSys), dp) * atomicWeights(iSys)%array(iAtom) * tmpLoss
          nValues = nValues + real(weighting(iSys), dp) * atomicWeights(iSys)%array(iAtom)
        end do
      end do
    end if

    loss = loss / nValues

  end function msLoss


  !> Calculates the mean squared logarithmic loss function.
  pure function mslLoss(predicts, globalTargets, atomicTargets, atomicWeights, weights) result(loss)

    !> neural network predictions
    type(TPredicts), intent(in) :: predicts

    !> system-wide target reference data
    type(TRealArray1D), intent(in) :: globalTargets(:)

    !> atomic target reference data
    type(TRealArray2D), intent(in) :: atomicTargets(:)

    !> contains atomic gradient weights
    type(TRealArray1D), intent(in) :: atomicWeights(:)

    !> optional weighting of individual datapoints
    integer, intent(in), optional :: weights(:)

    !! temporary real valued storage for summed up system-wide predictions
    real(dp), allocatable :: globalPredicts(:)

    !! summed mean absolute loss of predictions, in comparison to targets
    real(dp) :: loss

    !! temporary loss storage
    real(dp) :: tmpLoss

    !! weighting of individual datapoints
    integer, allocatable :: weighting(:)

    !! auxiliary variables
    integer :: iSys, iAtom
    real(dp) :: nValues

    allocate(weighting(size(predicts%sys)))

    if (present(weights)) then
      weighting(:) = weights
    else
      weighting(:) = 1
    end if

    loss = 0.0_dp
    nValues = 0.0_dp

    if (predicts%nGlobalTargets > 0) then
      do iSys = 1, predicts%nDatapoints
        globalPredicts = sum(predicts%sys(iSys)%array(1:predicts%nGlobalTargets, :), dim=2)
        tmpLoss = simpleMslLoss(globalPredicts, globalTargets(iSys)%array)
        ! multiply loss with datapoint weighting and number of atoms contributing to the prediction
        loss = loss + real(weighting(iSys), dp) * sum(atomicWeights(iSys)%array) * tmpLoss
        nValues = nValues + real(weighting(iSys), dp) * sum(atomicWeights(iSys)%array)
      end do
    end if

    if (predicts%nAtomicTargets > 0) then
      do iSys = 1, predicts%nDatapoints
        do iAtom = 1, size(predicts%sys(iSys)%array, dim=2)
          tmpLoss = simpleMslLoss(predicts%sys(iSys)%array(predicts%nGlobalTargets + 1:, iAtom),&
              & atomicTargets(iSys)%array(:, iAtom))
          loss = loss + real(weighting(iSys), dp) * atomicWeights(iSys)%array(iAtom) * tmpLoss
          nValues = nValues + real(weighting(iSys), dp) * atomicWeights(iSys)%array(iAtom)
        end do
      end do
    end if

    loss = loss / nValues

  end function mslLoss


  !> Calculates the root mean square loss function.
  pure function rmsLoss(predicts, globalTargets, atomicTargets, atomicWeights, weights) result(loss)

    !> neural network predictions
    type(TPredicts), intent(in) :: predicts

    !> system-wide target reference data
    type(TRealArray1D), intent(in) :: globalTargets(:)

    !> atomic target reference data
    type(TRealArray2D), intent(in) :: atomicTargets(:)

    !> contains atomic gradient weights
    type(TRealArray1D), intent(in) :: atomicWeights(:)

    !> optional weighting of individual datapoints
    integer, intent(in), optional :: weights(:)

    !! temporary real valued storage for summed up system-wide predictions
    real(dp), allocatable :: globalPredicts(:)

    !! summed mean absolute loss of predictions, in comparison to targets
    real(dp) :: loss

    !! temporary loss storage
    real(dp) :: tmpLoss

    !! weighting of individual datapoints
    integer, allocatable :: weighting(:)

    !! auxiliary variables
    integer :: iSys, iAtom
    real(dp) :: nValues

    allocate(weighting(size(predicts%sys)))

    if (present(weights)) then
      weighting(:) = weights
    else
      weighting(:) = 1
    end if

    loss = 0.0_dp
    nValues = 0.0_dp

    if (predicts%nGlobalTargets > 0) then
      do iSys = 1, predicts%nDatapoints
        globalPredicts = sum(predicts%sys(iSys)%array(1:predicts%nGlobalTargets, :), dim=2)
        tmpLoss = simpleRmsLoss(globalPredicts, globalTargets(iSys)%array)
        ! multiply loss with datapoint weighting and number of atoms contributing to the prediction
        loss = loss + real(weighting(iSys), dp) * sum(atomicWeights(iSys)%array) * tmpLoss
        nValues = nValues + real(weighting(iSys), dp) * sum(atomicWeights(iSys)%array)
      end do
    end if

    if (predicts%nAtomicTargets > 0) then
      do iSys = 1, predicts%nDatapoints
        do iAtom = 1, size(predicts%sys(iSys)%array, dim=2)
          tmpLoss = simpleRmsLoss(predicts%sys(iSys)%array(predicts%nGlobalTargets + 1:, iAtom),&
              & atomicTargets(iSys)%array(:, iAtom))
          loss = loss + real(weighting(iSys), dp) * atomicWeights(iSys)%array(iAtom) * tmpLoss
          nValues = nValues + real(weighting(iSys), dp) * atomicWeights(iSys)%array(iAtom)
        end do
      end do
    end if

    loss = loss / nValues

  end function rmsLoss


  !> Calculates the absolute minium deviation between predictions and targets.
  pure function minError(predicts, globalTargets, atomicTargets) result(err)

    !> neural network predictions
    type(TPredicts), intent(in) :: predicts

    !> system-wide target reference data
    type(TRealArray1D), intent(in) :: globalTargets(:)

    !> atomic target reference data
    type(TRealArray2D), intent(in) :: atomicTargets(:)

    !! temporary real valued storage for summed up system-wide predictions
    real(dp), allocatable :: globalPredicts(:)

    !! absolute minium deviation between predictions and targets
    real(dp) :: err

    !! temporary error storage
    real(dp) :: globalTargetsErr, atomicTargetsErr, tmpErr

    !! auxiliary variable
    integer :: iSys

    if (predicts%nGlobalTargets > 0) then
      globalPredicts = sum(predicts%sys(1)%array(1:predicts%nGlobalTargets, :), dim=2)
      globalTargetsErr = minval(abs(globalPredicts - globalTargets(1)%array))
      do iSys = 2, predicts%nDatapoints
        globalPredicts = sum(predicts%sys(iSys)%array(1:predicts%nGlobalTargets, :), dim=2)
        tmpErr = minval(abs(globalPredicts - globalTargets(iSys)%array))
        if (tmpErr < globalTargetsErr) globalTargetsErr = tmpErr
      end do
    end if

    if (predicts%nAtomicTargets > 0) then
      atomicTargetsErr = minval(abs(predicts%sys(1)%array(predicts%nGlobalTargets + 1:, :)&
          & - atomicTargets(1)%array))
      do iSys = 2, predicts%nDatapoints
        tmpErr = minval(abs(predicts%sys(iSys)%array(predicts%nGlobalTargets + 1:, :)&
            & - atomicTargets(iSys)%array))
        if (tmpErr < atomicTargetsErr) atomicTargetsErr = tmpErr
      end do
    end if

    if (predicts%nGlobalTargets > 0 .and. predicts%nAtomicTargets == 0) then
      err = globalTargetsErr
    elseif (predicts%nGlobalTargets == 0 .and. predicts%nAtomicTargets > 0) then
      err = atomicTargetsErr
    elseif (predicts%nGlobalTargets > 0 .and. predicts%nAtomicTargets > 0) then
      err = min(globalTargetsErr, atomicTargetsErr)
    else
      err = 0.0_dp
    end if

  end function minError


  !> Calculates the absolute maximum deviation between predictions and targets.
  pure function maxError(predicts, globalTargets, atomicTargets) result(err)

    !> neural network predictions
    type(TPredicts), intent(in) :: predicts

    !> system-wide target reference data
    type(TRealArray1D), intent(in) :: globalTargets(:)

    !> atomic target reference data
    type(TRealArray2D), intent(in) :: atomicTargets(:)

    !! temporary real valued storage for summed up system-wide predictions
    real(dp), allocatable :: globalPredicts(:)

    !! absolute maximum deviation between predictions and targets
    real(dp) :: err

    !! temporary error storage
    real(dp) :: globalTargetsErr, atomicTargetsErr, tmpErr

    !! auxiliary variable
    integer :: iSys

    if (predicts%nGlobalTargets > 0) then
      globalPredicts = sum(predicts%sys(1)%array(1:predicts%nGlobalTargets, :), dim=2)
      globalTargetsErr = maxval(abs(globalPredicts - globalTargets(1)%array))
      do iSys = 2, predicts%nDatapoints
        globalPredicts = sum(predicts%sys(iSys)%array(1:predicts%nGlobalTargets, :), dim=2)
        tmpErr = maxval(abs(globalPredicts - globalTargets(iSys)%array))
        if (tmpErr > globalTargetsErr) globalTargetsErr = tmpErr
      end do
    end if

    if (predicts%nAtomicTargets > 0) then
      atomicTargetsErr = maxval(abs(predicts%sys(1)%array(predicts%nGlobalTargets + 1:, :)&
          & - atomicTargets(1)%array))
      do iSys = 2, predicts%nDatapoints
        tmpErr = maxval(abs(predicts%sys(iSys)%array(predicts%nGlobalTargets + 1:, :)&
            & - atomicTargets(iSys)%array))
        if (tmpErr > atomicTargetsErr) atomicTargetsErr = tmpErr
      end do
    end if

    if (predicts%nGlobalTargets > 0 .and. predicts%nAtomicTargets == 0) then
      err = globalTargetsErr
    elseif (predicts%nGlobalTargets == 0 .and. predicts%nAtomicTargets > 0) then
      err = atomicTargetsErr
    elseif (predicts%nGlobalTargets > 0 .and. predicts%nAtomicTargets > 0) then
      err = max(globalTargetsErr, atomicTargetsErr)
    else
      err = 0.0_dp
    end if

  end function maxError

end module fnet_loss
