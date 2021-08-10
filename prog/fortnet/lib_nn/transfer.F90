!--------------------------------------------------------------------------------------------------!
!  FORTNET: A Behler-Parrinello-Neural-Network Implementation                                      !
!  Copyright (C) 2020 - 2021  T. W. van der Heide                                                  !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Provides various activation/transfer functions as well as their derivatives.
module fnet_transfer

  use dftbp_accuracy, only: dp

  implicit none

  private

  public :: transferFunc
  public :: gaussian, gaussianDeriv
  public :: relu, reluDeriv
  public :: sigmoid, sigmoidDeriv
  public :: heaviside, heavisideDeriv
  public :: tanhf, tanhDeriv
  public :: linear, linearDeriv

  interface

    pure function transferFunc(xx)

      use dftbp_accuracy, only: dp

      implicit none

      !> array to calculate the transfer function for
      real(dp), intent(in) :: xx(:)

      !> corresponding transfer function values
      real(dp) :: transferFunc(size(xx))

    end function transferFunc

  end interface


contains


  !> Calculates a gaussian transfer for given arguments.
  pure function gaussian(xx) result(res)

    !> array to calculate the transfer function for
    real(dp), intent(in) :: xx(:)

    !> corresponding transfer function values
    real(dp) :: res(size(xx))

    res(:) = exp(- xx**2)

  end function gaussian


  !> Calculates the derivatives of a gaussian transfer.
  pure function gaussianDeriv(xx) result(res)

    !> array to calculate the transfer function for
    real(dp), intent(in) :: xx(:)

    !> corresponding transfer function values
    real(dp) :: res(size(xx))

    res(:) = - 2.0_dp * xx * gaussian(xx)

  end function gaussianDeriv


  !> Calculates a ReLU transfer for given arguments.
  pure function relu(xx) result(res)

    !> array to calculate the transfer function for
    real(dp), intent(in) :: xx(:)

    !> corresponding transfer function values
    real(dp) :: res(size(xx))

    res(:) = max(0.0_dp, xx)

  end function relu


  !> Calculates the derivatives of a ReLU transfer.
  pure function reluDeriv(xx) result(res)

    !> array to calculate the transfer function for
    real(dp), intent(in) :: xx(:)

    !> corresponding transfer function values
    real(dp) :: res(size(xx))

    where (xx > 0.0_dp)
      res(:) = 1.0_dp
    elsewhere
      res(:) = 0.0_dp
    end where

  end function reluDeriv


  !> Calculates a logistic, sigmoidal transfer for given arguments.
  pure function sigmoid(xx) result(res)

    !> array to calculate the transfer function for
    real(dp), intent(in) :: xx(:)

    !> corresponding transfer function values
    real(dp) :: res(size(xx))

    res(:) = 1.0_dp / (1.0_dp + exp(- xx))

  endfunction sigmoid


  !> Calculates the derivatives of a logistic, sigmoidal transfer.
  pure function sigmoidDeriv(xx) result(res)

    !> array to calculate the transfer function for
    real(dp), intent(in) :: xx(:)

    !> corresponding transfer function values
    real(dp) :: res(size(xx))

    res(:) = sigmoid(xx) * (1.0_dp - sigmoid(xx))

  end function sigmoidDeriv


  !> Calculates a heaviside transfer for given arguments.
  pure function heaviside(xx) result(res)

    !> array to calculate the transfer function for
    real(dp), intent(in) :: xx(:)

    !> corresponding transfer function values
    real(dp) :: res(size(xx))

    where (xx > 0.0_dp)
      res(:) = 1.0_dp
    elsewhere
      res(:) = 0.0_dp
    end where

  end function heaviside


  !> Calculates the derivatives of a heaviside transfer.
  pure function heavisideDeriv(xx) result(res)

    !> array to calculate the transfer function for
    real(dp), intent(in) :: xx(:)

    !> corresponding transfer function values
    real(dp) :: res(size(xx))

    res(:) = 0.0_dp

  end function heavisideDeriv


  !> Calculates a hyperbolic tangent, sigmoidal transfer for given arguments.
  pure function tanhf(xx) result(res)

    !> array to calculate the transfer function for
    real(dp), intent(in) :: xx(:)

    !> corresponding transfer function values
    real(dp) :: res(size(xx))

    res(:) = tanh(xx)

  end function tanhf


  !> Calculates the derivatives of a hyperbolic tangent, sigmoidal transfer.
  pure function tanhDeriv(xx) result(res)

    !> array to calculate the transfer function for
    real(dp), intent(in) :: xx(:)

    !> corresponding transfer function values
    real(dp) :: res(size(xx))

    res(:) = 1.0_dp - tanh(xx)**2

  end function tanhDeriv


  !> Calculates a linear transfer for given arguments.
  pure function linear(xx) result(res)

    !> array to calculate the transfer function for
    real(dp), intent(in) :: xx(:)

    !> corresponding transfer function values
    real(dp) :: res(size(xx))

    res(:) = xx

  end function linear


  !> Calculates the derivatives of a linear transfer.
  pure function linearDeriv(xx) result(res)

    !> array to calculate the transfer function for
    real(dp), intent(in) :: xx(:)

    !> corresponding transfer function values
    real(dp) :: res(size(xx))

    res(:) = 1.0_dp

  end function linearDeriv

end module fnet_transfer
