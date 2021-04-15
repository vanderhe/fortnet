!--------------------------------------------------------------------------------------------------!
!  FORTNET: A Behler-Parrinello-Neural-Network Implementation                                      !
!  Copyright (C) 2020 - 2021  T. W. van der Heide                                                  !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

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


  pure function gaussian(xx) result(res)

    !> array to calculate the transfer function for
    real(dp), intent(in) :: xx(:)

    !> corresponding transfer function values
    real(dp) :: res(size(xx))

    res = exp(- xx**2)

  end function gaussian


  pure function gaussianDeriv(xx) result(res)

    !> array to calculate the transfer function for
    real(dp), intent(in) :: xx(:)

    !> corresponding transfer function values
    real(dp) :: res(size(xx))

    res = - 2.0_dp * xx * gaussian(xx)

  end function gaussianDeriv


  pure function relu(xx) result(res)

    !> array to calculate the transfer function for
    real(dp), intent(in) :: xx(:)

    !> corresponding transfer function values
    real(dp) :: res(size(xx))

    res = max(0.0_dp, xx)

  end function relu


  pure function reluDeriv(xx) result(res)

    !> array to calculate the transfer function for
    real(dp), intent(in) :: xx(:)

    !> corresponding transfer function values
    real(dp) :: res(size(xx))

    where (xx > 0.0_dp)
      res = 1.0_dp
    elsewhere
      res = 0.0_dp
    end where

  end function reluDeriv


  pure function sigmoid(xx) result(res)

    !> array to calculate the transfer function for
    real(dp), intent(in) :: xx(:)

    !> corresponding transfer function values
    real(dp) :: res(size(xx))

    res = 1.0_dp / (1.0_dp + exp(- xx))

  endfunction sigmoid


  pure function sigmoidDeriv(xx) result(res)

    !> array to calculate the transfer function for
    real(dp), intent(in) :: xx(:)

    !> corresponding transfer function values
    real(dp) :: res(size(xx))

    res = sigmoid(xx) * (1.0_dp - sigmoid(xx))

  end function sigmoidDeriv


  pure function heaviside(xx) result(res)

    !> array to calculate the transfer function for
    real(dp), intent(in) :: xx(:)

    !> corresponding transfer function values
    real(dp) :: res(size(xx))

    where (xx > 0.0_dp)
      res = 1.0_dp
    elsewhere
      res = 0.0_dp
    end where

  end function heaviside


  pure function heavisideDeriv(xx) result(res)

    !> array to calculate the transfer function for
    real(dp), intent(in) :: xx(:)

    !> corresponding transfer function values
    real(dp) :: res(size(xx))

    res = 0.0_dp

  end function heavisideDeriv


  pure function tanhf(xx) result(res)

    !> array to calculate the transfer function for
    real(dp), intent(in) :: xx(:)

    !> corresponding transfer function values
    real(dp) :: res(size(xx))

    res = tanh(xx)

  end function tanhf


  pure function tanhDeriv(xx) result(res)

    !> array to calculate the transfer function for
    real(dp), intent(in) :: xx(:)

    !> corresponding transfer function values
    real(dp) :: res(size(xx))

    res = 1.0_dp - tanh(xx)**2

  end function tanhDeriv


  pure function linear(xx) result(res)

    !> array to calculate the transfer function for
    real(dp), intent(in) :: xx(:)

    !> corresponding transfer function values
    real(dp) :: res(size(xx))

    res = xx

  end function linear


  pure function linearDeriv(xx) result(res)

    !> array to calculate the transfer function for
    real(dp), intent(in) :: xx(:)

    !> corresponding transfer function values
    real(dp) :: res(size(xx))

    res = 1.0_dp

  end function linearDeriv

end module fnet_transfer
