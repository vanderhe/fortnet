!--------------------------------------------------------------------------------------------------!
!  FORTNET: A Behler-Parrinello-Neural-Network Implementation                                      !
!  Copyright (C) 2020 - 2021  T. W. van der Heide                                                  !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

module fnet_random

  use dftbp_message, only : error
  use dftbp_constants, only: pi
  use dftbp_accuracy, only: dp
  use dftbp_ranlux, only : TRanlux, getRandom

  implicit none

  private

  public :: randn
  public :: normalXavier

  interface normalXavier
    module procedure :: normalXavier1d, normalXavier2d
  end interface normalXavier

  interface randn
    module procedure :: randn1d, randn2d
  end interface randn


contains


  subroutine normalXavier1d(rndGen, array, nLastLayer, nCurrentLayer, gain, min)

    !> luxury pseudorandom generator instance
    type(TRanlux), intent(inout) :: rndGen

    !> contains random values at exit
    real(dp), intent(inout) :: array(:)

    !> number of neurons of (current-1) layer
    integer, intent(in) :: nLastLayer

    !> number of neurons of current layer
    integer, intent(in) :: nCurrentLayer

    !> optional scaling of distribution
    real(dp), intent(in), optional :: gain

    !> optional minimal value specification
    real(dp), intent(in), optional :: min

    !> variance of normal distribution
    real(dp) :: var

    !> scaling of distribution (default: 1.0)
    real(dp) :: scale

    !> minimal value specification
    real(dp) :: minval

    !> calculated boundary for random numbers to fulfill minimal value specification
    real(dp) :: bound

    !> luxury random numbers
    real(dp) :: rnd(size(array))

    if (present(gain)) then
      scale = gain
    else
      scale = 1.0_dp
    end if

    if (present(min)) then
      if (min < 1e-10_dp) then
        call error('Minimum value for normal distribution set too low.')
      end if
      minval = min
    else
      minval = 0.1_dp
    end if

    call getRandom(rndGen, rnd)

    var = scale**2 * 2.0_dp / real(nLastLayer + nCurrentLayer, dp)
    bound = sqrt(-2.0_dp * var * log(minval * sqrt(2.0_dp * pi * var)))

    rnd = (rnd - 0.5_dp) * 2.0_dp * bound

    array(:) = exp(-0.5_dp * rnd**2 / var) / sqrt(2.0_dp * pi * var)

  end subroutine normalXavier1d


  subroutine normalXavier2d(rndGen, array, nLastLayer, nCurrentLayer, gain, min)

    !> luxury pseudorandom generator instance
    type(TRanlux), intent(inout) :: rndGen

    !> contains random values at exit
    real(dp), intent(inout) :: array(:,:)

    !> number of neurons of (current-1) layer
    integer, intent(in) :: nLastLayer

    !> number of neurons of current layer
    integer, intent(in) :: nCurrentLayer

    !> optional scaling of distribution
    real(dp), intent(in), optional :: gain

    !> optional minimal value specification
    real(dp), intent(in), optional :: min

    !> variance of normal distribution
    real(dp) :: var

    !> scaling of distribution (default: 1.0)
    real(dp) :: scale

    !> minimal value specification
    real(dp) :: minval

    !> calculated boundary for random numbers to fulfill minimal value specification
    real(dp) :: bound

    !> luxury random numbers
    real(dp) :: rnd(size(array, dim=1), size(array, dim=2))

    if (present(gain)) then
      scale = gain
    else
      scale = 1.0_dp
    end if

    if (present(min)) then
      if (min < 1e-10_dp) then
        call error('Minimum value for normal distribution set too low.')
      end if
      minval = min
    else
      minval = 0.1_dp
    end if

    call getRandom(rndGen, rnd)

    var = scale**2 * 2.0_dp / real(nLastLayer + nCurrentLayer, dp)
    bound = sqrt(-2.0_dp * var * log(minval * sqrt(2.0_dp * pi * var)))

    rnd = (rnd - 0.5_dp) * 2.0_dp * bound

    array(:,:) = exp(-0.5_dp * rnd**2 / var) / sqrt(2.0_dp * pi * var)

  end subroutine normalXavier2d


  function randn1d(rndGen, amount) result(rnd1)

    !> luxury pseudorandom generator instance
    type(TRanlux), intent(inout) :: rndGen

    ! amount of random numbers to generate
    integer, intent(in) :: amount

    !> generated random numbers
    real(dp) :: rnd1(amount), rnd2(amount)

    call getRandom(rndGen, rnd1)
    call getRandom(rndGen, rnd2)

    rnd1 = sqrt(- 2.0_dp * log(rnd1)) * cos(2.0_dp * pi * rnd2)

  end function randn1d


  function randn2d(rndGen, amount1, amount2) result(rnd1)

    !> luxury pseudorandom generator instance
    type(TRanlux), intent(inout) :: rndGen

    ! amount of random numbers to generate
    integer, intent(in) :: amount1, amount2

    !> generated random numbers
    real(dp) :: rnd1(amount1, amount2), rnd2(amount1, amount2)

    call getRandom(rndGen, rnd1)
    call getRandom(rndGen, rnd2)

    rnd1 = sqrt(- 2.0_dp * log(rnd1)) * cos(2.0_dp * pi * rnd2)

  end function randn2d

end module fnet_random
