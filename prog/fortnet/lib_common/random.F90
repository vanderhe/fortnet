!--------------------------------------------------------------------------------------------------!
!  FORTNET: A Behler-Parrinello-Neural-Network Implementation                                      !
!  Copyright (C) 2020 - 2021  T. W. van der Heide                                                  !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Provides utilities to generate random distributions for network initialization.
module fnet_random

  use dftbp_message, only : error
  use dftbp_constants, only: pi
  use dftbp_accuracy, only: dp
  use dftbp_ranlux, only : TRanlux, getRandom

  implicit none

  private

  public :: knuthShuffle
  public :: normalXavier

  interface normalXavier
    module procedure :: normalXavier1d, normalXavier2d
  end interface normalXavier


contains

  !> Implements the Knuth-Shuffle algorithm.
  subroutine knuthShuffle(rndGen, array)

    !> luxury pseudorandom generator instance
    type(TRanlux), intent(inout) :: rndGen

    !> array that will be shuffled on exit
    integer, intent(inout) :: array(:)

    !> temporarily stores a real-valued random number
    real(dp) :: rnd

    !> auxiliary variables
    integer :: ii, rndpos, tmp

    do ii = size(array), 2, -1
      call getRandom(rndGen, rnd)
      rndpos = int(rnd * ii) + 1
      tmp = array(rndpos)
      array(rndpos) = array(ii)
      array(ii) = tmp
    end do

  end subroutine knuthShuffle


  !> Generates a 1d truncated, normal distributed Xavier network initialization.
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


  !> Generates a 2d truncated, normal distributed Xavier network initialization.
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

end module fnet_random
