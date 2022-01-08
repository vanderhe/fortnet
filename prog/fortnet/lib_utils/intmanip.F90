!--------------------------------------------------------------------------------------------------!
!  FORTNET: A Behler-Parrinello-Neural-Network Implementation                                      !
!  Copyright (C) 2020 - 2022  T. W. van der Heide                                                  !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Provides several utilities to tackle integer based tasks.
module fnet_intmanip

  use dftbp_message, only : error

  implicit none

  private

  public :: getUniqueInt, getNumberOfUniqueInt
  public :: combinationsWithReplacement, factorial


contains

  !> Returns the unique integers of an array.
  pure subroutine getUniqueInt(array, unique)

    !> array to investigate
    integer, intent(in) :: array(:)

    !> number of unique entries
    integer, intent(out), allocatable :: unique(:)

    !> auxiliary variables
    integer :: ii, jj, kk, tmp(size(array))

    kk = 1
    tmp(1) = array(1)

    outer: do ii = 2, size(array)
      do jj = 1, kk
        if (tmp(jj) == array(ii)) then
          cycle outer
        end if
      end do
      kk = kk + 1
      tmp(kk) = array(ii)
    end do outer

    unique = tmp(1:kk)

  end subroutine getUniqueInt


  !> Returns the number of unique integers of an array.
  pure subroutine getNumberOfUniqueInt(array, nUnique)

    !> array to investigate
    integer, intent(in) :: array(:)

    !> number of unique entries
    integer, intent(out) :: nUnique

    !> auxiliary variables
    integer :: ii, jj, tmp(size(array))

    nUnique = 1
    tmp(1) = array(1)

    outer: do ii = 2, size(array)
      do jj = 1, nUnique
        if (tmp(jj) == array(ii)) then
          cycle outer
        end if
      end do
      nUnique = nUnique + 1
      tmp(nUnique) = array(ii)
    end do outer

  end subroutine getNumberOfUniqueInt


  !> Generates all combinations of multiple integers with replacement.
  subroutine combinationsWithReplacement(array, nChosen, comb)

    !> integer array to generate combinations for
    integer, intent(in) :: array(:)

    !> number of entries each tupel should contain
    integer, intent(in) :: nChosen

    !> container of all cominations obtained
    integer, intent(out), allocatable :: comb(:,:)

    !> pool length to draw numbers from
    integer :: nPool

    !> expected number of combinations
    integer :: nComb

    !> auxiliary variables
    integer :: ii, jj, ind
    integer, allocatable :: indices(:)

    nPool = size(array)

    nComb = factorial(nPool + nChosen - 1) / (factorial(nChosen) * factorial(nPool - 1))
    allocate(comb(nChosen, nComb))

    if (nPool <= 1) then
      call error('Cannot generate combinations for array of size equal or below 1.')
    end if

    if (nChosen <= 1) then
      call error('Cannot generate combinations for tupel size equal or below 1.')
    end if

    allocate(indices(nChosen))
    indices(:) = 1
    ind = 1

    do ii = 1, nChosen
      comb(ii, ind) = array(indices(ii))
    end do

    outer: do while (.true.)
      inner: do ii = nChosen, 1, -1
        if (indices(ii) /= nPool) then
          exit inner
        end if
      end do inner
      if (ii == 0) then
        return
      end if
      indices(ii:) = indices(ii) + 1
      ind = ind + 1
      do jj = 1, nChosen
        comb(jj, ind) = array(indices(jj))
      end do
    end do outer

  end subroutine combinationsWithReplacement


  !> Calculates the factorial of an integer.
  function factorial(xx) result(fact)

    !> integer to calculate the factorial for
    integer, intent(in) :: xx

    !> resulting factorial
    integer :: fact

    !> auxiliary variable
    integer :: ii

    if (xx < 0) then
      call error('Factorial is singular for negative integers.')
    end if

    fact = 1

    do ii = 2, xx
      fact = fact * ii
    end do

  end function factorial

end module fnet_intmanip
