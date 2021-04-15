!--------------------------------------------------------------------------------------------------!
!  FORTNET: A Behler-Parrinello-Neural-Network Implementation                                      !
!  Copyright (C) 2020 - 2021  T. W. van der Heide                                                  !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

module fnet_parallel

  use dftbp_accuracy, only: dp

  implicit none

  private

  public :: getStartAndEndIndex


contains

  pure subroutine getStartAndEndIndex(nSystems, nProcs, iProc, iStart, iEnd)

    !> array size to split
    integer, intent(in) :: nSystems

    !> number of available processes
    integer, intent(in) :: nProcs

    !> current process index
    integer, intent(in) :: iProc

    !> start and end index of current tile
    integer, intent(out) :: iStart, iEnd

    !> size of splitted index regions
    integer :: splitSize

    !> number of systems that exceeds integer times nProcs
    integer :: offset

    splitSize = nSystems / nProcs

    ! start and end indices assuming equal split sizes
    iStart = iProc * splitSize + 1
    iEnd = iStart + splitSize - 1

    ! distribute possible remainder to the tiles at the end
    offset = nProcs - mod(nSystems, nProcs)
    if (iProc + 1 > offset) then
      iStart = iStart + iProc - offset
      iEnd = iEnd + iProc - offset + 1
    end if

  end subroutine getStartAndEndIndex

end module fnet_parallel
