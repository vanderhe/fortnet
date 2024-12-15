!--------------------------------------------------------------------------------------------------!
!  FORTNET: A Behler-Parrinello-Neural-Network Implementation                                      !
!  Copyright (C) 2020 - 2025  T. W. van der Heide                                                  !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

module fnet_workarounds

  implicit none

  private

  public :: myFindloc


contains

  ! workaround: GCC7 and GCC8 lack findloc as partially implemented F08 standard
  ! besides that, it seems to be buggy for some intel compilers as well
  pure function myFindloc(chars, selector) result(iPos)

    !> character array to be searched
    character(len=*), intent(in) :: chars(:)

    !> selector to search for
    character(len=*), intent(in) :: selector

    !! position in character array
    integer :: iPos

    iPos = 1

    do while (iPos <= size(chars))
      if (chars(iPos) == selector) then
        exit
      else
        iPos = iPos + 1
      end if
    end do

    if (iPos > size(chars)) then
      iPos = 0
    end if

  end function myFindloc

end module fnet_workarounds
