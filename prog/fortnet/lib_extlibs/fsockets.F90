!--------------------------------------------------------------------------------------------------!
!  FORTNET: A Behler-Parrinello-Neural-Network Implementation                                      !
!  Copyright (C) 2020 - 2024  T. W. van der Heide                                                  !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Defines a socket interface library.
module fnet_fsockets

  use fsockets
  implicit none


#:if WITH_SOCKETS

  !> True, if code was built with socket support
  logical, parameter :: withSockets = .true.

#:else

  !> True, if code was built with socket support
  logical, parameter :: withSockets = .false.

#:endif

end module fnet_fsockets
