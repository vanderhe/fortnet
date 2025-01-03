!--------------------------------------------------------------------------------------------------!
!  FORTNET: A Behler-Parrinello-Neural-Network Implementation                                      !
!  Copyright (C) 2020 - 2025  T. W. van der Heide                                                  !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Exporting mpifx functionality if compiled with mpi support, otherwise empty.
module dftbp_mpifx

#:if WITH_MPI
  use libmpifx_module
#:endif

  implicit none
  public

end module dftbp_mpifx
