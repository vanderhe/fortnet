!--------------------------------------------------------------------------------------------------!
!  FORTNET: A Behler-Parrinello-Neural-Network Implementation                                      !
!  Copyright (C) 2020 - 2021  T. W. van der Heide                                                  !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

module fnet_precond

  use dftbp_accuracy, only: dp
  use dftbp_message, only : error

  implicit none

  private

  public :: readPreconditioning, writePreconditioning


contains

  !> read target preconditioning informations from disk
  subroutine readPreconditioning(fname, tReadNetStats, tZscore, nTargets, zPrec)

    !> filename to read preconditioning parameters from
    character(len=*), intent(in) :: fname

    !> resume from network status on disc?
    logical, intent(in) :: tReadNetStats

    !> true, if z-score standardization should be applied
    logical, intent(in) :: tZscore

    !> number of target values per atom (if tAtomicTargets = .true.) or system
    integer, intent(out) :: nTargets

    !> storage container of means and variances to calculate z-score, shape: [nTargets, 2]
    real(dp), intent(out), allocatable :: zPrec(:,:)

    !> file and error identifier
    integer :: fp, iErr

    !> true, if preconditioning file is in place
    logical :: tExist

    !> auxiliary variable
    integer :: iTarget

    if (tReadNetStats .and. tZscore) then
      inquire(file=fname, exist=tExist)
    end if

    if (tReadNetStats .and. tZscore .and. (.not. tExist)) then
      call error("Preconditioning file '" // fname // "' is absent.")
    end if

    if (tReadNetStats .and. tZscore) then

      open(newunit=fp, file=fname, form='formatted', status='old', action='read', iostat=iErr)

      if (iErr /= 0) then
        call error("Could not open file '" // fname // "' for direct reading.")
      end if

      read(fp, *, iostat=iErr) nTargets

      allocate(zPrec(nTargets, 2))
      do iTarget = 1, nTargets
        read(fp, *, iostat=iErr) zPrec(iTarget, 1), zPrec(iTarget, 2)
      end do

      if (iErr /= 0) then
        call error("Error during reading file '" // fname // "'.")
      end if

      close(fp)
      
    end if

  end subroutine readPreconditioning


  !> write target preconditioning informations to disk
  subroutine writePreconditioning(fname, zPrec)

    !> filename or path to save preconditioning parameters to
    character(len=*), intent(in) :: fname

    !> storage container of means and variances to calculate z-score, shape: [nTargets, 2]
    real(dp), intent(in) :: zPrec(:,:)

    !> file and error identifier
    integer :: fp, iErr

    !> auxiliary variable
    integer :: iTarget

    open(newunit=fp, file=fname, form='formatted', status='replace', action='write',&
        & iostat=iErr)

    write(fp, '(I26)', iostat=iErr) size(zPrec, dim=1)

    do iTarget = 1, size(zPrec, dim=1)
      write(fp, '(2ES26.16E3)', iostat=iErr) zPrec(iTarget, 1), zPrec(iTarget, 2)
    end do

    if (iErr /= 0) then
      call error("Error during writing file '" // fname // "'.")
    end if

    close(fp)

  end subroutine writePreconditioning

end module fnet_precond
