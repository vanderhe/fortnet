!--------------------------------------------------------------------------------------------------!
!  FORTNET: A Behler-Parrinello-Neural-Network Implementation                                      !
!  Copyright (C) 2020 - 2021  T. W. van der Heide                                                  !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

module fnet_iterout

  ! use dftbp_xmlf90
  ! use dftbp_hsdutils, only : writeChildValue
  use dftbp_message, only : error
  use dftbp_accuracy, only: dp
  ! use dftbp_charmanip, only : i2c

  ! use fnet_nestedtypes, only : TRealArray2D, TPredicts

  implicit none

  private

  public :: writeIterTrajToFile


contains

  !> Write obtained results to fnetout.xml file
  subroutine writeIterTrajToFile(fname, loss, validLoss, gradients)

    !> filename (will be iterout.dat)
    character(len=*), intent(in) :: fname

    !> loss trajectory during training
    real(dp), intent(in), optional, target :: loss(:)

    !> validation loss trajectory during training
    real(dp), intent(in), optional, target :: validLoss(:)

    !> gradient trajectory during training
    real(dp), intent(in), optional, target :: gradients(:)

    !> pointer holding the available data
    real(dp), pointer :: pData(:,:)

    !> unique fileunit
    integer :: fd

    !> auxiliary variables
    integer :: iLine, nLines, nColumns

    if (present(loss)) then
      nLines = size(loss)
    elseif (present(validLoss)) then
      nLines = size(validLoss)
    elseif (present(gradients)) then
      nLines = size(gradients)
    else
      call error('Empty list of entries obtained. Provide at least one array with data.')
    end if

    nColumns = 0
    allocate(pData(nColumns, nLines))

    if (present(loss)) then
      nColumns = nColumns + 1
    end if

    if (present(validLoss)) then
      nColumns = nColumns + 1
    end if

    if (present(gradients)) then
      nColumns = nColumns + 1
    end if

    allocate(pData(nColumns, nLines))
    nColumns = 1

    if (present(loss)) then
      pData(nColumns, 1:nLines) = loss
      nColumns = nColumns + 1
    end if

    if (present(gradients)) then
      pData(nColumns, 1:nLines) = gradients
      nColumns = nColumns + 1
    end if

    if (present(validLoss)) then
      pData(nColumns, 1:nLines) = validLoss
    end if

    open(newunit=fd, file=fname, form='formatted', status='replace', action='write')

    do iLine = 1, nLines
      write(fd, '(I0,3ES26.16E3)') iLine, pData(:, iLine)
    end do

    close(fd)

  end subroutine writeIterTrajToFile

end module fnet_iterout
