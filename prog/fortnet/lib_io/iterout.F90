!--------------------------------------------------------------------------------------------------!
!  FORTNET: A Behler-Parrinello-Neural-Network Implementation                                      !
!  Copyright (C) 2020 - 2025  T. W. van der Heide                                                  !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Provides utilities to write out more detailed information about the loss course.
module fnet_iterout

  use dftbp_accuracy, only: dp
  use dftbp_message, only : error
  use dftbp_charmanip, only : i2c

  implicit none

  private

  public :: writeIterTrajToFile


contains

  !> Writes obtained loss/v-loss/gradients to file.
  subroutine writeIterTrajToFile(fname, trainLoss, validLoss, gradients)

    !> filename (will be iterout.dat)
    character(len=*), intent(in) :: fname

    !> loss trajectory during training
    real(dp), intent(in), optional, target :: trainLoss(:)

    !> validation loss trajectory during training
    real(dp), intent(in), optional, target :: validLoss(:)

    !> gradient trajectory during training
    real(dp), intent(in), optional, target :: gradients(:)

    !! pointer holding the available data
    real(dp), pointer :: pData(:,:)

    !! formatting of output file
    character(len=:), allocatable :: fmt

    !! unique fileunit
    integer :: fd

    !! auxiliary variables
    integer :: iLine, nLines, nColumns

    if (present(trainLoss)) then
      nLines = size(trainLoss)
    elseif (present(validLoss)) then
      nLines = size(validLoss)
    elseif (present(gradients)) then
      nLines = size(gradients)
    else
      call error('Empty list of entries obtained. Provide at least one array with data.')
    end if

    fmt = '(I' // i2c(len(i2c(nLines))) // ',3ES26.16E3)'

    nColumns = 0
    allocate(pData(nColumns, nLines))

    if (present(trainLoss)) then
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

    if (present(trainLoss)) then
      pData(nColumns, 1:nLines) = trainLoss
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
      write(fd, fmt) iLine, pData(:, iLine)
    end do

    close(fd)

  end subroutine writeIterTrajToFile

end module fnet_iterout
