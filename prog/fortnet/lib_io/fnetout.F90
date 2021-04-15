!--------------------------------------------------------------------------------------------------!
!  FORTNET: A Behler-Parrinello-Neural-Network Implementation                                      !
!  Copyright (C) 2020 - 2021  T. W. van der Heide                                                  !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

module fnet_fnetout

  use dftbp_xmlf90
  use dftbp_hsdutils, only : writeChildValue
  use dftbp_message, only : error
  use dftbp_accuracy, only: dp
  use dftbp_charmanip, only : i2c

  use fnet_nestedtypes, only : TRealArray2D, TPredicts

  implicit none

  private

  public :: writeFnetout


contains

  !> Write obtained results to fnetout.xml file
  subroutine writeFnetout(fname, mode, targets, output, tAtomicTargets)

    !> filename (will be fnetout.xml)
    character(len=*), intent(in) :: fname

    !> mode of calculation (train, validate, run)
    character(len=*), intent(in) :: mode

    !> target values for training
    type(TRealArray2D), intent(in) :: targets(:)

    !> obtained output values of network
    type(TPredicts), intent(in) :: output

    !> true, if targets are atomic properties
    logical, intent(in) :: tAtomicTargets

    !> xml file instance
    type(xmlf_t) :: xf

    !> real valued onedimensional buffer
    real(dp), allocatable :: bufferRealR1(:)

    !> number of total datapoints/structures
    integer :: nDatapoints

    !> number of target values per atom (if tAtomicTargets = .true.) or system
    integer :: nTargets

    !> number of predictions per atom (if tAtomicTargets = .true.) or system
    integer :: nOutputs

    !> auxiliary variables
    integer :: iSys, iAtom

    if (mode /= 'validate' .and. mode /= 'predict') then
      call error('Invalid program running mode selected.')
    end if

    nDatapoints = size(output%sys)
    nOutputs = size(output%sys(1)%array, dim=1)

    select case(mode)
    case('validate')
      nTargets = size(targets(1)%array, dim=1)
    case('predict')
      nTargets = 0
    end select

    call xml_OpenFile(fname, xf, indent=.true.)
    call xml_ADDXMLDeclaration(xf)

    call xml_NewElement(xf, 'fnetout')

    call writeChildValue(xf, 'mode', [mode])

    call xml_NewElement(xf, 'output')

    call writeChildValue(xf, 'ndatapoints', nDatapoints)
    call writeChildValue(xf, 'atomic', tAtomicTargets)

    select case(mode)
    case('validate')
      call writeChildValue(xf, 'ntargets', nTargets)
    case('predict')
      call writeChildValue(xf, 'ntargets', nOutputs)
    end select

    allocate(bufferRealR1(nOutputs + nTargets))

    do iSys = 1, nDatapoints

      select case(mode)
      case('validate')
        if (tAtomicTargets) then
          call xml_NewElement(xf, 'datapoint' // trim(i2c(iSys)))
          do iAtom = 1, size(output%sys(iSys)%array, dim=2)
            bufferRealR1(:) = [output%sys(iSys)%array(:, iAtom), targets(iSys)%array(:, iAtom)]
            call writeChildValue(xf, 'atom' // trim(i2c(iAtom)), bufferRealR1)
          end do
          call xml_EndElement(xf, 'datapoint' // trim(i2c(iSys)))
        else
          bufferRealR1(:) = [output%sys(iSys)%array(:, 1), targets(iSys)%array(:, 1)]
          call writeChildValue(xf, 'datapoint' // trim(i2c(iSys)), bufferRealR1)
        end if
      case('predict')
        if (tAtomicTargets) then
          call xml_NewElement(xf, 'datapoint' // trim(i2c(iSys)))
          do iAtom = 1, size(output%sys(iSys)%array, dim=2)
            bufferRealR1(:) = output%sys(iSys)%array(:, iAtom)
            call writeChildValue(xf, 'atom' // trim(i2c(iAtom)), bufferRealR1)
          end do
          call xml_EndElement(xf, 'datapoint' // trim(i2c(iSys)))
        else
          bufferRealR1(:) = output%sys(iSys)%array(:, 1)
          call writeChildValue(xf, 'datapoint' // trim(i2c(iSys)), bufferRealR1)
        end if
      end select

    end do

    call xml_EndElement(xf, 'output')
    call xml_EndElement(xf, 'fnetout')
    call xml_Close(xf)

  end subroutine writeFnetout

end module fnet_fnetout
