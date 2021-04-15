!--------------------------------------------------------------------------------------------------!
!  FORTNET: A Behler-Parrinello-Neural-Network Implementation                                      !
!  Copyright (C) 2020 - 2021  T. W. van der Heide                                                  !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

module fnet_fnetdata

  use xmlf90_strings
  use xmlf90_flib_dom, only : fnode

  use dftbp_linkedlist
  use dftbp_accuracy, only: dp
  use dftbp_typegeometry, only : TGeometry, normalize
  use dftbp_typegeometryhsd, only : readTGeometryHSD
  use dftbp_unitconversion, only : lengthUnits, unit
  use dftbp_hsdutils, only : getChild, getChildValue, setChildValue, detailedError
  use dftbp_hsdutils2, only : getModifierIndex
  use dftbp_simplealgebra, only : invert33, determinant33
  use dftbp_charmanip, only : tolower, i2c

  implicit none

  private

  public :: readFnetdataGeometry, readFnetdataTargets


contains

  !> interpret the geometry information stored in fnetdata.xml files
  subroutine readFnetdataGeometry(fnetdata, geo)

    !> pointer to the node, containing the data
    type(fnode), pointer :: fnetdata

    !> contains the geometry information on exit
    type(TGeometry), intent(out) :: geo

    !> node to get geometry from
    type(fnode), pointer :: geonode

    !> nodes, containing the informations
    type(fnode), pointer :: child, typesAndCoords

    !> temporary string buffer
    type(TListString) :: strBuffer

    !> temporary integer and real valued buffer
    type(TListIntR1) :: intBuffer
    type(TListRealR1) :: realBuffer

    !> lattice vectors and determinant
    real(dp) :: latvec(9), det

    !> auxiliary variable
    integer, allocatable :: tmpInt(:,:)

    ! read geometry information
    call getChild(fnetdata, 'geometry', geonode)
    call getChildValue(geonode, 'periodic', geo%tPeriodic, default=.false.)

    ! no support for helical boundary conditions
    geo%tHelical = .false.

    call init(strBuffer)

    call getChildValue(geonode, 'TypeNames', strBuffer)
    geo%nSpecies = len(strBuffer)

    if (geo%nSpecies == 0) then
      call detailedError(geonode, 'Missing species names.')
    end if

    allocate(geo%speciesNames(geo%nSpecies))

    call asArray(strBuffer, geo%speciesNames)
    call destruct(strBuffer)

    call getChildValue(geonode, 'fractional', geo%tFracCoord)

    if (geo%tFracCoord .and. (.not. geo%tPeriodic)) then
      call detailedError(typesAndCoords, 'Fractional coordinates are only allowed for periodic&
          & systems.')
    end if

    call init(intBuffer)
    call init(realBuffer)

    call getChildValue(geonode, 'TypesAndCoordinates', 1, intBuffer, 3, realBuffer,&
        & child=typesAndCoords)

    geo%nAtom = len(intBuffer)

    if (geo%nAtom == 0) then
      call detailedError(typesAndCoords, 'Missing coordinates')
    end if

    allocate(geo%species(geo%nAtom))
    allocate(geo%coords(3, geo%nAtom))
    allocate(tmpInt(1, geo%nAtom))

    call asArray(intBuffer, tmpInt)
    call destruct(intBuffer)

    geo%species(:) = tmpInt(1,:)
    deallocate(tmpInt)

    !! Check validity of species
    if (any(geo%species < 1 .or. geo%species > geo%nSpecies)) then
      call detailedError(typesAndCoords, 'Type index must be between 1 and ' // i2c(geo%nSpecies)&
          & // '.')
    end if

    call asArray(realBuffer, geo%coords)
    call destruct(realBuffer)

    if (geo%tPeriodic) then
      allocate(geo%origin(3))
      call getChildValue(geonode, 'CoordinateOrigin', geo%origin, [0.0_dp, 0.0_dp, 0.0_dp])
      geo%coords(:,:) = geo%coords - spread(geo%origin, 2, geo%nAtom)
      allocate(geo%latVecs(3,3))
      call getChildValue(geonode, 'LatticeVectors', latvec, child=child)
      geo%latVecs(:,:) = reshape(latvec, [3, 3])
      if (geo%tFracCoord) then
        geo%coords(:,:) = matmul(geo%latVecs, geo%coords)
        geo%origin(:) = matmul(geo%latVecs, geo%origin)
      end if
      allocate(geo%recVecs2p(3, 3))
      det = determinant33(geo%latVecs)
      if (abs(det) < 1e-12_dp) then
        call detailedError(child, 'Dependent lattice vectors.')
      end if
      call invert33(geo%recVecs2p, geo%latVecs, det)
      geo%recVecs2p(:,:) = reshape(geo%recVecs2p, [3, 3], order=[2, 1])
    end if

    call normalize(geo)

  end subroutine readFnetdataGeometry


  !> interpret the target information stored in fnetdata.xml files
  subroutine readFnetdataTargets(fnetdata, nAtom, tAtomicTargets, targets)

    !> pointer to the node, containing the data
    type(fnode), pointer :: fnetdata

    !> number of atoms of current geometry
    integer, intent(in) :: nAtom

    !> true, if targets for every atom provided
    logical, intent(out) :: tAtomicTargets

    !> contains the target information on exit
    real(dp), intent(out), allocatable :: targets(:,:)

    !> node to get targets from
    type(fnode), pointer :: trainnode

    !> node containing target informations
    type(fnode), pointer :: child

    !> temporary string buffer
    type(string) :: modifier

    !> temporary target storage container
    real(dp), allocatable :: tmpTargets(:)

    !> number of targets per atom or system
    integer :: nTargets

    ! read target information
    call getChild(fnetdata, 'training', trainnode)
    call getChildValue(trainnode, 'atomic', tAtomicTargets)
    call getChildValue(trainnode, 'ntargets', nTargets)

    if (tAtomicTargets) then
      allocate(tmpTargets(nTargets * nAtom))
    else
      allocate(tmpTargets(nTargets))
    end if

    call getChildValue(trainnode, 'targets', tmpTargets, modifier=modifier, child=child)

    if (tAtomicTargets) then
      targets = reshape(tmpTargets, [nTargets, nAtom])
    else
      targets = reshape(tmpTargets, [nTargets, 1])
    end if

  end subroutine readFnetdataTargets

end module fnet_fnetdata
