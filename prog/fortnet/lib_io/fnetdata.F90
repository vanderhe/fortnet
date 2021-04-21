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

  use fnet_nestedtypes, only : TRealArray2D

  implicit none

  private

  public :: readFnetdataGeometry, readContiguousFnetdataGeometries
  public :: readFnetdataTargets, readContiguousFnetdataTargets


contains

  !> interpret the geometry information stored in fnetdata.xml files
  subroutine readFnetdataGeometry(root, geo)

    !> pointer to the node, containing the data
    type(fnode), pointer :: root

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
    call getChild(root, 'geometry', geonode)
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


  !> interpret the geometry information stored in a contiguous fnetdata.xml file
  subroutine readContiguousFnetdataGeometries(root, geo)

    !> pointer to the node, containing the data
    type(fnode), pointer :: root

    !> contains the geometry information on exit
    type(TGeometry), intent(out), allocatable :: geo(:)

    !> node to get dataset size from
    type(fnode), pointer :: datasetnode

    !> temporary pointer to node, containing information
    type(fnode), pointer :: tmp

    !> number of datapoints contained in the dataset file
    integer :: nDatapoints

    !> auxiliary variable
    integer :: iGeo

    ! read dataset size
    call getChild(root, 'dataset', datasetnode)
    call getChildValue(datasetnode, 'ndatapoints', nDatapoints)

    allocate(geo(nDatapoints))

    do iGeo = 1, nDatapoints
      call getChild(root, 'datapoint' // i2c(iGeo), tmp)
      call readFnetdataGeometry(tmp, geo(iGeo))
    end do

  end subroutine readContiguousFnetdataGeometries


  !> interpret the target information of a training xml-tree
  subroutine readTargets(root, nTargets, nAtom, tAtomicTargets, targets)

    !> pointer to the node, containing the data
    type(fnode), pointer :: root

    !> number of targets per atom or system
    integer, intent(in) :: nTargets

    !> number of atoms of the corresponding geometry
    integer, intent(in) :: nAtom

    !> true, if targets for every atom provided
    logical, intent(in) :: tAtomicTargets

    !> contains the target information on exit
    real(dp), intent(out), allocatable :: targets(:,:)

    !> node containing target informations
    type(fnode), pointer :: child

    !> temporary target storage container
    real(dp), allocatable :: tmpTargets(:)

    if (tAtomicTargets) then
      allocate(tmpTargets(nTargets * nAtom))
    else
      allocate(tmpTargets(nTargets))
    end if

    call getChildValue(root, 'targets', tmpTargets, child=child)

    if (tAtomicTargets) then
      targets = reshape(tmpTargets, [nTargets, nAtom])
    else
      targets = reshape(tmpTargets, [nTargets, 1])
    end if

  end subroutine readTargets


  !> interpret the target information stored in fnetdata.xml files
  subroutine readFnetdataTargets(fnetdata, nAtom, targets, tAtomicTargets)

    !> pointer to the node, containing the data
    type(fnode), pointer :: fnetdata

    !> number of atoms of current geometry
    integer, intent(in) :: nAtom

    !> contains the target information on exit
    real(dp), intent(out), allocatable :: targets(:,:)

    !> true, if targets for every atom provided
    logical, intent(out) :: tAtomicTargets

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

    call readTargets(trainnode, nTargets, nAtom, tAtomicTargets, targets)

  end subroutine readFnetdataTargets


  !> interpret the geometry information stored in a contiguous fnetdata.xml file
  subroutine readContiguousFnetdataTargets(root, geo, targets, nTargets, tAtomicTargets)

    !> pointer to the node, containing the data
    type(fnode), pointer :: root

    !> contains corresponding geometry information
    type(TGeometry), intent(in) :: geo(:)

    !> target values for training
    type(TRealArray2D), intent(out), allocatable :: targets(:)

    !> number of target values per atom (if tAtomicTargets = .true.) or system
    integer, intent(out) :: nTargets

    !> true, if targets for every atom provided
    logical, intent(out) :: tAtomicTargets

    !> node to get target properties from
    type(fnode), pointer :: trainnode

    !> temporary pointer to node, containing information
    type(fnode), pointer :: tmp

    !> auxiliary variable
    integer :: iTarget

    ! read target information
    call getChild(root, 'training', trainnode)
    call getChildValue(trainnode, 'atomic', tAtomicTargets)
    call getChildValue(trainnode, 'ntargets', nTargets)

    allocate(targets(size(geo)))

    do iTarget = 1, size(geo)
      call getChild(root, 'datapoint' // i2c(iTarget), tmp)
      call getChild(tmp, 'training', trainnode)
      call readTargets(trainnode, nTargets, geo(iTarget)%nAtom, tAtomicTargets,&
          & targets(iTarget)%array)
    end do

  end subroutine readContiguousFnetdataTargets

end module fnet_fnetdata
