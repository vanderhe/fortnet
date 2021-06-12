!--------------------------------------------------------------------------------------------------!
!  FORTNET: A Behler-Parrinello-Neural-Network Implementation                                      !
!  Copyright (C) 2020 - 2021  T. W. van der Heide                                                  !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains the routines for initialising Fortnet from HSD input.
module fnet_initprogram

  use dftbp_xmlf90
  use dftbp_linkedlist
  use dftbp_assert
  use dftbp_hsdutils, only : getChild, setChild, getChildValue, setChildValue
  use dftbp_hsdutils, only : convRangeToInt, detailedError
  use dftbp_hsdutils2, only : readHSDAsXML, setUnprocessed, warnUnprocessedNodes
  use dftbp_charmanip, only : unquote, tolower, i2c
  use dftbp_accuracy, only: dp, lc
  use dftbp_message, only : error, warning
  use dftbp_steepdesc, only : init
  use dftbp_conjgrad, only : init
  use dftbp_lbfgs, only : TLbfgs_init
  use dftbp_fire, only : TFire_init
  use dftbp_ranlux, only : TRanlux, init
  use dftbp_constants, only : AA__Bohr
  use dftbp_typegeometry, only : TGeometry
  use dftbp_globalenv, only : stdOut
  use dftbp_hsdparser, only : parseHSD, dumpHSD

  use fnet_precond, only : readPreconditioning, writePreconditioning
  use fnet_optimizers, only : TOptimizer, optimizerTypes, init, reset
  use fnet_fnetdata, only : readFnetdataWeight, readContiguousFnetdataWeights
  use fnet_fnetdata, only : readFnetdataGeometry, readContiguousFnetdataGeometries
  use fnet_fnetdata, only : readFnetdataTargets, readContiguousFnetdataTargets
  use fnet_fnetdata, only : readFnetdataAtomIdentifier, readContiguousFnetdataAtomIdentifier
  use fnet_fnetdata, only : readFnetdataFeatures, readContiguousFnetdataFeatures
  use fnet_fnetdata, only : inquireFeatures
  use fnet_acsf, only : TAcsf
  use fnet_nestedtypes, only : TIntArray1D, TRealArray1D, TRealArray2D, TEnv
  use fnet_nestedtypes, only : TWrapSteepDesc, TWrapConjGrad, TWrapLbfgs, TWrapFire
  use fnet_workarounds, only : myFindloc
  use fnet_loss, only : lossFunc, maLoss, msLoss, mslLoss, rmsLoss

#:if WITH_MPI
  use fnet_mpifx
#:endif

  implicit none

  private
  save

  public :: TProgramVariables, TProgramVariables_init
  public :: TFeatures, TFeatures_init, TFeatures_collect
  public :: TArch, TData, TEnv, TExternal, TOption
  public :: initOptimizer, readAcsfFromFile


  !> Data type containing variables of the Network block
  type TArch

    !> number of nodes per hidden layer, expected shape: [nHiddenLayer]
    integer, allocatable :: hidden(:)

    !> number of nodes per layer, including in- and output, expected shape: [nHiddenLayer + 2]
    integer, allocatable :: allDims(:)

    !> number of hidden layers
    integer :: nHiddenLayer

    !> number of network paramaters (weights + biases) per sub-nn
    integer :: nSubNnParams

    !> type of activation functions
    character(len=:), allocatable :: activation

    !> architecture type (currently, only the BPNN topology is available)
    character(len=:), allocatable :: type

  end type TArch


  !> Data type containing variables of the Mapping block
  type TMapping

    !> number of radial symmetry functions
    integer :: nRadial

    !> number of angular symmetry functions
    integer :: nAngular

    !> cutoff radius, defining the sphere to search for neighboring atoms
    real(dp) :: rCut

    !> true, if standardization of input features is desired
    logical :: tStandardize

    !> structural mapping type (currently only ACSF's are available)
    character(len=:), allocatable :: type

  end type TMapping


  !> Data type containing the collected features
  type TFeatures

    !> pointer to the different input features for training
    type(TRealArray2D), allocatable :: features(:)

    !> pointer to the different input features for validation
    type(TRealArray2D), allocatable :: validFeatures(:)

#:if WITH_MPI
  contains
    procedure :: sync => TFeatures_sync
#:endif

  end type TFeatures


  !> Data type containing variables of the Training block
  type TTraining

    !> integer ID of specified optimizert
    integer :: iOptimizer

    !> general function optimizer
    type(TOptimizer), allocatable :: pOptimizer(:)

    !> maximum number of training iterations
    integer :: nTrainIt

    !> printout loss/gradient information every nPrintOut steps
    integer :: nPrintOut

    !> save network status every nSaveNet steps
    integer :: nSaveNet

    !> type of loss function to use during the training
    character(len=:), allocatable :: lossType

    !> procedure, pointing to the choosen loss function
    procedure(lossFunc), pointer, nopass :: loss => null()

    !> gradient threshold where to stop the training, if provided
    real(dp) :: threshold

    !> learning rate for steepest descent optimizer
    real(dp) :: learningRate

    !> minimal displacement in at least one component in one step
    real(dp) :: minDisplacement

    !> maximum displacement in at least one component in one step
    real(dp) :: maxDisplacement

    !> number of past iterations which will be kept in memory by L-BFGS
    integer :: mem

    !> wether the maximum step size is considered for the quasi-Newton direction
    logical :: maxForQNDisplacement

    !> wether a line search should be used along the quasi-Newton direction
    logical :: tLinesearch

  contains

    procedure :: setLossFunc => TTraining_setLossFunc

  end type TTraining


  !> Data type containing variables of the Data block
  type TData

    !> prefix of paths to saved species networks
    character(len=1024) :: prefix

    !> suffix of paths to saved species networks
    character(len=1024) :: suffix

    !> path to file containing paths to training datapoints
    character(len=:), allocatable :: datapath

    !> path to file containing paths to validation datapoints
    character(len=:), allocatable :: validpath

    !> contains paths to fnetdata.xml files for training
    character(len=1024), allocatable :: datapaths(:)

    !> contains paths to fnetdata.xml files for validation
    character(len=1024), allocatable :: validpaths(:)

    !> contains paths to netstat files
    character(len=1024), allocatable :: netstatNames(:)

    !> number of systems in training dataset
    integer :: nDatapoints

    !> number of systems in validation dataset
    integer :: nValidDatapoints

    !> number of dataset species
    integer :: nSpecies

    !> total number of atoms in the dataset
    integer :: nTotalAtoms

    !> total number of available external, atomic features in dataset
    integer :: nExtFeatures

    !> number of targets per network parameter
    real(dp) :: nTargetsPerParam

    !> number of target values per atom (if tAtomicTargets = .true.) or system
    integer :: nTargets

    !> number of validation target values per atom (if tAtomicTargets = .true.) or system
    integer :: nValidTargets

    !> true, if targets are atomic properties
    logical :: tAtomicTargets

    !> true, if a contiguous dataset file is provided
    logical :: tContiguous

    !> true, if a contiguous validation dataset file is provided
    logical :: tValidContiguous

    !> true, if validation monitoring is desired
    logical :: tMonitorValid

    !> target values for training
    type(TRealArray2D), allocatable :: targets(:)

    !> target values for validation
    type(TRealArray2D), allocatable :: validTargets(:)

    !> standardized target values for training
    type(TRealArray2D), allocatable :: zTargets(:)

    !> standardized validation target values
    type(TRealArray2D), allocatable :: zValidTargets(:)

    !> storage container of means and variances to calculate z-score, shape: [nTargets, 2]
    real(dp), allocatable :: zPrec(:,:)

    !> true, if z-score standardization should be applied
    logical :: tZscore

    !> contains obtained datapoint weights
    integer, allocatable :: weights(:)

    !> contains obtained dataset geometries
    type(TGeometry), allocatable :: geos(:)

    !> contains obtained validation geometries
    type(TGeometry), allocatable :: validGeos(:)

    !> contains (unique) species of all dataset geometries
    character(len=50), allocatable :: globalSpNames(:)

    !> index mapping local species --> global species index
    type(TIntArray1D), allocatable :: localSpToGlobalSp(:)

    !> index mapping local validation species --> global species index
    type(TIntArray1D), allocatable :: localValidSpToGlobalSp(:)

    !> index mapping local atom --> global species index
    type(TIntArray1D), allocatable :: localAtToGlobalSp(:)

    !> index mapping local validation atom --> global species index
    type(TIntArray1D), allocatable :: localValidAtToGlobalSp(:)

  contains

    procedure :: getIndexMappings => TData_getIndexMappings
    procedure :: getTargetMeansAndVariances => TData_getTargetMeansAndVariances
    procedure :: applyZscore => TData_applyZscore

  end type TData


  !> Data type containing variables of the Option block
  type TExternal

    !> real valued species identifier, expected shape: [nSpecies]
    real(dp), allocatable :: speciesIds(:)

    !> index of external atomic features to use as atom identifier
    integer :: atomIdIndex

    !> real valued atom identifier for training, expected shape: [nDatapoints]
    type(TRealArray1D), allocatable :: atomIds(:)

    !> real valued atom identifier for validation, expected shape: [nValidDatapoints]
    type(TRealArray1D), allocatable :: validAtomIds(:)

    !> true, if external training features are provided and selected
    logical :: tExtFeatures

    !> true, if external validation features are provided and selected
    logical :: tExtValidFeatures

    !> number of external, atomic features
    integer :: nExtFeatures

    !> indices of external, atomic features to use
    integer, allocatable :: extFeaturesInd(:)

    !> additional external, atomic features of training dataset
    type(TRealArray2D), allocatable :: extFeatures(:)

    !> additional external, atomic features of validation dataset
    type(TRealArray2D), allocatable :: extValidFeatures(:)

  end type TExternal


  !> Data type containing variables of the Option block
  type TOption

    !> wether to resume from existing netstat files on disk
    logical :: tReadNetStats

    !> wether to write loss and gradients for all training iterations to disk
    logical :: tWriteIterTraj

    !> mode of current run (train, validate, predict)
    character(len=:), allocatable :: mode

    !> (user defined) random seed of the run
    integer :: seed

  end type TOption


  !> Data type containing the program variables
  type TProgramVariables

    !> contains mpi communicator, if compiled with mpi enabled
    type(TEnv) :: env

    !> data of Network block
    type(TArch) :: arch

    !> data of Mapping block
    type(TMapping) :: mapping

    !> acsf value instance
    type(TAcsf) :: acsf

    !> acsf value instance (for optional validation)
    type(TAcsf) :: validAcsf

    !> data of Training block
    type(TTraining) :: train

    !> data of Data block
    type(TData) :: data

    !> contains external feature and identifier information
    type(TExternal) :: ext

    !> collected features of data and mapping block
    type(TFeatures) :: features

    !> data of Option block
    type(TOption) :: option

    !> luxury pseudorandom generator instance
    type(TRanlux) :: rndGen

  contains

    procedure :: checkInputConsistency => TProgramVariables_checkInputConsistency

  end type TProgramVariables


  !> program version
  character(len=*), parameter :: version =  '0.1'

  !> copyright year
  integer, parameter :: copyrightYear = 2021

  !> root node name of the input tree
  character(len=*), parameter :: rootTag = 'fortnet'

  !> input file name
  character(len=*), parameter :: hsdInput = 'fortnet_in.hsd'

  !> parsed output name
  character(len=*), parameter :: hsdParsedInput = 'fortnet_pin.hsd'

  !> file name of generic Fortnet datapoint
  character(len=*), parameter :: fnetdataFile = 'fnetdata.xml'

  !> file name of generic Fortnet datapoint
  character(len=*), parameter :: fnetvdataFile = 'fnetvdata.xml'

  !> file name for preconditioning informations
  character(len=*), parameter :: precFile = 'precond.out'

  !> version of the input document
  integer, parameter :: parserVersion = 1


contains

  !> Initializes all necessary program variables depending on the user input.
  subroutine TProgramVariables_init(this)

    !> container of program variables
    type(TProgramVariables), intent(inout) :: this

    !> input tree node pointers
    type(fnode), pointer :: root, hsdTree, tmp

    !> string buffer instance
    type(string) :: strBuffer

    !> input version number
    integer :: inputVersion

    ! write standard output header
    call printFortnetHeader(version, copyrightYear)
    call printDateAndTime()

    ! read user input file as HSD
    call parseHSD(rootTag, hsdInput, hsdTree)
    call getChild(hsdTree, rootTag, root)

    write(stdout, '(A)') "Interpreting input file '" // hsdInput // "'"

    ! check if input version is the one, which Fortnet can handle
    call getChildValue(root, "InputVersion", inputVersion, parserVersion)
    if (inputVersion /= parserVersion) then
      call error("Version of input (" // i2c(inputVersion) // ") and parser (" &
          &// i2c(parserVersion) // ") do not match")
    end if

    ! read options
    call getChild(root, 'Options', tmp)
    call readOptions(this%option, tmp)

    ! define state of luxury random number generator
    call init(this%rndGen, luxlev=3, initSeed=this%option%seed)

    ! read data informations
    call getChild(root, 'Data', tmp)
    call readPrec(this%data, tmp, this%option%tReadNetStats)
    call readData(this%data, this%option, tmp)

    ! if present, read external features and identifiers
    call readExternal(this%ext, this%data, root)

    select case(this%option%mode)

    case ('train')

      if (this%option%tReadNetStats) then

        ! in the case of resumed training, read network informations from netstat files
        call readFromNetstats(this%arch, this%data)

      else

        ! read network informations
        call getChildValue(root, 'Network', tmp)
        call getNodeName(tmp, strBuffer)
        call readNetwork(this%arch, tmp, trim(char(strBuffer)))

        ! read mapping informations
        call getChildValue(root, 'Mapping', tmp)
        call getNodeName(tmp, strBuffer)
        call readMapping(this%mapping, this%arch, tmp, trim(char(strBuffer)), this%data%nTargets,&
            & this%ext%nExtFeatures)

      end if

      ! read training informations
      call getChildValue(root, 'Training', tmp)
      call getNodeName(tmp, strBuffer)
      call readTraining(this%train, tmp, trim(char(strBuffer)))

    case ('validate', 'predict')

      ! read informations from acsf and netstat files
      call readFromNetstats(this%arch, this%data)

    end select

    ! calculate the number of targets per sub-network parameter
    call calcTargetsPerParam(this%data, this%arch%nSubNnParams)

    ! issue warning about unprocessed nodes
    call warnUnprocessedNodes(root, .true.)

    ! check the rudimentary plausibility of the imported parameters
    call this%checkInputConsistency()

    ! finish parsing, dump parsed and processed input
    call dumpHSD(hsdTree, hsdParsedInput)

    write(stdout, '(A,/)') "Processed input written as HSD to '" // hsdParsedInput //"'"
    write(stdout, '(A,/)') repeat('-', 80)

    call destroyNode(tmp)
    call destroyNode(hsdTree)

  end subroutine TProgramVariables_init


  !> Initializes basic ACSF parameters from a given file.
  subroutine readAcsfFromFile(mapping, acsf, filename, globalSpNames, nSpecies)

    !> representation of mapping informations
    type(TMapping), intent(inout) :: mapping

    !> representation of ACSF mappings
    type(TAcsf), intent(inout) :: acsf

    !> filename or path to save acsf parameters to
    character(len=*), intent(in) :: filename

    !> contains (unique) species of all dataset geometries
    character(len=50), intent(in) :: globalSpNames(:)

    !> number of dataset species
    integer, intent(in) :: nSpecies

    mapping%type = 'acsf'

    write(stdOut, '(A)', advance='no') 'reading ACSF from file...'
    call acsf%fromFile(filename, globalSpNames)
    @:ASSERT(size(acsf%speciesIds) == nSpecies)
    write(stdout, '(A,/)') 'done'    

  end subroutine readAcsfFromFile


  !> Initializes neural sub-networks from netstat files.
  subroutine readFromNetstats(arch, data)

    !> representation of network architecture information
    type(TArch), intent(inout) :: arch

    !> representation of dataset information
    type(TData), intent(inout) :: data

    !> unique fileunit
    integer :: fd

    !> auxiliary variables
    integer :: nLayers, iSpecies
    character(len=50) :: archType, targetType
    character(len=50) :: spName
    character(len=1024) :: activation

    do iSpecies = 1, data%nSpecies

      open(newunit=fd, file=data%netstatNames(iSpecies), form='formatted', status='old',&
          & action='read')

      read(fd, *) archType, targetType
      arch%type = tolower(trim(archType))
      if (tolower(trim(targetType)) == 'atomic') then
        data%tAtomicTargets = .true.
      elseif (tolower(trim(targetType)) == 'global') then
        data%tAtomicTargets = .false.
      else
        call error("Unrecognized target type in file '" // data%netstatNames(iSpecies) // "'.")
      end if

      read(fd, *) spName
      @:ASSERT(tolower(data%globalSpNames(iSpecies)) == tolower(spName))

      read(fd, *) nLayers
      @:ASSERT(nLayers > 2)
      if (allocated(arch%allDims)) deallocate(arch%allDims)
      allocate(arch%allDims(nLayers))
      read(fd, *) arch%allDims

      arch%hidden = arch%allDims(2:size(arch%allDims) - 1)
      arch%nHiddenLayer = size(arch%hidden)

      read(fd, *) activation
      arch%activation = tolower(trim(activation))

      close(fd)

    end do

    arch%nSubNnParams = getNumberOfParameters(arch%allDims)

  end subroutine readFromNetstats


  !> Interprets the Network HSD block.
  subroutine readNetwork(arch, node, case)

    !> representation of network architecture information
    type(TArch), intent(inout) :: arch

    !> node containig the information
    type(fnode), pointer :: node

    !> type of neural network
    character(len=*), intent(in) :: case

    !> list of integers to parse hidden layer configuration
    type(TListInt) :: integerList

    !> string buffer instance to parse type of transfer function
    type(string) :: strBuffer

    select case (tolower(case))

    case ('bpnn')

      arch%type = 'bpnn'

      call init(integerList)
      call getChildValue(node, 'Hidden', integerList)
      allocate(arch%hidden(len(integerList)))
      call asArray(integerList, arch%hidden)
      call destruct(integerList)
      arch%nHiddenLayer = size(arch%hidden)

      call getChildValue(node, 'Activation', strBuffer)
      arch%activation = tolower(trim(unquote(char(strBuffer))))

    case default

      call detailedError(node, 'Invalid network type')

    end select

  end subroutine readNetwork


  !> Interprets the Mapping HSD block.
  subroutine readMapping(mapping, arch, node, case, nTargets, nExtFeatures)

    !> representation of mapping informations
    type(TMapping), intent(inout) :: mapping

    !> representation of network architecture information
    type(TArch), intent(inout) :: arch

    !> node containing the information
    type(fnode), intent(in), pointer :: node

    !> type of neural network
    character(len=*), intent(in) :: case

    !> number of target values per atom (if tAtomicTargets = .true.) or system
    integer, intent(in) :: nTargets

    !> number of external, atomic features
    integer, intent(in) :: nExtFeatures

    call getChildValue(node, 'Standardization', mapping%tStandardize, .true.)

    select case (tolower(case))

    case ('acsf')

      call getChildValue(node, 'NRadial', mapping%nRadial)
      call getChildValue(node, 'NAngular', mapping%nAngular)
      call getChildValue(node, 'RCut', mapping%rCut)

      ! convert Angstrom to Bohr
      mapping%rCut = mapping%rCut * AA__Bohr

    case default

      call detailedError(node, 'Invalid mapping type')

    end select

    arch%allDims = [mapping%nRadial + mapping%nAngular + nExtFeatures, arch%hidden, nTargets]

    arch%nSubNnParams = getNumberOfParameters(arch%allDims)

  end subroutine readMapping


  !> Interprets the External HSD block.
  subroutine readExternal(ext, data, node)

    !> representation of external feature information
    type(TExternal), intent(inout) :: ext

    !> representation of dataset information
    type(TData), intent(inout) :: data

    !> node containig the information
    type(fnode), intent(in), pointer :: node

    !> node containing informations
    type(fnode), pointer :: xml, rootNode, extnode, featuresnode, child, child1, child2

    !> string buffer instances
    type(string) :: buffer

    !> filename of .xml files
    character(len=:), allocatable :: filename

    !> auxiliary variables
    integer :: iSp, iSys, nTmpExtFeatures

    ! inquire external training dataset features
    if (data%tContiguous) then
      call readHSDAsXML(data%datapaths(1), xml)
      call getChild(xml, 'fnetdata', rootNode)
      call inquireFeatures(rootNode, data%nExtFeatures)
      call destroyNode(xml)
    else
      data%nExtFeatures = 0
      nTmpExtFeatures = data%nExtFeatures
      do iSys = 1, data%nDatapoints
        filename = trim(data%datapaths(iSys)) // '/' // fnetdataFile
        call readHSDAsXML(filename, xml)
        call getChild(xml, 'fnetdata', rootNode)
        call inquireFeatures(rootNode, data%nExtFeatures)
        if ((nTmpExtFeatures /= data%nExtFeatures) .and. (iSys /= 1)) then
          call error('Inconsistency in number of external features of dataset found.')
        end if
        nTmpExtFeatures = data%nExtFeatures
        call destroyNode(xml)
      end do
    end if

    ! inquire external validation dataset features
    if (data%tMonitorValid) then
      if (data%tValidContiguous) then
        nTmpExtFeatures = data%nExtFeatures
        call readHSDAsXML(data%validpaths(1), xml)
        call getChild(xml, 'fnetdata', rootNode)
        call inquireFeatures(rootNode, data%nExtFeatures)
        call destroyNode(xml)
        if (nTmpExtFeatures /= data%nExtFeatures) then
          call error('Inconsistency in number of external features between datasets found.')
        end if
        else
          nTmpExtFeatures = data%nExtFeatures
          do iSys = 1, data%nValidDatapoints
            filename = trim(data%validpaths(iSys)) // '/' // fnetdataFile
            call readHSDAsXML(filename, xml)
            call getChild(xml, 'fnetdata', rootNode)
            call inquireFeatures(rootNode, data%nExtFeatures)
            if (nTmpExtFeatures /= data%nExtFeatures) then
              call error('Inconsistency in number of external features between datasets found.')
            end if
            nTmpExtFeatures = data%nExtFeatures
            call destroyNode(xml)
          end do
        end if
      end if

    call getChild(node, 'External', extnode, requested=.false.)

    if (associated(extnode)) then

      ! read species identifier
      allocate(ext%speciesIds(data%nSpecies))
      call getChild(extnode, 'SpeciesID', child1, requested=.false.)
      if (associated(child1)) then
        do iSp = 1, data%nSpecies
          call getChildValue(child1, data%globalSpNames(iSp), ext%speciesIds(iSp), 1.0_dp)
          if (ext%speciesIds(iSp) <= 0.0_dp) then
            call warning("Obtained potentially dangerous '" // trim(data%globalSpNames(iSp))&
                & // "' species" // NEW_LINE('A') //&
                & '   identifier equal or below zero, watch out.')
          end if
        end do
      else
        ext%speciesIds(:) = 1.0_dp
      end if

      call getChild(extnode, 'Features', featuresnode, requested=.false.)
      if (associated(featuresnode)) then
        call getChildValue(featuresnode, '', buffer, child=child, multiple=.true.)
        call convRangeToInt(char(buffer), featuresnode, ext%extFeaturesInd, data%nExtFeatures)
        call setChildValue(child, '', ext%extFeaturesInd, replace=.true.)
        ext%nExtFeatures = size(ext%extFeaturesInd)
        if (size(ext%extFeaturesInd) > 0) then
          ext%tExtFeatures = .true.
        else
          ext%tExtFeatures = .false.
          ext%nExtFeatures = 0
        end if
      else
        ext%tExtFeatures = .false.
        ext%nExtFeatures = 0
      end if

      ! read external features from generic fnetdata.xml file(s)
      if (data%tContiguous .and. ext%tExtFeatures) then
        call readHSDAsXML(data%datapaths(1), xml)
        call getChild(xml, 'fnetdata', rootNode)
        call readContiguousFnetdataFeatures(rootNode, data%geos, ext%extFeatures,&
            & inds=ext%extFeaturesInd)
        call destroyNode(xml)
      elseif ((.not. data%tContiguous) .and. ext%tExtFeatures) then
        allocate(ext%extFeatures(data%nDatapoints))
        do iSys = 1, data%nDatapoints
          filename = trim(data%datapaths(iSys)) // '/' // fnetdataFile
          call readHSDAsXML(filename, xml)
          call getChild(xml, 'fnetdata', rootNode)
          call readFnetdataFeatures(rootNode, data%geos(iSys)%nAtom, ext%extFeatures(iSys)%array,&
              & inds=ext%extFeaturesInd)
          call destroyNode(xml)
        end do
      end if

      ! if present, read external validation features from generic fnetdata.xml file(s)
      if (data%tMonitorValid .and. ext%tExtFeatures .and. data%tValidContiguous) then
        call readHSDAsXML(data%validpaths(1), xml)
        call getChild(xml, 'fnetdata', rootNode)
        call readContiguousFnetdataFeatures(rootNode, data%validGeos, ext%extValidFeatures,&
            & inds=ext%extFeaturesInd)
        call destroyNode(xml)
      elseif (data%tMonitorValid .and. ext%tExtFeatures .and. (.not. data%tValidContiguous)) then
        allocate(ext%extValidFeatures(data%nValidDatapoints))
        do iSys = 1, data%nValidDatapoints
          filename = trim(data%validpaths(iSys)) // '/fnetdata.xml'
          call readHSDAsXML(filename, xml)
          call getChild(xml, 'fnetdata', rootNode)
          call readFnetdataFeatures(rootNode, data%validGeos(iSys)%nAtom,&
              & ext%extValidFeatures(iSys)%array, inds=ext%extFeaturesInd)
          call destroyNode(xml)
        end do
      end if

      ! read atom identifier
      call getChild(extnode, 'AtomID', child2, requested=.false.)
      if (associated(child2)) then
        call getChildValue(child2, '', ext%atomIdIndex)
        if (ext%atomIdIndex <= 0) then
          call warning('Atom identifier index less or equal zero specified, will be ignored.')
        elseif (ext%atomIdIndex > data%nExtFeatures) then
          call error('Atom identifier index exceeds number of external features in dataset.')
        end if
      else
        ext%atomIdIndex = 0
      end if

    else

      allocate(ext%speciesIds(data%nSpecies))
      ext%speciesIds(:) = 1.0_dp
      ext%atomIdIndex = 0
      ext%tExtFeatures = .false.
      ext%nExtFeatures = 0

    end if

    if (ext%atomIdIndex > 0) then
      if (data%tContiguous) then
        call readHSDAsXML(data%datapaths(1), xml)
        call getChild(xml, 'fnetdata', rootNode)
        call readContiguousFnetdataAtomIdentifier(rootNode, data%geos, ext%atomIdIndex, ext%atomIds)
        call destroyNode(xml)
      else
        allocate(ext%atomIds(data%nDatapoints))
        do iSys = 1, data%nDatapoints
          filename = trim(data%datapaths(iSys)) // '/' // fnetdataFile
          call readHSDAsXML(filename, xml)
          call getChild(xml, 'fnetdata', rootNode)
          call readFnetdataAtomIdentifier(rootNode, data%geos(iSys)%nAtom, ext%atomIdIndex,&
              & ext%atomIds(iSys)%array)
          call destroyNode(xml)
        end do
      end if
    else
      allocate(ext%atomIds(data%nDatapoints))
      do iSys = 1, data%nDatapoints
        allocate(ext%atomIds(iSys)%array(data%geos(iSys)%nAtom))
        ext%atomIds(iSys)%array(:) = 1.0_dp
      end do
    end if

    if (data%tMonitorValid) then
      if (ext%atomIdIndex > 0) then
        if (data%tValidContiguous) then
          call readHSDAsXML(data%validpaths(1), xml)
          call getChild(xml, 'fnetdata', rootNode)
          call readContiguousFnetdataAtomIdentifier(rootNode, data%validGeos, ext%atomIdIndex,&
              & ext%validAtomIds)
          call destroyNode(xml)
        else
          allocate(ext%validAtomIds(data%nValidDatapoints))
          do iSys = 1, data%nValidDatapoints
            filename = trim(data%validpaths(iSys)) // '/' // fnetdataFile
            call readHSDAsXML(filename, xml)
            call getChild(xml, 'fnetdata', rootNode)
            call readFnetdataAtomIdentifier(rootNode, data%validGeos(iSys)%nAtom, ext%atomIdIndex,&
                & ext%validAtomIds(iSys)%array)
            call destroyNode(xml)
          end do
        end if
      else
        allocate(ext%validAtomIds(data%nValidDatapoints))
        do iSys = 1, data%nValidDatapoints
          allocate(ext%validAtomIds(iSys)%array(data%validGeos(iSys)%nAtom))
          ext%validAtomIds(iSys)%array(:) = 1.0_dp
        end do
      end if
    end if

  end subroutine readExternal


  !> Interprets the Training HSD block.
  subroutine readTraining(train, node, case)

    !> representation of training/optimizer information
    type(TTraining), intent(inout) :: train

    !> node containig the information
    type(fnode), intent(in), pointer :: node

    !> type of neural network
    character(len=*), intent(in) :: case

    !> string buffer instance
    type(string) :: strBuffer

    select case (tolower(case))

    case ('sd')

      train%iOptimizer = optimizerTypes%steepDesc
      call getChildValue(node, 'LearningRate', train%learningRate, 0.01_dp)

    case ('cg')

      train%iOptimizer = optimizerTypes%conjGrad

    case ('lbfgs')

      train%iOptimizer = optimizerTypes%lbfgs
      call getChildValue(node, 'MaxForQNDisplacement', train%maxForQNDisplacement, .false.)
      call getChildValue(node, 'Linemin', train%tLinesearch, .true.)
      call getChildValue(node, 'Memory', train%mem, 1000)

    case ('fire')

      train%iOptimizer = optimizerTypes%fire

    case default

      call detailedError(node, 'Invalid training algorithm')

    end select

    call getChildValue(node, 'MinDisplacement', train%minDisplacement, 1e-06_dp)
    call getChildValue(node, 'MaxDisplacement', train%maxDisplacement, 1e+04_dp)

    call getChildValue(node, 'NIterations', train%nTrainIt, huge(0))
    call getChildValue(node, 'Threshold', train%threshold, tiny(0.0_dp))

    call getChildValue(node, 'NPrintout', train%nPrintOut, 10)
    call getChildValue(node, 'NSaveNet', train%nSaveNet, 100)

    call getChildValue(node, 'Loss', strBuffer, 'rms')
    call train%setLossFunc(tolower(trim(unquote(char(strBuffer)))))

  end subroutine readTraining


  !> Interprets the Options HSD block.
  subroutine readOptions(option, node)

    !> representation of user specified options
    type(TOption), intent(inout) :: option

    !> node containig the information
    type(fnode), pointer :: node

    !> string buffer instance
    type(string) :: strBuffer

    !> temporary random integer
    integer :: tmpIntSeed

    !> auxiliary variable
    real(dp) :: tmpRealSeed

    call getChildValue(node, 'ReadNetStats', option%tReadNetStats, .false.)

    call getChildValue(node, 'WriteIterationTrajectory', option%tWriteIterTraj, .false.)

    call getChildValue(node, 'Mode', strBuffer)
    option%mode = tolower(unquote(char(strBuffer)))

    if ((trim(option%mode) == 'validate' .or. trim(option%mode) == 'predict') .and.&
        & (option%tWriteIterTraj .eqv. .true.)) then
      write(stdout, '(A)') ''
      call warning('Running in validation or prediction mode does not produce an iteration'&
          & //NEW_LINE('A')//'   trajectory. Overwriting user input...')
      write(stdout, '(A)') ''
      option%tWriteIterTraj = .false.
      call setChildValue(node, 'WriteIterationTrajectory', option%tWriteIterTraj, replace=.true.)
    end if

    if ((trim(option%mode) == 'validate' .or. trim(option%mode) == 'predict') .and.&
        & (option%tReadNetStats .eqv. .false.)) then
      write(stdout, '(A)') ''
      call warning('Running in validation or prediction mode without initialising from'&
          & //NEW_LINE('A')//'   existing netstat files is not possible. Overwriting user input...')
      write(stdout, '(A)') ''
      option%tReadNetStats = .true.
      call setChildValue(node, 'ReadNetStats', option%tReadNetStats, replace=.true.)
    end if

    call random_number(tmpRealSeed)
    tmpIntSeed = floor(tmpRealSeed * huge(0) + tiny(0.0))

    call getChildValue(node, 'RandomSeed', option%seed, tmpIntSeed)

    if (option%seed < 0) then
      call detailedError(node, 'Random seed must be greater or equal zero.')
    end if

  end subroutine readOptions


  !> Interprets and reads preconditioning of the Data HSD block.
  subroutine readPrec(data, node, tReadNetStats)

    !> representation of dataset information
    type(TData), intent(inout) :: data

    !> node containig the information
    type(fnode), pointer :: node, child

    !> wether to resume from existing netstat files on disk
    logical, intent(in) :: tReadNetStats

    !> true, if preconditioning file is in place
    logical :: tExist

    call getChildValue(node, 'Standardization', data%tZscore, .false., child=child)

    inquire(file=precFile, exist=tExist)

    if (.not. data%tZscore .and. tExist) then
      call warning('User input manually deactivated target standardization,'&
          & //NEW_LINE('A')// '   but preconditioning parameters are present.')
    end if

    call readPreconditioning(precFile, tReadNetStats, data%tZscore, data%nTargets, data%zPrec)

  end subroutine readPrec


  !> Reads datapoint paths from file.
  subroutine readDatapathsFromFile(fname, node, nDatapoints, datapaths)

    !> name of file that contains the datapoint paths
    character(len=*), intent(in) :: fname

    !> node containig the information
    type(fnode), pointer :: node

    !> total number of datapoints found
    integer, intent(out) :: nDatapoints

    !> paths to .xml files with geometry and target informations
    character(len=1024), intent(out), allocatable :: datapaths(:)

    !> file and error identifier
    integer :: fp, iErr

    !> auxiliary variable
    integer :: iSys

    open(newunit=fp, file=fname, form='formatted', status='old', action='read', iostat=iErr)

    if (iErr /= 0) then
      call detailedError(node, "Could not open file '" // fname // "' for direct reading" )
    end if

    read(fp, *, iostat=iErr) nDatapoints

    allocate(datapaths(nDatapoints))
    do iSys = 1, nDatapoints
      read(fp, '(A)', iostat=iErr) datapaths(iSys)
    end do

    if (iErr /= 0) then
      call detailedError(node, "Error during direct reading '" // fname // "'")
    end if

    close(fp)

  end subroutine readDatapathsFromFile


  !> Interprets the Data HSD block.
  subroutine readData(data, option, node)

    !> representation of dataset information
    type(TData), intent(inout) :: data

    !> representation of user specified options
    type(TOption), intent(in) :: option

    !> node containig the information
    type(fnode), pointer :: node

    !> string buffer instances
    type(string) :: strBuffer, buffer1, buffer2

    !> nodes containig the information
    type(fnode), pointer :: child1, child2, value1, xml, rootNode

    !> list of strings
    type(TListString) :: lStr

    !> characters
    character(lc) :: prefix, suffix, elem, strTmp

    !> filename of .xml files
    character(len=:), allocatable :: filename

    !> auxiliary variables
    integer :: iSys, iSp, ii

    !> if the name of the types should be converted to lower case
    !> (otherwise they are used in the same way, specified in the geometry input)
    logical :: tLower

    !> true, if current netstat file is in place
    logical :: tExist

    call getChildValue(node, 'Dataset', strBuffer)
    data%datapath = trim(unquote(char(strBuffer)))

    if (data%datapath(len(data%datapath)-len(fnetdataFile)+1:len(data%datapath))&
        & == fnetdataFile) then
      data%tContiguous = .true.
      allocate(data%datapaths(1))
      data%datapaths(1) = data%datapath
    else
      data%tContiguous = .false.
      call readDatapathsFromFile(data%datapath, node, data%nDatapoints, data%datapaths)
    end if

    call getChildValue(node, 'Validset', strBuffer, default='')

    select case (char(strBuffer))
    case ('')
      data%tMonitorValid = .false.
    case ('none')
      data%tMonitorValid = .false.
    case default
      data%tMonitorValid = .true.
      data%validpath = trim(unquote(char(strBuffer)))
      if (data%validpath(len(data%validpath)-len(fnetvdataFile)+1:len(data%validpath))&
          & == fnetvdataFile) then
        data%tValidContiguous = .true.
        allocate(data%validpaths(1))
        data%validpaths(1) = data%validpath
      else
        data%tValidContiguous = .false.
        call readDatapathsFromFile(data%validpath, node, data%nValidDatapoints, data%validpaths)
      end if
    end select

    ! read geometries from generic fnetdata.xml file(s)
    if (data%tContiguous) then
      call readHSDAsXML(data%datapaths(1), xml)
      call getChild(xml, 'fnetdata', rootNode)
      call readContiguousFnetdataWeights(rootNode, data%weights)
      call readContiguousFnetdataGeometries(rootNode, data%geos)
      data%nDatapoints = size(data%geos)
      select case (option%mode)
      case('train', 'validate')
        call readContiguousFnetdataTargets(rootNode, data%geos, data%targets, data%nTargets,&
            & data%tAtomicTargets)
      end select
      call destroyNode(xml)
    else
      allocate(data%weights(data%nDatapoints))
      allocate(data%geos(data%nDatapoints))
      select case (option%mode)
      case('train', 'validate')
        allocate(data%targets(data%nDatapoints))
      end select
      do iSys = 1, data%nDatapoints
        filename = trim(data%datapaths(iSys)) // '/' // fnetdataFile
        call readHSDAsXML(filename, xml)
        call getChild(xml, 'fnetdata', rootNode)
        call readFnetdataWeight(rootNode, data%weights(iSys))
        call readFnetdataGeometry(rootNode, data%geos(iSys))
        select case (option%mode)
        case('train', 'validate')
          call readFnetdataTargets(rootNode, data%geos(iSys)%nAtom, data%targets(iSys)%array,&
              & data%tAtomicTargets)
          data%nTargets = size(data%targets(iSys)%array, dim=1)
        end select
        call destroyNode(xml)
      end do
    end if

    data%nTotalAtoms = getTotalNumberOfAtoms(data%geos)

    ! if present, read validation geometries from generic fnetdata.xml file(s)
    if (data%tMonitorValid) then
      if (data%tValidContiguous) then
        call readHSDAsXML(data%validpaths(1), xml)
        call getChild(xml, 'fnetdata', rootNode)
        call readContiguousFnetdataGeometries(rootNode, data%validGeos)
        data%nValidDatapoints = size(data%validGeos)
        select case (option%mode)
        case('train', 'validate')
          call readContiguousFnetdataTargets(rootNode, data%validGeos, data%validTargets,&
              & data%nValidTargets, data%tAtomicTargets)
        end select
        call destroyNode(xml)
      else
        allocate(data%validGeos(data%nValidDatapoints))
        select case (option%mode)
        case('train', 'validate')
          allocate(data%validTargets(data%nValidDatapoints))
        end select
        do iSys = 1, data%nValidDatapoints
          filename = trim(data%validpaths(iSys)) // '/fnetdata.xml'
          call readHSDAsXML(filename, xml)
          call getChild(xml, 'fnetdata', rootNode)
          call readFnetdataGeometry(rootNode, data%validGeos(iSys))
          select case (option%mode)
          case('train', 'validate')
            call readFnetdataTargets(rootNode, data%validGeos(iSys)%nAtom,&
                & data%validTargets(iSys)%array, data%tAtomicTargets)
            data%nValidTargets = size(data%validTargets(iSys)%array, dim=1)
          end select
          call destroyNode(xml)
        end do
      end if
    end if

    ! apply z-score standardization, if desired
    if (data%tZscore) then
      if (.not. allocated(data%zPrec)) then
        call data%getTargetMeansAndVariances()
        call writePreconditioning(precFile, data%zPrec)
      end if
      call data%applyZscore()
    elseif (.not. option%mode == 'predict') then
      data%zTargets = data%targets
      if (data%tMonitorValid) then
        data%zValidTargets = data%validTargets
      end if
    end if

    ! establish convenient index mappings and global species list
    call data%getIndexMappings()
    allocate(data%netstatNames(data%nSpecies))

    ! netstat file names
    call getChildValue(node, 'NetstatFiles', value1, child=child1)
    call getNodeName(value1, buffer1)

    select case(char(buffer1))

    case ('type2filenames')

      call getChildValue(value1, 'Prefix', buffer2, '')
      prefix = unquote(char(buffer2))
      call getChildValue(value1, 'Suffix', buffer2, '')
      suffix = unquote(char(buffer2))

      call getChildValue(value1, 'LowerCaseTypeName', tLower, .false.)

      do iSp = 1, data%nSpecies
        if (tLower) then
          elem = tolower(data%globalSpNames(iSp))
        else
          elem = data%globalSpNames(iSp)
        end if

        strTmp = trim(prefix) // trim(elem) // trim(suffix)
        data%netstatNames(iSp) = strTmp

        if (option%tReadNetStats) then
          inquire(file=strTmp, exist=tExist)
        end if

        if (option%tReadNetStats .and. (.not. tExist)) then
          call detailedError(value1, "Netstat file with generated name '" // trim(strTmp) //&
              & "' does not exist.")
        end if

      end do

    case default

      call setUnprocessed(value1)

      do iSp = 1, data%nSpecies

          call init(lStr)
          call getChildValue(child1, trim(data%globalSpNames(iSp)), lStr, child=child2)

          data%netstatNames(iSp) = trim(data%globalSpNames(iSp))

          if (len(lStr) /= data%nSpecies) then
            call detailedError(child2, "Incorrect number of netstat files")
          end if

          do ii = 1, len(lStr)
            call get(lStr, data%globalSpNames(iSp), ii)

            if (option%tReadNetStats) then
              inquire(file=strTmp, exist=tExist)
            end if

            if (option%tReadNetStats .and. (.not. tExist)) then
              call detailedError(child2, "Netstat file '" // trim(data%globalSpNames(iSp))&
                  & // "' does not exist'")
            end if
          end do

          call destruct(lStr)

        end do

      end select

  end subroutine readData


  !> Calculates means and variances for z-score standardization.
  subroutine TData_getTargetMeansAndVariances(this)

    !> target values for training
    class(TData), intent(inout) :: this

    !> auxiliary variables
    integer :: iSys, iAtom, iTarget, nTotAtoms

    if (allocated(this%zPrec)) then
      call error('Container for target means and variances is already allocated.')
    else
      allocate(this%zPrec(this%nTargets, 2))
      this%zPrec(:,:) = 0.0_dp    
    end if

    nTotAtoms = 0

    ! calculate means of the different outputs
    do iSys = 1, size(this%targets)
      do iAtom = 1, size(this%targets(iSys)%array, dim=2)
        nTotAtoms = nTotAtoms + this%weights(iSys)
        do iTarget = 1, this%nTargets
          this%zPrec(iTarget, 1) = this%zPrec(iTarget, 1) + real(this%weights(iSys), dp)&
              & * this%targets(iSys)%array(iTarget, iAtom)
        end do
      end do
    end do

    this%zPrec(:, 1) = this%zPrec(:, 1) / real(nTotAtoms, dp)

    ! calculate variances of the different outputs
    do iSys = 1, size(this%targets)
      do iAtom = 1, size(this%targets(iSys)%array, dim=2)
        do iTarget = 1, this%nTargets
          this%zPrec(iTarget, 2) = this%zPrec(iTarget, 2) + real(this%weights(iSys), dp)&
              & * (this%targets(iSys)%array(iTarget, iAtom) - this%zPrec(iTarget, 1))**2
        end do
      end do
    end do

    this%zPrec(:, 2) = sqrt(this%zPrec(:, 2) / real(nTotAtoms, dp))

    ! correct for vanishing variances (in those cases, effectively do nothing)
    do iTarget = 1, this%nTargets
      if (this%zPrec(iTarget, 2) < 1e-08) then
        this%zPrec(iTarget, 1) = 0.0_dp
        this%zPrec(iTarget, 2) = 1.0_dp
      end if
    end do

  end subroutine TData_getTargetMeansAndVariances


  !> applies z-score standardization to output data
  subroutine TData_applyZscore(this)

    !> target values for training
    class(TData), intent(inout) :: this

    !> auxiliary variables
    integer :: iSys, iAtom

    if (.not. allocated(this%zPrec)) then
      call error('Cannot apply acsf z-score standardization without means and variances.')   
    end if

    allocate(this%zTargets(size(this%targets)))
    if (this%tMonitorValid) then
      allocate(this%zValidTargets(size(this%validTargets)))
    end if

    ! apply z-score standardization to training outputs
    do iSys = 1, size(this%targets)
      allocate(this%zTargets(iSys)%array(this%nTargets, size(this%targets(iSys)%array, dim=2)))
      do iAtom = 1, size(this%targets(iSys)%array, dim=2)
        this%zTargets(iSys)%array(:, iAtom) = (this%targets(iSys)%array(:, iAtom) -&
            & this%zPrec(:, 1)) / this%zPrec(:, 2)
      end do
    end do

    ! if desired, apply z-score standardization to validation outputs
    if (this%tMonitorValid) then
      do iSys = 1, size(this%validTargets)
        allocate(this%zValidTargets(iSys)%array(this%nTargets,&
            & size(this%validTargets(iSys)%array, dim=2)))
        do iAtom = 1, size(this%validTargets(iSys)%array, dim=2)
          this%zValidTargets(iSys)%array(:, iAtom) =&
              & (this%validTargets(iSys)%array(:, iAtom) - this%zPrec(:, 1)) / this%zPrec(:, 2)
        end do
      end do
    end if

  end subroutine TData_applyZscore


  !> interpret the Training block
  subroutine initOptimizer(train, initialGuess)

    !> representation of training/optimizer information
    type(TTraining), intent(inout) :: train

    !> starting points, shape: [nPoints, nStrucs]
    real(dp), intent(in) :: initialGuess(:,:)

    !> steepest descent optimizer
    type(TWrapSteepDesc), allocatable :: wrapSteepDesc(:)

    !> conjugate gradients optimizer
    type(TWrapConjGrad), allocatable :: wrapConjGrad(:)

    !> limited memory bfgs optimizer
    type(TWrapLbfgs), allocatable :: wrapLbfgs(:)

    !> fire optimizer
    type(TWrapFire), allocatable :: wrapFire(:)

    !> weights for steepest descent optimizer
    real(dp), allocatable :: weights(:)

    !> auxiliary variable
    integer :: iStruc, nStrucs, nValues

    nValues = size(initialGuess, dim=1)
    nStrucs = size(initialGuess, dim=2)

    allocate(train%pOptimizer(nStrucs))

    select case (train%iOptimizer)

    case(optimizerTypes%steepDesc)

      allocate(wrapSteepDesc(nStrucs))
      allocate(weights(nValues))
      weights(:) = train%learningRate

      do iStruc = 1, nStrucs
        allocate(wrapSteepDesc(iStruc)%pSteepDesc)
        call init(wrapSteepDesc(iStruc)%pSteepDesc, nValues, train%threshold,&
            & train%maxDisplacement, weights)
        call init(wrapSteepDesc(iStruc)%pSteepDesc, train%pOptimizer(iStruc))
      end do

      deallocate(weights)

    case (optimizerTypes%conjGrad)

      allocate(wrapConjGrad(nStrucs))

      do iStruc = 1, nStrucs
        allocate(wrapConjGrad(iStruc)%pConjGrad)
        call init(wrapConjGrad(iStruc)%pConjGrad, nValues, train%threshold, train%maxDisplacement)
        call init(wrapConjGrad(iStruc)%pConjGrad, train%pOptimizer(iStruc))
      end do

    case (optimizerTypes%lbfgs)

      allocate(wrapLbfgs(nStrucs))

      do iStruc = 1, nStrucs
        allocate(wrapLbfgs(iStruc)%pLbfgs)
        call TLbfgs_init(wrapLbfgs(iStruc)%pLbfgs, nValues, train%threshold, train%minDisplacement,&
            & train%maxDisplacement, train%mem, train%tLinesearch, .false.,&
            & train%maxForQNDisplacement)
        call init(wrapLbfgs(iStruc)%pLbfgs, train%pOptimizer(iStruc))
      end do

    case (optimizerTypes%fire)

      allocate(wrapFire(nStrucs))

      do iStruc = 1, nStrucs
        allocate(wrapFire(iStruc)%pFire)
        call TFire_init(wrapFire(iStruc)%pFire, nValues, train%threshold, train%maxDisplacement)
        call init(wrapFire(iStruc)%pFire, train%pOptimizer(iStruc))
      end do

    end select

    do iStruc = 1, nStrucs
      call reset(train%pOptimizer(iStruc), initialGuess(:, iStruc))
    end do

  end subroutine initOptimizer


  subroutine TFeatures_init(features, nFeatures, data, ext, acsf, validAcsf)

    !> collected features of data and mapping block
    type(TFeatures), intent(inout) :: features

    !> total number of input features
    integer, intent(out) :: nFeatures

    !> representation of dataset informations
    type(TData), intent(in) :: data

    !> representation of external informations
    type(TExternal), intent(in) :: ext

    !> representation of acsf mapping informations
    type(TAcsf), intent(in) :: acsf

    !> representation of acsf mapping informations
    type(TAcsf), intent(in) :: validAcsf

    !> auxiliary variables
    integer :: iGeo, nAcsf

    nAcsf = acsf%nRadial + acsf%nAngular
    nFeatures = nAcsf + ext%nExtFeatures
    allocate(features%features(data%nDatapoints))

    do iGeo = 1, data%nDatapoints
      allocate(features%features(iGeo)%array(nFeatures, data%geos(iGeo)%nAtom))
      features%features(iGeo)%array(:,:) = 0.0_dp
    end do

    if (data%tMonitorValid) then

      nAcsf = validAcsf%nRadial + validAcsf%nAngular
      nFeatures = nAcsf + ext%nExtFeatures
      allocate(features%validFeatures(data%nValidDatapoints))

      do iGeo = 1, size(data%validGeos)
        allocate(features%validFeatures(iGeo)%array(nFeatures, data%validGeos(iGeo)%nAtom))
        features%validFeatures(iGeo)%array(:,:) = 0.0_dp
      end do
    end if

  end subroutine TFeatures_init


  subroutine TFeatures_collect(features, data, ext, acsf, validAcsf)

    !> collected features of data and mapping block
    type(TFeatures), intent(inout) :: features

    !> representation of dataset informations
    type(TData), intent(in) :: data

    !> representation of external informations
    type(TExternal), intent(in) :: ext

    !> representation of acsf mapping informations
    type(TAcsf), intent(in) :: acsf

    !> representation of acsf mapping informations
    type(TAcsf), intent(in) :: validAcsf

    !> auxiliary variables
    integer :: iGeo, nAcsf, nFeatures

    nAcsf = acsf%nRadial + acsf%nAngular
    nFeatures = nAcsf + ext%nExtFeatures

    if (ext%tExtFeatures) then
      do iGeo = 1, data%nDatapoints
        features%features(iGeo)%array(1:nAcsf, :) = acsf%vals%vals(iGeo)%array
        features%features(iGeo)%array(nAcsf+1:nFeatures, :) = ext%extFeatures(iGeo)%array
      end do
    else
      features%features = acsf%vals%vals(:)
    end if

    if (data%tMonitorValid) then
      nAcsf = validAcsf%nRadial + validAcsf%nAngular
      nFeatures = nAcsf + ext%nExtFeatures
      if (ext%tExtFeatures) then
        do iGeo = 1, size(data%validGeos)
          features%validFeatures(iGeo)%array(1:nAcsf, :) = validAcsf%vals%vals(iGeo)%array
          features%validFeatures(iGeo)%array(nAcsf+1:nFeatures, :) =&
              & ext%extValidFeatures(iGeo)%array
        end do
      else
        features%validFeatures = validAcsf%vals%vals(:)
      end if
    end if

  end subroutine TFeatures_collect


#:if WITH_MPI
  !> Synchronizes the collected features between the MPI nodes.
  subroutine TFeatures_sync(this, comm)

    !> collected features of data and mapping block
    class(TFeatures), intent(in) :: this

    !> mpi communicator with some additional information
    type(mpifx_comm), intent(in) :: comm

    !> auxiliary variable
    integer :: iGeo

    do iGeo = 1, size(this%features)
      call mpifx_bcast(comm, this%features(iGeo)%array)
    end do

    if (allocated(this%validFeatures)) then
      do iGeo = 1, size(this%validFeatures)
        call mpifx_bcast(comm, this%validFeatures(iGeo)%array)
      end do
    end if

  end subroutine TFeatures_sync
#:endif


  !> Calculates the total number of parameters per sub-network.
  pure function getNumberOfParameters(allDims) result(nSubNnParams)

    !> number of nodes per layer, including in- and output, expected shape: [nHiddenLayer + 2]
    integer, intent(in) :: allDims(:)

    !> number of paramaters (weights + biases) per sub-nn
    integer :: nSubNnParams

    !> auxiliary variable
    integer :: iLayer

    nSubNnParams = sum(allDims(2:))

    do iLayer = 1, size(allDims) - 1
      nSubNnParams = nSubNnParams + allDims(iLayer) * allDims(iLayer + 1)
    end do

  end function getNumberOfParameters


  !> Calculates the total number of atoms in a list of geometry instances.
  pure function getTotalNumberOfAtoms(geos) result(nTotalAtoms)

    !> contains the geometry information
    type(TGeometry), intent(in) :: geos(:)

    !> total number of atoms in the geometry list
    integer :: nTotalAtoms

    !> auxiliary variable
    integer :: iGeo

    nTotalAtoms = 0

    do iGeo = 1, size(geos)
      nTotalAtoms = nTotalAtoms + geos(iGeo)%nAtom
    end do

  end function getTotalNumberOfAtoms


  !> Calculates the number of targets per sub-network parameter.
  subroutine calcTargetsPerParam(data, nSubNnParams)

    !> representation of dataset information
    type(TData), intent(inout) :: data

    !> number of network paramaters (weights + biases) per sub-nn
    integer, intent(in) :: nSubNnParams

    ! calculate number of targets per network parameter
    if (data%tAtomicTargets) then
      data%nTargetsPerParam = real(data%nTargets * data%nTotalAtoms, dp) /&
          & real(data%nSpecies * nSubNnParams , dp)
    else
      data%nTargetsPerParam = real(data%nTargets * data%nDatapoints, dp) /&
          & real(data%nSpecies * nSubNnParams, dp)
    end if

  end subroutine calcTargetsPerParam


  !> Set the user defined loss function type.
  subroutine TTraining_setLossFunc(this, descriptor)

    !> data type containing variables of the Training block
    class(TTraining), intent(inout) :: this

    !> descriptor of loss function to use
    character(len=*), intent(in) :: descriptor

    select case(trim(descriptor))

      case('rms')
        this%loss => rmsLoss
        this%lossType = 'rms'

      case('mae')
        this%loss => maLoss
        this%lossType = 'mae'

      case('mse')
        this%loss => msLoss
        this%lossType = 'mse'

      case('msle')
        this%loss => mslLoss
        this%lossType = 'msle'

      case default
        this%loss => rmsLoss
        this%lossType = 'rms'

    end select

  end subroutine TTraining_setLossFunc


  !> Establishes convenient index mappings.
  subroutine TData_getIndexMappings(this)

    !> representation of dataset information
    class(TData), intent(inout) :: this

    !> auxiliary variables
    integer :: iGeo, iSp, iAtom, localSp, ind

    allocate(this%globalSpNames(1))
    this%globalSpNames(1) = this%geos(1)%speciesNames(1)

    allocate(this%localSpToGlobalSp(this%nDatapoints))
    if (this%tMonitorValid) then

    end if

    ind = 1

    do iGeo = 1, this%nDatapoints
      allocate(this%localSpToGlobalSp(iGeo)%array(this%geos(iGeo)%nSpecies))
      do iSp = 1, this%geos(iGeo)%nSpecies
        if (.not. any(this%globalSpNames .eq. this%geos(iGeo)%speciesNames(iSp))) then
          this%globalSpNames = [this%globalSpNames, this%geos(iGeo)%speciesNames(iSp)]
          ind = ind + 1
        end if
        ! work around the findloc intrinsic
        ! this%localSpToGlobalSp(iGeo)%array(iSp) = findloc(this%globalSpNames,&
        !     & this%geos(iGeo)%speciesNames(iSp), dim=1)
        this%localSpToGlobalSp(iGeo)%array(iSp) = myFindloc(this%globalSpNames,&
            & this%geos(iGeo)%speciesNames(iSp))
      end do
    end do

    this%nSpecies = ind

    if (this%tMonitorValid) then
      allocate(this%localValidSpToGlobalSp(this%nValidDatapoints))
      do iGeo = 1, this%nValidDatapoints
        allocate(this%localValidSpToGlobalSp(iGeo)%array(this%validGeos(iGeo)%nSpecies))
        do iSp = 1, this%validGeos(iGeo)%nSpecies
          ! work around the findloc intrinsic
          ! this%localValidSpToGlobalSp(iGeo)%array(iSp) = findloc(this%globalSpNames,&
          !     & this%validGeos(iGeo)%speciesNames(iSp), dim=1)
          this%localValidSpToGlobalSp(iGeo)%array(iSp) = myFindloc(this%globalSpNames,&
              & this%validGeos(iGeo)%speciesNames(iSp))
        end do
      end do
    end if

    allocate(this%localAtToGlobalSp(this%nDatapoints))

    do iGeo = 1, this%nDatapoints
      allocate(this%localAtToGlobalSp(iGeo)%array(this%geos(iGeo)%nAtom))
      do iAtom = 1, this%geos(iGeo)%nAtom
        localSp = this%geos(iGeo)%species(iAtom)
        this%localAtToGlobalSp(iGeo)%array(iAtom) = this%localSpToGlobalSp(iGeo)%array(localSp)
      end do
    end do

    if (this%tMonitorValid) then
      allocate(this%localValidAtToGlobalSp(this%nValidDatapoints))
      do iGeo = 1, this%nValidDatapoints
        allocate(this%localValidAtToGlobalSp(iGeo)%array(this%validGeos(iGeo)%nAtom))
        do iAtom = 1, this%validGeos(iGeo)%nAtom
          localSp = this%validGeos(iGeo)%species(iAtom)
          this%localValidAtToGlobalSp(iGeo)%array(iAtom) =&
              & this%localValidSpToGlobalSp(iGeo)%array(localSp)
        end do
      end do
    end if

  end subroutine TData_getIndexMappings


  !> Check input consistency by evaluating several assertions.
  subroutine TProgramVariables_checkInputConsistency(this)

    !> instance containing program variables
    class(TProgramVariables), target :: this

    write(stdOut, '(A)', advance='no') 'Checking Input Consistency...'

    @:ASSERT(size(this%arch%hidden) > 0)
    @:ASSERT(minval(this%arch%hidden) > 0)

    @:ASSERT(size(this%arch%allDims) > 2)
    @:ASSERT(minval(this%arch%allDims) > 0)

    @:ASSERT(this%arch%nHiddenLayer >= 0)
    select case (tolower(this%option%mode))
    case('train')
      @:ASSERT(this%arch%nSubNnParams > 0)
    end select
    @:ASSERT(tolower(this%arch%type) == 'bpnn')

    associate (descriptor => this%arch%activation)
      @:ASSERT(descriptor == 'gaussian' .or. descriptor == 'relu' .or. descriptor == 'sigmoid'&
          & .or. descriptor == 'heaviside' .or. descriptor == 'tanh' .or. descriptor == 'linear')
    end associate

    if (tolower(this%option%mode) == 'train') then
      associate (descriptor => this%train%lossType)
        @:ASSERT(descriptor == 'mse' .or. descriptor == 'rms' .or. descriptor == 'mae'&
            & .or. descriptor == 'msle')
      end associate
    end if

    if (tolower(this%option%mode) == 'train') then
      @:ASSERT(this%train%nTrainIt >= 0)
      @:ASSERT(this%train%nPrintOut >= 1)
      @:ASSERT(this%train%nSaveNet >= 1)
    end if

    if ((tolower(this%option%mode) == 'train') .and. (.not. this%option%tReadNetStats)) then
      @:ASSERT(this%mapping%nRadial >= 0)
      @:ASSERT(this%mapping%nAngular >= 0)
      @:ASSERT(minval([this%mapping%nRadial, this%mapping%nAngular]) > 0)
      @:ASSERT(this%mapping%rCut > 0.0_dp)
    end if

    @:ASSERT(this%data%nDatapoints >= 1)
    @:ASSERT(this%data%nSpecies >= 1)
    @:ASSERT(this%data%nExtFeatures >= 0)
    @:ASSERT(size(this%data%geos) >= 1)
    @:ASSERT(size(this%data%netstatNames) == this%data%nSpecies)
    @:ASSERT(size(this%data%globalSpNames) == this%data%nSpecies)
    @:ASSERT(minval(this%data%weights) >= 1)
    @:ASSERT(size(this%data%weights) == size(this%data%geos))
    @:ASSERT(size(this%data%localSpToGlobalSp) == size(this%data%geos))
    @:ASSERT(size(this%data%localAtToGlobalSp) == size(this%data%geos))

    if (this%data%tContiguous) then
      @:ASSERT(size(this%data%datapaths) == 1)
    else
      @:ASSERT(size(this%data%datapaths) == size(this%data%geos))
    end if

    if (this%ext%tExtFeatures) then
      @:ASSERT(this%ext%nExtFeatures >= 0)
      @:ASSERT(this%ext%nExtFeatures == size(this%ext%extFeaturesInd))
      @:ASSERT(minval(this%ext%extFeaturesInd) >= 0)
    end if

    if ((tolower(this%option%mode) == 'train') .and. this%data%tMonitorValid) then
      @:ASSERT(this%data%nTargets == this%data%nValidTargets)
      @:ASSERT(this%data%nValidDatapoints >= 1)
      @:ASSERT(size(this%data%validGeos) >= 1)
      @:ASSERT(size(this%data%localValidSpToGlobalSp) == size(this%data%validGeos))
      @:ASSERT(size(this%data%localValidAtToGlobalSp) == size(this%data%validGeos))
      if (this%data%tValidContiguous) then
        @:ASSERT(size(this%data%validpaths) == 1)
      else
        @:ASSERT(size(this%data%validpaths) == size(this%data%validGeos))
      end if
    end if

    @:ASSERT(this%option%seed >= 0)

    associate (descriptor => this%option%mode)
      @:ASSERT(descriptor == 'train' .or. descriptor == 'validate' .or. descriptor == 'predict')
    end associate

    write(stdOut, '(A)') 'passed'

  end subroutine TProgramVariables_checkInputConsistency


  !> Prints the greeting message of Fortnet to standard output.
  subroutine printFortnetHeader(release, year)

    !> release version of the code
    character(len=*), intent(in) :: release

    !> release year
    integer, intent(in) :: year

    character, parameter :: vBar = '|'
    character, parameter :: hBar = '='
    integer, parameter :: headerWidth = 80

    write(stdOut, '(3A)') vBar, repeat(hBar, headerWidth - 2), vBar

    write(stdOut, '(5A)') vBar, '  Fortnet - A BPNN Implementation, Version ', trim(release),&
        & repeat(' ', 32), vBar
    write(stdOut, '(3A)') vBar, repeat(' ', headerWidth - 2), vBar
    write(stdOut, '(2A,I0,3A)') vBar, '  Copyright (C) 2020 - ', year,&
        & '  T. W. van der Heide', repeat(' ', 30), vBar
    write(stdOut, '(3A,/)') vBar, repeat(hBar, headerWidth - 2), vBar

  end subroutine printFortnetHeader


  !> Prints date and time information to standard output.
  subroutine printDateAndTime()

    !> Current date
    character(len=8) :: date

    !> Current time
    character(len=10) :: time

    !> Current time zone
    character(len=5) :: zone

    call date_and_time(date=date, time=time, zone=zone)

    write(stdout, '(6A)') 'date: ', date(7:8), '.', date(5:6), '.', date(1:4)
    write(stdout, '(8A)') 'time: ', time(1:2), ':', time(3:4), ':', date(5:6), ', ', zone

  end subroutine printDateAndTime

end module fnet_initprogram
