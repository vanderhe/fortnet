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

#:if WITH_MPI
  public :: syncFeatures
#:endif

  public :: TProgramVariables, TProgramVariables_init
  public :: TFeatures, TFeatures_init, TFeatures_collect
  public :: TArch, TData, TEnv, TExternal
  public :: initOptimizer, readAcsfFromFile


  !> data type containing variables of the Network block
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

    !> architecture type (currently only the BPNN topology is available)
    character(len=:), allocatable :: type

  end type TArch


  !> data type containing variables of the Mapping block
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


  !> data type containing the collected features
  type TFeatures

    !> pointer to the different input features for training
    type(TRealArray2D), allocatable :: features(:)

    !> pointer to the different input features for validation
    type(TRealArray2D), allocatable :: validFeatures(:)

  end type TFeatures


  !> data type containing variables of the Training block
  type TTraining

    !> integer ID of specified optimizertExtFeatures
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

    !> procedure, pointing to the choosed loss function
    procedure(lossFunc), pointer, nopass :: loss => null()

    !> gradient threshold where to stop the training, if given
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


  !> data type containing variables of the Data block
  type TData

    !> prefix of paths to saved species networks
    character(len=1024) :: prefix

    !> suffix of paths to saved species networks
    character(len=1024) :: suffix

    !> path to file containing paths to dataset geometries/targets
    character(len=:), allocatable :: datapath

    !> path to file containing paths to validation geometries/targets
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

    procedure :: getTargetMeansAndVariances => TData_getTargetMeansAndVariances
    procedure :: applyZscore => TData_applyZscore

  end type TData


  !> data type containing variables of the Option block
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


  !> data type containing variables of the Option block
  type TOption

    !> wether to resume from existing netstat files on disk
    logical :: tReadNetStats

    !> mode of current run (train, validate, predict)
    character(len=:), allocatable :: mode

    !> (user defined) random seed of the run
    integer :: seed

  end type TOption


  !> data type containing the program variables
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
    procedure :: getIndexMappings => TProgramVariables_getIndexMappings

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

  !> initialise program variables
  subroutine TProgramVariables_init(this)

    !> container of program variables
    type(TProgramVariables), intent(inout) :: this

    !> input tree node pointers
    type(fnode), pointer :: root, hsdTree, tmp

    !> string buffer instance
    type(string) :: strBuffer

    !> input version number
    integer :: inputVersion

    ! write header
    call printFortnetHeader(version, copyrightYear)
    call printDateAndTime()

    ! read in input file as HSD
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
    call readOptions(this, tmp)

    ! read data informations
    call getChild(root, 'Data', tmp)
    call readPrec(this, tmp)
    call readData(this, tmp)

    ! if present, read external features and identifiers
    call readExternal(this, root)

    select case(this%option%mode)

    case ('train')

      if (this%option%tReadNetStats) then

        ! in the case of resumed training, read network informations from netstat files
        call readFromNetstats(this)

      else

        ! read network informations
        call getChildValue(root, 'Network', tmp)
        call getNodeName(tmp, strBuffer)
        call readNetwork(this, tmp, trim(char(strBuffer)))

        ! read mapping informations
        call getChildValue(root, 'Mapping', tmp)
        call getNodeName(tmp, strBuffer)
        call readMapping(this, tmp, trim(char(strBuffer)))

      end if

      ! read training informations
      call getChildValue(root, 'Training', tmp)
      call getNodeName(tmp, strBuffer)
      call readTraining(this, tmp, trim(char(strBuffer)))

    case ('validate', 'predict')

      ! read informations from acsf and netstat files
      call readFromNetstats(this)

    end select

    ! calculate the number of targets per sub-network parameters
    call calcTargetsPerParam(this)

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


  !> initialize ACSF parameters from given file
  subroutine readAcsfFromFile(this, filename)

    !> container of program variables
    class(TProgramVariables), intent(inout) :: this

    !> filename or path to save acsf parameters to
    character(len=*), intent(in) :: filename

    this%mapping%type = 'acsf'

    write(stdOut, '(A)', advance='no') 'reading ACSF from file...'
    call this%acsf%fromFile(filename, this%data%globalSpNames)
    @:ASSERT(size(this%acsf%speciesIds) == this%data%nSpecies)
    write(stdout, '(A,/)') 'done'    

  end subroutine readAcsfFromFile


  !> initialize networks from netstat files
  subroutine readFromNetstats(this)

    !> container of program variables
    type(TProgramVariables), intent(inout) :: this

    !> unique fileunit
    integer :: fd

    !> auxiliary variables
    integer :: nLayers, iSpecies
    character(len=50) :: archType, targetType
    character(len=50) :: spName
    character(len=1024) :: activation

    do iSpecies = 1, this%data%nSpecies

      open(newunit=fd, file=this%data%netstatNames(iSpecies), form='formatted', status='old',&
          & action='read')

      read(fd, *) archType, targetType
      this%arch%type = tolower(trim(archType))
      if (tolower(trim(targetType)) == 'atomic') then
        this%data%tAtomicTargets = .true.
      elseif (tolower(trim(targetType)) == 'global') then
        this%data%tAtomicTargets = .false.
      else
        call error("Unrecognized target type in file '" // this%data%netstatNames(iSpecies) // "'.")
      end if

      read(fd, *) spName
      @:ASSERT(tolower(this%data%globalSpNames(iSpecies)) == tolower(spName))

      read(fd, *) nLayers
      @:ASSERT(nLayers > 2)
      if (allocated(this%arch%allDims)) deallocate(this%arch%allDims)
      allocate(this%arch%allDims(nLayers))
      read(fd, *) this%arch%allDims

      this%arch%hidden = this%arch%allDims(2:size(this%arch%allDims) - 1)
      this%arch%nHiddenLayer = size(this%arch%hidden)

      read(fd, *) activation
      this%arch%activation = tolower(trim(activation))

      close(fd)

    end do

    this%arch%nSubNnParams = getNumberOfParameters(this%arch%allDims)

  end subroutine readFromNetstats


  !> interpret the Network block
  subroutine readNetwork(this, node, case)

    !> container of program variables
    type(TProgramVariables), intent(inout) :: this

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

      this%arch%type = 'bpnn'

      call init(integerList)
      call getChildValue(node, 'Hidden', integerList)
      allocate(this%arch%hidden(len(integerList)))
      call asArray(integerList, this%arch%hidden)
      call destruct(integerList)
      this%arch%nHiddenLayer = size(this%arch%hidden)

      call getChildValue(node, 'Activation', strBuffer)
      this%arch%activation = tolower(trim(unquote(char(strBuffer))))

    case default

      call detailedError(node, 'Invalid network type')

    end select

  end subroutine readNetwork


  !> interpret the Mapping block
  subroutine readMapping(this, node, case)

    !> container of program variables
    type(TProgramVariables), intent(inout) :: this

    !> node containig the information
    type(fnode), intent(in), pointer :: node

    !> type of neural network
    character(len=*), intent(in) :: case

    call getChildValue(node, 'Standardization', this%mapping%tStandardize, .true.)

    select case (tolower(case))

    case ('acsf')

      call getChildValue(node, 'NRadial', this%mapping%nRadial)
      call getChildValue(node, 'NAngular', this%mapping%nAngular)
      call getChildValue(node, 'RCut', this%mapping%rCut)

      ! convert Angstrom to Bohr
      this%mapping%rCut = this%mapping%rCut * AA__Bohr

    case default

      call detailedError(node, 'Invalid mapping type')

    end select

    this%arch%allDims = [this%mapping%nRadial + this%mapping%nAngular + this%ext%nExtFeatures,&
        & this%arch%hidden, this%data%nTargets]

    this%arch%nSubNnParams = getNumberOfParameters(this%arch%allDims)

  end subroutine readMapping


  !> interpret the External block
  subroutine readExternal(this, node)

    !> container of program variables
    type(TProgramVariables), intent(inout) :: this

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
    if (this%data%tContiguous) then
      call readHSDAsXML(this%data%datapaths(1), xml)
      call getChild(xml, 'fnetdata', rootNode)
      call inquireFeatures(rootNode, this%data%nExtFeatures)
      call destroyNode(xml)
    else
      this%data%nExtFeatures = 0
      nTmpExtFeatures = this%data%nExtFeatures
      do iSys = 1, this%data%nDatapoints
        filename = trim(this%data%datapaths(iSys)) // '/' // fnetdataFile
        call readHSDAsXML(filename, xml)
        call getChild(xml, 'fnetdata', rootNode)
        call inquireFeatures(rootNode, this%data%nExtFeatures)
        if ((nTmpExtFeatures /= this%data%nExtFeatures) .and. (iSys /= 1)) then
          call error('Inconsistency in number of external features of dataset found.')
        end if
        nTmpExtFeatures = this%data%nExtFeatures
        call destroyNode(xml)
      end do
    end if

    ! inquire external validation dataset features
    if (this%data%tMonitorValid) then
      if (this%data%tValidContiguous) then
        nTmpExtFeatures = this%data%nExtFeatures
        call readHSDAsXML(this%data%validpaths(1), xml)
        call getChild(xml, 'fnetdata', rootNode)
        call inquireFeatures(rootNode, this%data%nExtFeatures)
        call destroyNode(xml)
        if (nTmpExtFeatures /= this%data%nExtFeatures) then
          call error('Inconsistency in number of external features between datasets found.')
        end if
        else
          nTmpExtFeatures = this%data%nExtFeatures
          do iSys = 1, this%data%nValidDatapoints
            filename = trim(this%data%validpaths(iSys)) // '/' // fnetdataFile
            call readHSDAsXML(filename, xml)
            call getChild(xml, 'fnetdata', rootNode)
            call inquireFeatures(rootNode, this%data%nExtFeatures)
            if (nTmpExtFeatures /= this%data%nExtFeatures) then
              call error('Inconsistency in number of external features between datasets found.')
            end if
            nTmpExtFeatures = this%data%nExtFeatures
            call destroyNode(xml)
          end do
        end if
      end if

    call getChild(node, 'External', extnode, requested=.false.)

    if (associated(extnode)) then

      ! read species identifier
      allocate(this%ext%speciesIds(this%data%nSpecies))
      call getChild(extnode, 'SpeciesID', child1, requested=.false.)
      if (associated(child1)) then
        do iSp = 1, this%data%nSpecies
          call getChildValue(child1, this%data%globalSpNames(iSp), this%ext%speciesIds(iSp),&
              & 1.0_dp)
          if (this%ext%speciesIds(iSp) <= 0.0_dp) then
            call warning("Obtained potentially dangerous '" // trim(this%data%globalSpNames(iSp))&
                & // "' species" // NEW_LINE('A') //&
                & '   identifier equal or below zero, watch out.')
          end if
        end do
      else
        this%ext%speciesIds(:) = 1.0_dp
      end if

      call getChild(extnode, 'Features', featuresnode, requested=.false.)
      if (associated(featuresnode)) then
        call getChildValue(featuresnode, '', buffer, child=child, multiple=.true.)
        call convRangeToInt(char(buffer), featuresnode, this%ext%extFeaturesInd,&
            & this%data%nExtFeatures)
        call setChildValue(child, '', this%ext%extFeaturesInd, replace=.true.)
        this%ext%nExtFeatures = size(this%ext%extFeaturesInd)
        if (size(this%ext%extFeaturesInd) > 0) then
          this%ext%tExtFeatures = .true.
        else
          this%ext%tExtFeatures = .false.
          this%ext%nExtFeatures = 0
        end if
      else
        this%ext%tExtFeatures = .false.
        this%ext%nExtFeatures = 0
      end if

      ! read external features from generic fnetdata.xml file(s)
      if (this%data%tContiguous .and. this%ext%tExtFeatures) then
        call readHSDAsXML(this%data%datapaths(1), xml)
        call getChild(xml, 'fnetdata', rootNode)
        call readContiguousFnetdataFeatures(rootNode, this%data%geos, this%ext%extFeatures,&
            & inds=this%ext%extFeaturesInd)
        call destroyNode(xml)
      elseif ((.not. this%data%tContiguous) .and. this%ext%tExtFeatures) then
        allocate(this%ext%extFeatures(this%data%nDatapoints))
        do iSys = 1, this%data%nDatapoints
          filename = trim(this%data%datapaths(iSys)) // '/' // fnetdataFile
          call readHSDAsXML(filename, xml)
          call getChild(xml, 'fnetdata', rootNode)
          call readFnetdataFeatures(rootNode, this%data%geos(iSys)%nAtom,&
              & this%ext%extFeatures(iSys)%array, inds=this%ext%extFeaturesInd)
          call destroyNode(xml)
        end do
      end if

      ! if present, read external validation features from generic fnetdata.xml file(s)
      if (this%data%tMonitorValid .and. this%ext%tExtFeatures .and. this%data%tValidContiguous) then
        call readHSDAsXML(this%data%validpaths(1), xml)
        call getChild(xml, 'fnetdata', rootNode)
        call readContiguousFnetdataFeatures(rootNode, this%data%validGeos,&
            & this%ext%extValidFeatures, inds=this%ext%extFeaturesInd)
        call destroyNode(xml)
      elseif (this%data%tMonitorValid .and. this%ext%tExtFeatures .and.&
          & (.not. this%data%tValidContiguous)) then
        allocate(this%ext%extValidFeatures(this%data%nValidDatapoints))
        do iSys = 1, this%data%nValidDatapoints
          filename = trim(this%data%validpaths(iSys)) // '/fnetdata.xml'
          call readHSDAsXML(filename, xml)
          call getChild(xml, 'fnetdata', rootNode)
          call readFnetdataFeatures(rootNode, this%data%validGeos(iSys)%nAtom,&
              & this%ext%extValidFeatures(iSys)%array, inds=this%ext%extFeaturesInd)
          call destroyNode(xml)
        end do
      end if

      ! read atom identifier
      call getChild(extnode, 'AtomID', child2, requested=.false.)
      if (associated(child2)) then
        call getChildValue(child2, '', this%ext%atomIdIndex)
        if (this%ext%atomIdIndex <= 0) then
          call warning('Atom identifier index less or equal zero specified, will be ignored.')
        elseif (this%ext%atomIdIndex > this%data%nExtFeatures) then
          call error('Atom identifier index exceeds number of external features in dataset.')
        end if
      else
        this%ext%atomIdIndex = 0
      end if

    else

      allocate(this%ext%speciesIds(this%data%nSpecies))
      this%ext%speciesIds(:) = 1.0_dp
      this%ext%atomIdIndex = 0
      this%ext%tExtFeatures = .false.
      this%ext%nExtFeatures = 0

    end if

    if (this%ext%atomIdIndex > 0) then
      if (this%data%tContiguous) then
        call readHSDAsXML(this%data%datapaths(1), xml)
        call getChild(xml, 'fnetdata', rootNode)
        call readContiguousFnetdataAtomIdentifier(rootNode, this%data%geos,&
            & this%ext%atomIdIndex, this%ext%atomIds)
        call destroyNode(xml)
      else
        allocate(this%ext%atomIds(this%data%nDatapoints))
        do iSys = 1, this%data%nDatapoints
          filename = trim(this%data%datapaths(iSys)) // '/' // fnetdataFile
          call readHSDAsXML(filename, xml)
          call getChild(xml, 'fnetdata', rootNode)
          call readFnetdataAtomIdentifier(rootNode, this%data%geos(iSys)%nAtom,&
              & this%ext%atomIdIndex, this%ext%atomIds(iSys)%array)
          call destroyNode(xml)
        end do
      end if
    else
      allocate(this%ext%atomIds(this%data%nDatapoints))
      do iSys = 1, this%data%nDatapoints
        allocate(this%ext%atomIds(iSys)%array(this%data%geos(iSys)%nAtom))
        this%ext%atomIds(iSys)%array(:) = 1.0_dp
      end do
    end if

    if (this%data%tMonitorValid) then
      if (this%ext%atomIdIndex > 0) then
        if (this%data%tValidContiguous) then
          call readHSDAsXML(this%data%validpaths(1), xml)
          call getChild(xml, 'fnetdata', rootNode)
          call readContiguousFnetdataAtomIdentifier(rootNode, this%data%validGeos,&
              & this%ext%atomIdIndex, this%ext%validAtomIds)
          call destroyNode(xml)
        else
          allocate(this%ext%validAtomIds(this%data%nValidDatapoints))
          do iSys = 1, this%data%nValidDatapoints
            filename = trim(this%data%validpaths(iSys)) // '/' // fnetdataFile
            call readHSDAsXML(filename, xml)
            call getChild(xml, 'fnetdata', rootNode)
            call readFnetdataAtomIdentifier(rootNode, this%data%validGeos(iSys)%nAtom,&
                & this%ext%atomIdIndex, this%ext%validAtomIds(iSys)%array)
            call destroyNode(xml)
          end do
        end if
      else
        allocate(this%ext%validAtomIds(this%data%nValidDatapoints))
        do iSys = 1, this%data%nValidDatapoints
          allocate(this%ext%validAtomIds(iSys)%array(this%data%validGeos(iSys)%nAtom))
          this%ext%validAtomIds(iSys)%array(:) = 1.0_dp
        end do
      end if
    end if

  end subroutine readExternal


  !> interpret the Training block
  subroutine readTraining(this, node, case)

    !> container of program variables
    type(TProgramVariables), intent(inout) :: this

    !> node containig the information
    type(fnode), intent(in), pointer :: node

    !> type of neural network
    character(len=*), intent(in) :: case

    !> string buffer instance
    type(string) :: strBuffer

    select case (tolower(case))

    case ('sd')

      this%train%iOptimizer = optimizerTypes%steepDesc
      call getChildValue(node, 'LearningRate', this%train%learningRate, 0.01_dp)

    case ('cg')

      this%train%iOptimizer = optimizerTypes%conjGrad

    case ('lbfgs')

      this%train%iOptimizer = optimizerTypes%lbfgs
      call getChildValue(node, 'MaxForQNDisplacement', this%train%maxForQNDisplacement, .false.)
      call getChildValue(node, 'Linemin', this%train%tLinesearch, .true.)
      call getChildValue(node, 'Memory', this%train%mem, 1000)

    case ('fire')

      this%train%iOptimizer = optimizerTypes%fire

    case default

      call detailedError(node, 'Invalid training algorithm')

    end select

    call getChildValue(node, 'MinDisplacement', this%train%minDisplacement, 1e-06_dp)
    call getChildValue(node, 'MaxDisplacement', this%train%maxDisplacement, 1e+04_dp)

    call getChildValue(node, 'NIterations', this%train%nTrainIt, huge(0))
    call getChildValue(node, 'Threshold', this%train%threshold, tiny(0.0_dp))

    call getChildValue(node, 'NPrintout', this%train%nPrintOut, 10)
    call getChildValue(node, 'NSaveNet', this%train%nSaveNet, 100)

    call getChildValue(node, 'Loss', strBuffer, 'rms')
    call this%train%setLossFunc(tolower(trim(unquote(char(strBuffer)))))

  end subroutine readTraining


  !> interpret the Options block
  subroutine readOptions(this, node)

    !> container of program variables
    type(TProgramVariables), intent(inout) :: this

    !> node containig the information
    type(fnode), pointer :: node

    !> string buffer instance
    type(string) :: strBuffer

    !> temporary random integer
    integer :: tmpIntSeed

    !> auxiliary variable
    real(dp) :: tmpRealSeed

    call getChildValue(node, 'ReadNetStats', this%option%tReadNetStats, .false.)

    call getChildValue(node, 'Mode', strBuffer)
    this%option%mode = tolower(unquote(char(strBuffer)))

    if ((trim(this%option%mode) == 'validate' .or. trim(this%option%mode) == 'predict') .and.&
        & (this%option%tReadNetStats .eqv. .false.)) then
      write(stdout, '(A)') ''
      call warning('Running in validation or prediction mode without initialising from'&
          & //NEW_LINE('A')//'   existing netstat files is not possible. Overwriting user input...')
      write(stdout, '(A)') ''
      this%option%tReadNetStats = .true.
    end if

    call random_number(tmpRealSeed)
    tmpIntSeed = floor(tmpRealSeed * huge(0) + tiny(0.0))

    call getChildValue(node, 'RandomSeed', this%option%seed, tmpIntSeed)

    if (this%option%seed < 0) then
      call detailedError(node, 'Random seed must be greater or equal zero.')
    end if

    call init(this%rndGen, luxlev=3, initSeed=this%option%seed)

  end subroutine readOptions


  !> interpret preconditioning of the Data block
  subroutine readPrec(this, node)

    !> container of program variables
    type(TProgramVariables), intent(inout) :: this

    !> node containig the information
    type(fnode), pointer :: node, child

    !> true, if preconditioning file is in place
    logical :: tExist

    call getChildValue(node, 'Standardization', this%data%tZscore, .false., child=child)

    inquire(file=precFile, exist=tExist)

    if (.not. this%data%tZscore .and. tExist) then
      call warning('User input manually deactivated target standardization,'&
          & //NEW_LINE('A')// '   but preconditioning parameters are present.')
    end if

    call readPreconditioning(precFile, this%option%tReadNetStats, this%data%tZscore,&
        & this%data%nTargets, this%data%zPrec)

  end subroutine readPrec


  !> read datapoint paths from file
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


  !> interpret the Data block
  subroutine readData(this, node)

    !> container of program variables
    type(TProgramVariables), intent(inout) :: this

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
    this%data%datapath = trim(unquote(char(strBuffer)))

    if (this%data%datapath(len(this%data%datapath)-len(fnetdataFile)+1:len(this%data%datapath))&
        & == fnetdataFile) then
      this%data%tContiguous = .true.
      allocate(this%data%datapaths(1))
      this%data%datapaths(1) = this%data%datapath
    else
      this%data%tContiguous = .false.
      call readDatapathsFromFile(this%data%datapath, node, this%data%nDatapoints,&
          & this%data%datapaths)
    end if

    call getChildValue(node, 'Validset', strBuffer, default='')

    select case (char(strBuffer))
    case ('')
      this%data%tMonitorValid = .false.
    case ('none')
      this%data%tMonitorValid = .false.
    case default
      this%data%tMonitorValid = .true.
      this%data%validpath = trim(unquote(char(strBuffer)))
      if (this%data%validpath(len(this%data%validpath)-len(fnetvdataFile)+1&
          & :len(this%data%validpath)) == fnetvdataFile) then
        this%data%tValidContiguous = .true.
        allocate(this%data%validpaths(1))
        this%data%validpaths(1) = this%data%validpath
      else
        this%data%tValidContiguous = .false.
        call readDatapathsFromFile(this%data%validpath, node, this%data%nValidDatapoints,&
            & this%data%validpaths)
      end if
    end select

    ! read geometries from generic fnetdata.xml file(s)
    if (this%data%tContiguous) then
      call readHSDAsXML(this%data%datapaths(1), xml)
      call getChild(xml, 'fnetdata', rootNode)
      call readContiguousFnetdataWeights(rootNode, this%data%weights)
      call readContiguousFnetdataGeometries(rootNode, this%data%geos)
      this%data%nDatapoints = size(this%data%geos)
      select case (this%option%mode)
      case('train', 'validate')
        call readContiguousFnetdataTargets(rootNode, this%data%geos, this%data%targets,&
            & this%data%nTargets, this%data%tAtomicTargets)
      end select
      call destroyNode(xml)
    else
      allocate(this%data%weights(this%data%nDatapoints))
      allocate(this%data%geos(this%data%nDatapoints))
      select case (this%option%mode)
      case('train', 'validate')
        allocate(this%data%targets(this%data%nDatapoints))
      end select
      do iSys = 1, this%data%nDatapoints
        filename = trim(this%data%datapaths(iSys)) // '/' // fnetdataFile
        call readHSDAsXML(filename, xml)
        call getChild(xml, 'fnetdata', rootNode)
        call readFnetdataWeight(rootNode, this%data%weights(iSys))
        call readFnetdataGeometry(rootNode, this%data%geos(iSys))
        select case (this%option%mode)
        case('train', 'validate')
          call readFnetdataTargets(rootNode, this%data%geos(iSys)%nAtom,&
              & this%data%targets(iSys)%array, this%data%tAtomicTargets)
          this%data%nTargets = size(this%data%targets(iSys)%array, dim=1)
        end select
        call destroyNode(xml)
      end do
    end if

    this%data%nTotalAtoms = getTotalNumberOfAtoms(this%data%geos)

    ! if present, read validation geometries from generic fnetdata.xml file(s)
    if (this%data%tMonitorValid) then
      if (this%data%tValidContiguous) then
        call readHSDAsXML(this%data%validpaths(1), xml)
        call getChild(xml, 'fnetdata', rootNode)
        call readContiguousFnetdataGeometries(rootNode, this%data%validGeos)
        this%data%nValidDatapoints = size(this%data%validGeos)
        select case (this%option%mode)
        case('train', 'validate')
          call readContiguousFnetdataTargets(rootNode, this%data%validGeos, this%data%validTargets,&
              & this%data%nValidTargets, this%data%tAtomicTargets)
        end select
        call destroyNode(xml)
      else
        allocate(this%data%validGeos(this%data%nValidDatapoints))
        select case (this%option%mode)
        case('train', 'validate')
          allocate(this%data%validTargets(this%data%nValidDatapoints))
        end select
        do iSys = 1, this%data%nValidDatapoints
          filename = trim(this%data%validpaths(iSys)) // '/fnetdata.xml'
          call readHSDAsXML(filename, xml)
          call getChild(xml, 'fnetdata', rootNode)
          call readFnetdataGeometry(rootNode, this%data%validGeos(iSys))
          select case (this%option%mode)
          case('train', 'validate')
            call readFnetdataTargets(rootNode, this%data%validGeos(iSys)%nAtom,&
                & this%data%validTargets(iSys)%array, this%data%tAtomicTargets)
            this%data%nValidTargets = size(this%data%validTargets(iSys)%array, dim=1)
          end select
          call destroyNode(xml)
        end do
      end if
    end if

    ! apply z-score standardization, if desired
    if (this%data%tZscore) then
      if (.not. allocated(this%data%zPrec)) then
        call this%data%getTargetMeansAndVariances()
        call writePreconditioning(precFile, this%data%zPrec)
      end if
      call this%data%applyZscore()
    elseif (.not. this%option%mode == 'predict') then
      this%data%zTargets = this%data%targets
      if (this%data%tMonitorValid) then
        this%data%zValidTargets = this%data%validTargets
      end if
    end if

    ! establish convenient index mappings and global species list
    call this%getIndexMappings()
    allocate(this%data%netstatNames(this%data%nSpecies))

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

      do iSp = 1, this%data%nSpecies
        if (tLower) then
          elem = tolower(this%data%globalSpNames(iSp))
        else
          elem = this%data%globalSpNames(iSp)
        end if

        strTmp = trim(prefix) // trim(elem) // trim(suffix)
        this%data%netstatNames(iSp) = strTmp

        if (this%option%tReadNetStats) then
          inquire(file=strTmp, exist=tExist)
        end if

        if (this%option%tReadNetStats .and. (.not. tExist)) then
          call detailedError(value1, "Netstat file with generated name '" // trim(strTmp) //&
              & "' does not exist.")
        end if

      end do

    case default

      call setUnprocessed(value1)

      do iSp = 1, this%data%nSpecies

          call init(lStr)
          call getChildValue(child1, trim(this%data%globalSpNames(iSp)), lStr, child=child2)

          this%data%netstatNames(iSp) = trim(this%data%globalSpNames(iSp))

          if (len(lStr) /= this%data%nSpecies) then
            call detailedError(child2, "Incorrect number of netstat files")
          end if

          do ii = 1, len(lStr)
            call get(lStr, this%data%globalSpNames(iSp), ii)

            if (this%option%tReadNetStats) then
              inquire(file=strTmp, exist=tExist)
            end if

            if (this%option%tReadNetStats .and. (.not. tExist)) then
              call detailedError(child2, "Netstat file '" // trim(this%data%globalSpNames(iSp))&
                  & // "' does not exist'")
            end if
          end do

          call destruct(lStr)

        end do

      end select

  end subroutine readData


  !> calculates z-score standardization means and variances
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
  subroutine initOptimizer(this, initialPoints)

    !> container of program variables
    type(TProgramVariables), intent(inout) :: this

    !> starting points, shape: [nPoints, nStrucs]
    real(dp), intent(in) :: initialPoints(:,:)

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

    nValues = size(initialPoints, dim=1)
    nStrucs = size(initialPoints, dim=2)

    allocate(this%train%pOptimizer(nStrucs))

    select case (this%train%iOptimizer)

    case(optimizerTypes%steepDesc)

      allocate(wrapSteepDesc(nStrucs))
      allocate(weights(nValues))
      weights(:) = this%train%learningRate

      do iStruc = 1, nStrucs
        allocate(wrapSteepDesc(iStruc)%pSteepDesc)
        call init(wrapSteepDesc(iStruc)%pSteepDesc, nValues, this%train%threshold,&
            & this%train%maxDisplacement, weights)
        call init(wrapSteepDesc(iStruc)%pSteepDesc, this%train%pOptimizer(iStruc))
      end do

      deallocate(weights)

    case (optimizerTypes%conjGrad)

      allocate(wrapConjGrad(nStrucs))

      do iStruc = 1, nStrucs
        allocate(wrapConjGrad(iStruc)%pConjGrad)
        call init(wrapConjGrad(iStruc)%pConjGrad, nValues, this%train%threshold,&
            & this%train%maxDisplacement)
        call init(wrapConjGrad(iStruc)%pConjGrad, this%train%pOptimizer(iStruc))
      end do

    case (optimizerTypes%lbfgs)

      allocate(wrapLbfgs(nStrucs))

      do iStruc = 1, nStrucs
        allocate(wrapLbfgs(iStruc)%pLbfgs)
        call TLbfgs_init(wrapLbfgs(iStruc)%pLbfgs, nValues, this%train%threshold,&
            & this%train%minDisplacement, this%train%maxDisplacement, this%train%mem,&
            & this%train%tLinesearch, .false., this%train%maxForQNDisplacement)
        call init(wrapLbfgs(iStruc)%pLbfgs, this%train%pOptimizer(iStruc))
      end do

    case (optimizerTypes%fire)

      allocate(wrapFire(nStrucs))

      do iStruc = 1, nStrucs
        allocate(wrapFire(iStruc)%pFire)
        call TFire_init(wrapFire(iStruc)%pFire, nValues, this%train%threshold,&
            & this%train%maxDisplacement)
        call init(wrapFire(iStruc)%pFire, this%train%pOptimizer(iStruc))
      end do

    end select

    do iStruc = 1, nStrucs
      call reset(this%train%pOptimizer(iStruc), initialPoints(:, iStruc))
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
  subroutine syncFeatures(features, comm)

    !> collected features of data and mapping block
    type(TFeatures), intent(in) :: features

    !> mpi communicator with some additional information
    type(mpifx_comm), intent(in) :: comm

    !> auxiliary variable
    integer :: iGeo

    do iGeo = 1, size(features%features)
      call mpifx_bcast(comm, features%features(iGeo)%array)
    end do

    if (allocated(features%validFeatures)) then
      do iGeo = 1, size(features%validFeatures)
        call mpifx_bcast(comm, features%validFeatures(iGeo)%array)
      end do
    end if

  end subroutine syncFeatures
#:endif


  !> calculates the total number of network fitting parameters per sub-nn
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


  !> calculates the total number of atoms in a list of geometries
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


  !> calculates the number of targets per sub-network parameters
  subroutine calcTargetsPerParam(this)

    !> container of program variables
    class(TProgramVariables), intent(inout) :: this

    ! calculate number of targets per network parameter
    if (this%data%tAtomicTargets) then
      this%data%nTargetsPerParam = real(this%data%nTargets * this%data%nTotalAtoms, dp) /&
          & real(this%data%nSpecies * this%arch%nSubNnParams , dp)
    else
      this%data%nTargetsPerParam = real(this%data%nTargets * this%data%nDatapoints, dp) /&
          & real(this%data%nSpecies * this%arch%nSubNnParams, dp)
    end if

  end subroutine calcTargetsPerParam


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


  !> establish convenient index mappings
  subroutine TProgramVariables_getIndexMappings(this)

    !> container of program variables
    class(TProgramVariables), intent(inout) :: this

    !> auxiliary variables
    integer :: iGeo, iSp, iAtom, localSp, ind

    allocate(this%data%globalSpNames(1))
    this%data%globalSpNames(1) = this%data%geos(1)%speciesNames(1)

    allocate(this%data%localSpToGlobalSp(this%data%nDatapoints))
    if (this%data%tMonitorValid) then

    end if

    ind = 1

    do iGeo = 1, this%data%nDatapoints
      allocate(this%data%localSpToGlobalSp(iGeo)%array(this%data%geos(iGeo)%nSpecies))
      do iSp = 1, this%data%geos(iGeo)%nSpecies
        if (.not. any(this%data%globalSpNames .eq. this%data%geos(iGeo)%speciesNames(iSp))) then
          this%data%globalSpNames = [this%data%globalSpNames,&
              & this%data%geos(iGeo)%speciesNames(iSp)]
          ind = ind + 1
        end if
        ! work around the findloc intrinsic
        ! this%data%localSpToGlobalSp(iGeo)%array(iSp) = findloc(this%data%globalSpNames,&
        !     & this%data%geos(iGeo)%speciesNames(iSp), dim=1)
        this%data%localSpToGlobalSp(iGeo)%array(iSp) = myFindloc(this%data%globalSpNames,&
            & this%data%geos(iGeo)%speciesNames(iSp))
      end do
    end do

    this%data%nSpecies = ind

    if (this%data%tMonitorValid) then
      allocate(this%data%localValidSpToGlobalSp(this%data%nValidDatapoints))
      do iGeo = 1, this%data%nValidDatapoints
        allocate(this%data%localValidSpToGlobalSp(iGeo)%array(this%data%validGeos(iGeo)%nSpecies))
        do iSp = 1, this%data%validGeos(iGeo)%nSpecies
          ! work around the findloc intrinsic
          ! this%data%localValidSpToGlobalSp(iGeo)%array(iSp) = findloc(this%data%globalSpNames,&
          !     & this%data%validGeos(iGeo)%speciesNames(iSp), dim=1)
          this%data%localValidSpToGlobalSp(iGeo)%array(iSp) = myFindloc(this%data%globalSpNames,&
              & this%data%validGeos(iGeo)%speciesNames(iSp))
        end do
      end do
    end if

    allocate(this%data%localAtToGlobalSp(this%data%nDatapoints))

    do iGeo = 1, this%data%nDatapoints
      allocate(this%data%localAtToGlobalSp(iGeo)%array(this%data%geos(iGeo)%nAtom))
      do iAtom = 1, this%data%geos(iGeo)%nAtom
        localSp = this%data%geos(iGeo)%species(iAtom)
        this%data%localAtToGlobalSp(iGeo)%array(iAtom) =&
            & this%data%localSpToGlobalSp(iGeo)%array(localSp)
      end do
    end do

    if (this%data%tMonitorValid) then
      allocate(this%data%localValidAtToGlobalSp(this%data%nValidDatapoints))
      do iGeo = 1, this%data%nValidDatapoints
        allocate(this%data%localValidAtToGlobalSp(iGeo)%array(this%data%validGeos(iGeo)%nAtom))
        do iAtom = 1, this%data%validGeos(iGeo)%nAtom
          localSp = this%data%validGeos(iGeo)%species(iAtom)
          this%data%localValidAtToGlobalSp(iGeo)%array(iAtom) =&
              & this%data%localValidSpToGlobalSp(iGeo)%array(localSp)
        end do
      end do
    end if

  end subroutine TProgramVariables_getIndexMappings


  !> several assertions to check input consistency
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

    if (tolower(this%option%mode) == 'train' .and. (.not. this%option%tReadNetStats)) then
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

    if (this%data%tMonitorValid) then
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


  !> prints the greeting message of Fortnet to stdout
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
