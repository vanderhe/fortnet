!--------------------------------------------------------------------------------------------------!
!  FORTNET: A Behler-Parrinello-Neural-Network Implementation                                      !
!  Copyright (C) 2020 - 2021  T. W. van der Heide                                                  !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains the routines for initialising Fortnet from an HSD input.
module fnet_initprogram

  use dftbp_xmlf90
  use dftbp_linkedlist
  use dftbp_assert
  use dftbp_hsdutils, only : getChild, getChildren, setChild, getChildValue, setChildValue
  use dftbp_hsdutils, only : convRangeToInt, detailedError
  use dftbp_hsdutils2, only : readHSDAsXML, setUnprocessed, warnUnprocessedNodes
  use dftbp_charmanip, only : unquote, tolower, i2c
  use dftbp_accuracy, only: dp
  use dftbp_message, only : error, warning
  use dftbp_steepdesc, only : init
  use dftbp_conjgrad, only : init
  use dftbp_lbfgs, only : TLbfgs_init
  use dftbp_fire, only : TFire_init
  use dftbp_ranlux, only : TRanlux, init
  use dftbp_constants, only : AA__Bohr
  use dftbp_globalenv, only : stdOut
  use dftbp_hsdparser, only : parseHSD, dumpHSD

  use fnet_fnetdata, only : inquireExtFeatures, inquireTargets, inquireStructures, readHdfDataset,&
      & TDataset, checkDatasetCompatibility
  use fnet_netstat, only : createNetstat, readSubnetArchitecture
  use fnet_acsf, only : TGFunction_init, TGFunctions
  use fnet_features, only : TFeatures, TFeaturesBlock
  use fnet_nestedtypes, only : TEnv, TWrapSteepDesc, TWrapConjGrad, TWrapLbfgs, TWrapFire
  use fnet_optimizers, only : TOptimizer, init, reset, optimizerTypes
  use fnet_loss, only : lossFunc, lossGradientFunc, reguFunc, maLoss, mapLoss,&
      & msLoss, rmsLoss, nullReguLoss, lassoReguLoss, ridgeReguLoss, elasticNetReguLoss,&
      & maGradients, mapGradients, msGradients, rmsGradients, TRegularizationBlock
  use fnet_intmanip, only : combinationsWithReplacement, factorial

#:if WITH_MPI
  use fnet_mpifx
  use fnet_fnetdata, only : syncDataset
#:endif

  implicit none

  private
  save

  public :: TProgramVariables, TProgramVariables_init
  public :: TTraining_initOptimizer
  public :: TNetworkBlock, TTrainingBlock, TDataBlock, TOptionBlock


  !> Data type containing variables of the Network block
  type TNetworkBlock

    !> number of hidden layers
    integer :: nHiddenLayer

    !> number of network paramaters (weights + biases) per sub-nn
    integer :: nSubNnParams

    !> number of nodes per hidden layer, expected shape: [nHiddenLayer]
    integer, allocatable :: hidden(:)

    !> number of nodes per layer, including in- and output, expected shape: [nHiddenLayer + 2]
    integer, allocatable :: allDims(:)

    !> architecture type (currently, only the BPNN topology is available)
    character(len=:), allocatable :: type

    !> type of activation functions
    character(len=:), allocatable :: activation

  end type TNetworkBlock


  !> Data type containing variables of the Training block
  type TTrainingBlock

    !> maximum number of training iterations
    integer :: nTrainIt

    !> printout loss/gradient information every nPrintOut steps
    integer :: nPrintOut

    !> save network status every nSaveNet steps
    integer :: nSaveNet

    !> integer ID of specified optimizert
    integer :: iOptimizer

    !> type of loss function to use during the training
    character(len=:), allocatable :: lossType

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

    !> wether a Knuth-shuffle should be applied to the gradient calculation of datapoints
    logical :: tShuffle

    !> true, if loss-based regularization is requested
    logical :: tRegularization

    !> data type containing variables of the Regularization block
    type(TRegularizationBlock) :: regu

  end type TTrainingBlock


  !> Data type containing variables of the Data block
  type TDataBlock

    !> true, if validation monitoring is desired
    logical :: tMonitorValid

    !> path to file containing paths to training datapoints
    character(len=:), allocatable :: trainpath

    !> path to file containing paths to validation datapoints
    character(len=:), allocatable :: validpath

    !> path to netstat file defining the program state
    character(len=:), allocatable :: netstatpath

  end type TDataBlock


  !> Data type containing variables of the Option block
  type TOptionBlock

    !> wether to resume from existing netstat files on disk
    logical :: tReadNetStats

    !> wether to write loss and gradients for all training iterations to disk
    logical :: tWriteIterTraj

    !> (user defined) random seed of the run
    integer :: seed

    !> mode of current run (train, validate, predict)
    character(len=:), allocatable :: mode

  end type TOptionBlock


  !> Data type containing training related information
  type TTraining

    !> general function optimizer
    type(TOptimizer), allocatable :: pOptimizer(:)

    !> procedure, pointing to the choosen loss function
    procedure(lossFunc), pointer, nopass :: loss => null()

    !> procedure, pointing to the choosen loss function gradient
    procedure(lossGradientFunc), pointer, nopass :: lossgrad => null()

    !> procedure, pointing to the choosen regularization function
    procedure(reguFunc), pointer, nopass :: reguLoss => null()

  contains

    procedure :: setLossFunc => TTraining_setLossFunc
    procedure :: setReguFunc => TTraining_setReguFunc

  end type TTraining


  !> Data type containing all variable parsed form HSD input
  type TInput

    !> data type containing variables of the Network block
    type(TNetworkBlock) :: network

    !> data type containing variables of the Features block
    type(TFeaturesBlock) :: features

    !> data type containing variables of the Training block
    type(TTrainingBlock) :: training

    !> data type containing variables of the Data block
    type(TDataBlock) :: data

    !> data type containing variables of the Option block
    type(TOptionBlock) :: option

  contains

    procedure :: checkInputConsistency => TInput_checkInputConsistency

  end type TInput


  !> Data type containing the program variables
  type TProgramVariables

    !> contains the interpreted input
    type(TInput) :: inp

    !> contains mpi communicator, if compiled with mpi enabled
    type(TEnv) :: env

    !> training dataset representation
    type(TDataset) :: trainDataset

    !> validation dataset representation
    type(TDataset) :: validDataset

    !> training related information and functions
    type(TTraining) :: train

    !> collected features, including external and ACSF features
    type(TFeatures) :: features

    !> luxury pseudorandom generator instance
    type(TRanlux) :: rndGen

  end type TProgramVariables


  !> program version
  character(len=*), parameter :: version = '0.3'

  !> copyright year
  integer, parameter :: copyrightYear = 2021

  !> root node name of the input tree
  character(len=*), parameter :: rootTag = 'fortnet'

  !> input file name
  character(len=*), parameter :: hsdInput = 'fortnet_in.hsd'

  !> parsed output name
  character(len=*), parameter :: hsdParsedInput = 'fortnet_pin.hsd'

  !> version of the input document
  integer, parameter :: parserVersion = 1


contains

  !> Initializes program variables depending on the user input.
  subroutine TProgramVariables_init(this)

    !> container of program variables
    type(TProgramVariables), intent(inout) :: this

    !> input tree node pointers
    type(fnode), pointer :: root, hsdTree, tmp

    !> string buffer instance
    type(string) :: strBuffer

    !> input version number
    integer :: inputVersion

    !> true, if current mpi process is the lead
    logical :: tLead

  #:if WITH_MPI
    tLead = this%env%globalMpiComm%lead
  #:else
    tLead = .true.
  #:endif

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

    ! read options block
    call getChild(root, 'Options', tmp)
    call readOptionBlock(this%inp%option, tmp)

    ! define state of luxury random number generator
    call init(this%rndGen, luxlev=3, initSeed=this%inp%option%seed)

    ! read data informations
    call getChild(root, 'Data', tmp)
    call readDataBlock(this%inp%data, tmp, this%inp%option%mode)
    if (tLead) then
      call readHdfDataset(this%inp%data%trainpath, this%trainDataset)
      if (this%inp%data%tMonitorValid) then
        call readHdfDataset(this%inp%data%validpath, this%validDataset)
        call checkDatasetCompatibility(this%trainDataset, this%validDataset)
      end if
    end if

  #:if WITH_MPI
    call syncDataset(this%trainDataset, this%env%globalMpiComm)
    if (this%inp%data%tMonitorValid) then
      call syncDataset(this%validDataset, this%env%globalMpiComm)
    end if
  #:endif

    select case(this%inp%option%mode)

    case ('train')

      if (this%inp%option%tReadNetStats) then

        ! in the case of resumed training, read BPNN from netstat file
        call readFromNetstat(this%inp%data%netstatpath, this%inp%network)

      else

        ! read feature information
        call getChild(root, 'Features', tmp)
        call readFeaturesBlock(this%inp%features, tmp, this%trainDataset%nExtFeatures)
        call processFeaturesBlock(this%inp%features, this%trainDataset%atomicNumbers)

        ! read network informations
        call getChildValue(root, 'Network', tmp)
        call getNodeName(tmp, strBuffer)
        call readNetworkBlock(this%inp%network, tmp, trim(char(strBuffer)),&
            & this%inp%features%nFeatures, this%trainDataset%nTargets)

        ! create empty netstat file
        if (tLead) then
          call createNetstat(this%inp%data%netstatpath)
        end if

      end if

      ! read training informations
      call getChildValue(root, 'Training', tmp)
      call getNodeName(tmp, strBuffer)
      call readTrainingBlock(this%inp%training, tmp, trim(char(strBuffer)))
      call TTraining_init(this%train, this%inp%training)

    case ('validate', 'predict')

      ! read information from netstat file
      call readFromNetstat(this%inp%data%netstatpath, this%inp%network)

    end select

    ! calculate the number of targets per sub-network parameter
    call calcTargetsPerParam(this%trainDataset, this%inp%network%nSubNnParams)

    ! issue warning about unprocessed nodes
    call warnUnprocessedNodes(root, .true.)

    ! check the rudimentary plausibility of the imported parameters
    call this%inp%checkInputConsistency()

    ! finish parsing, dump parsed and processed input
    call dumpHSD(hsdTree, hsdParsedInput)

    write(stdout, '(A,/)') "Processed input written as HSD to '" // hsdParsedInput //"'"
    write(stdout, '(A,/)') repeat('-', 80)

    call destroyNode(tmp)
    call destroyNode(hsdTree)

  end subroutine TProgramVariables_init


  !> Calculates the number of targets per sub-network parameter.
  subroutine calcTargetsPerParam(dataset, nSubNnParams)

    !> representation of a dataset
    type(TDataset), intent(inout) :: dataset

    !> number of network paramaters (weights + biases) per sub-nn
    integer, intent(in) :: nSubNnParams

    ! calculate number of targets per network parameter
    if (dataset%tAtomicTargets) then
      dataset%nTargetsPerParam = real(dataset%nTargets * dataset%nTotalAtoms, dp) /&
          & real(dataset%nSpecies * nSubNnParams , dp)
    else
      dataset%nTargetsPerParam = real(dataset%nTargets * dataset%nDatapoints, dp) /&
          & real(dataset%nSpecies * nSubNnParams, dp)
    end if

  end subroutine calcTargetsPerParam


  !> Reads neural sub-network architecture from netstat file.
  subroutine readFromNetstat(fname, network)

    !> name of the netstat file
    character(len=*), intent(in) :: fname

    !> data type containing variables of the Network block
    type(TNetworkBlock), intent(out) :: network

    call readSubnetArchitecture(fname, network%type, network%activation, network%allDims)

    network%hidden = network%allDims(2:size(network%allDims)-1)
    network%nHiddenLayer = size(network%hidden)

    network%nSubNnParams = getNumberOfParameters(network%allDims)

  end subroutine readFromNetstat


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


  !> Interprets the Options HSD block.
  subroutine readOptionBlock(option, node)

    !> representation of user specified options
    type(TOptionBlock), intent(out) :: option

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

  end subroutine readOptionBlock


  !> Interprets the Data HSD block.
  subroutine readDataBlock(data, node, mode)

    !> representation of dataset information
    type(TDataBlock), intent(out) :: data

    !> node containig the information
    type(fnode), pointer :: node

    !> mode of current run (train, validate, predict)
    character(len=*), intent(in) :: mode

    !> string buffer instances
    type(string) :: strBuffer

    !> true, if files exist
    logical :: tExist

    !> true, if dataset holds corresponding information
    logical :: tStructures, tTargets, tExtFeatures

    call getChildValue(node, 'NetstatFile', strBuffer, default='fortnet.hdf5')

    ! netstat file must be present in validation or prediction mode
    if ((mode == 'validate') .or. (mode == 'predict')) then
      inquire(file=trim(unquote(char(strBuffer))), exist=tExist)
      if (.not. tExist) call detailedError(node, 'Specified netstat file is not present.')
    end if
    data%netstatpath = trim(unquote(char(strBuffer)))

    call getChildValue(node, 'Dataset', strBuffer)
    inquire(file=trim(unquote(char(strBuffer))), exist=tExist)
    if (tExist) then
      data%trainpath = trim(unquote(char(strBuffer)))
    else
      call detailedError(node, 'Specified dataset file is not present.')
    end if

    ! dataset must contain targets in training or validation mode
    if ((mode == 'train') .or. (mode == 'validate')) then
      call inquireTargets(data%trainpath, tTargets)
      if (.not. tTargets) then
        call error('Selected running mode requires the dataset to hold target information.')
      end if
    end if

    ! dataset must always hold some feature information
    call inquireStructures(data%trainpath, tStructures)
    call inquireExtFeatures(data%trainpath, tExtFeatures)
    if ((.not. tStructures) .and. (.not. tExtFeatures)) then
      call error('Neither mappings nor features provided by the dataset.')
    end if

    ! validation dataset only relevant for training runs
    if (mode == 'train') then

      call getChildValue(node, 'Validset', strBuffer, default='')

      select case (char(strBuffer))
      case ('')
        data%tMonitorValid = .false.
      case ('none')
        data%tMonitorValid = .false.
      case default
        inquire(file=trim(unquote(char(strBuffer))), exist=tExist)
        if (tExist) then
          data%tMonitorValid = .true.
          data%validpath = trim(unquote(char(strBuffer)))
        else
          call detailedError(node, 'Specified validset file is not present.')
        end if
      end select

    else

      data%tMonitorValid = .false.

    end if

  end subroutine readDataBlock


  !> Interprets the Network HSD block.
  subroutine readNetworkBlock(network, node, case, nFeatures, nTargets)

    !> representation of network architecture information
    type(TNetworkBlock), intent(out) :: network

    !> node containig the information
    type(fnode), pointer :: node

    !> type of neural network
    character(len=*), intent(in) :: case

    !> number of input features and training targets
    integer, intent(in) :: nFeatures, nTargets

    !> list of integers to parse hidden layer configuration
    type(TListInt) :: integerList

    !> string buffer instance to parse type of transfer function
    type(string) :: strBuffer

    select case (tolower(case))

    case ('bpnn')

      network%type = 'bpnn'

      call init(integerList)
      call getChildValue(node, 'Hidden', integerList)
      allocate(network%hidden(len(integerList)))
      call asArray(integerList, network%hidden)
      call destruct(integerList)
      network%nHiddenLayer = size(network%hidden)

      call getChildValue(node, 'Activation', strBuffer)
      network%activation = tolower(trim(unquote(char(strBuffer))))

    case default

      call detailedError(node, 'Invalid network type.')

    end select

    allocate(network%allDims(network%nHiddenLayer + 2))
    network%allDims(1) = nFeatures
    network%allDims(2:size(network%allDims)-1) = network%hidden
    network%allDims(size(network%allDims)) = nTargets
    network%nSubNnParams = getNumberOfParameters(network%allDims)

  end subroutine readNetworkBlock


  !> Interprets the Training HSD block.
  subroutine readTrainingBlock(train, node, case)

    !> representation of training/optimizer information
    type(TTrainingBlock), intent(out) :: train

    !> node containig the information
    type(fnode), intent(in), pointer :: node

    !> node containing requested regularization
    type(fnode), pointer :: regunode, tmp

    !> type of neural network
    character(len=*), intent(in) :: case

    !> string buffer instances
    type(string) :: strBuffer1, strBuffer2

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

    call getChildValue(node, 'Loss', strBuffer1, 'mse')
    train%lossType = tolower(trim(unquote(char(strBuffer1))))

    call getChildValue(node, 'Shuffle', train%tShuffle, .false.)

    ! read optional regularization technique
    call getChild(node, 'Regularization', child=regunode, requested=.false.)

    if (associated(regunode)) then
      call getChildValue(node, 'Regularization', tmp)
      call getNodeName(tmp, strBuffer2)
      select case (tolower(trim(char(strBuffer2))))
      case ('lasso')
        train%regu%type = 'lasso'
        train%regu%alpha = 1.0_dp
        call getChildValue(tmp, 'Strength', train%regu%strength)
      case ('ridge')
        train%regu%type = 'ridge'
        train%regu%alpha = 0.0_dp
        call getChildValue(tmp, 'Strength', train%regu%strength)
      case ('elasticnet')
        train%regu%type = 'elasticnet'
        call getChildValue(tmp, 'Strength', train%regu%strength)
        call getChildValue(tmp, 'Alpha', train%regu%alpha)
        if ((train%regu%alpha < 0.0_dp) .or. (train%regu%alpha > 1.0_dp)) then
          call detailedError(tmp, 'Elastic net regularization alpha out of bounds [0, 1].')
        end if
      case default
        call detailedError(regunode, 'Invalid regularization type.')
      end select
      if (train%regu%strength < 0.0_dp) then
        call detailedError(tmp, 'Regularzation strength must be zero or positive.')
      elseif (train%regu%strength > 0.0_dp) then
        train%tRegularization = .true.
      else
        train%tRegularization = .false.
        train%regu%type = ''
        train%regu%strength = 0.0_dp
        train%regu%alpha = 0.0_dp
      end if
    else
      train%tRegularization = .false.
      train%regu%type = ''
      train%regu%strength = 0.0_dp
      train%regu%alpha = 0.0_dp
    end if

  end subroutine readTrainingBlock


  !> Initializes the training instance.
  subroutine TTraining_init(this, train)

    !> data type containing variables of the Training block
    type(TTraining), intent(out) :: this

    !> representation of training input block
    type(TTrainingBlock), intent(in) :: train

    call this%setLossFunc(train%lossType)
    call this%setReguFunc(train%regu)

  end subroutine TTraining_init


  !> Initializes the gradient optimizer.
  subroutine TTraining_initOptimizer(this, train, initialGuess)

    !> representation of training/optimizer information
    type(TTraining), intent(inout) :: this

    !> representation of training input block
    type(TTrainingBlock), intent(in) :: train

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

    allocate(this%pOptimizer(nStrucs))

    select case (train%iOptimizer)

    case(optimizerTypes%steepDesc)

      allocate(wrapSteepDesc(nStrucs))
      allocate(weights(nValues))
      weights(:) = train%learningRate

      do iStruc = 1, nStrucs
        allocate(wrapSteepDesc(iStruc)%pSteepDesc)
        call init(wrapSteepDesc(iStruc)%pSteepDesc, nValues, train%threshold,&
            & train%maxDisplacement, weights)
        call init(wrapSteepDesc(iStruc)%pSteepDesc, this%pOptimizer(iStruc))
      end do

      deallocate(weights)

    case (optimizerTypes%conjGrad)

      allocate(wrapConjGrad(nStrucs))

      do iStruc = 1, nStrucs
        allocate(wrapConjGrad(iStruc)%pConjGrad)
        call init(wrapConjGrad(iStruc)%pConjGrad, nValues, train%threshold, train%maxDisplacement)
        call init(wrapConjGrad(iStruc)%pConjGrad, this%pOptimizer(iStruc))
      end do

    case (optimizerTypes%lbfgs)

      allocate(wrapLbfgs(nStrucs))

      do iStruc = 1, nStrucs
        allocate(wrapLbfgs(iStruc)%pLbfgs)
        call TLbfgs_init(wrapLbfgs(iStruc)%pLbfgs, nValues, train%threshold, train%minDisplacement,&
            & train%maxDisplacement, train%mem, train%tLinesearch, .false.,&
            & train%maxForQNDisplacement)
        call init(wrapLbfgs(iStruc)%pLbfgs, this%pOptimizer(iStruc))
      end do

    case (optimizerTypes%fire)

      allocate(wrapFire(nStrucs))

      do iStruc = 1, nStrucs
        allocate(wrapFire(iStruc)%pFire)
        call TFire_init(wrapFire(iStruc)%pFire, nValues, train%threshold, train%maxDisplacement)
        call init(wrapFire(iStruc)%pFire, this%pOptimizer(iStruc))
      end do

    end select

    do iStruc = 1, nStrucs
      call reset(this%pOptimizer(iStruc), initialGuess(:, iStruc))
    end do

  end subroutine TTraining_initOptimizer


  !> Sets the user defined loss function type as well as its derivative.
  subroutine TTraining_setLossFunc(this, descriptor)

    !> data type containing variables of the Training block
    class(TTraining), intent(inout) :: this

    !> descriptor of loss function to use
    character(len=*), intent(in) :: descriptor

    select case(trim(descriptor))

      case('rms')
        this%loss => rmsLoss
        this%lossgrad => rmsGradients

      case('mae')
        this%loss => maLoss
        this%lossgrad => maGradients

      case('mape')
        this%loss => mapLoss
        this%lossgrad => mapGradients

      case('mse')
        this%loss => msLoss
        this%lossgrad => msGradients

      case default
        this%loss => msLoss
        this%lossgrad => msGradients

    end select

  end subroutine TTraining_setLossFunc


  !> Sets the user defined regularization type as well as its derivative.
  subroutine TTraining_setReguFunc(this, regu)

    !> data type containing variables of the Training block
    class(TTraining), intent(inout) :: this

    !> data type containing variables of the Regularization block
    type(TRegularizationBlock), intent(in) :: regu

    select case(trim(regu%type))

      case('lasso')
        this%reguLoss => lassoReguLoss

      case('ridge')
        this%reguLoss => ridgeReguLoss

      case('elasticnet')
        this%reguLoss => elasticNetReguLoss

      case default
        this%reguLoss => nullReguLoss

    end select

  end subroutine TTraining_setReguFunc


  !> Reads the Features HSD block.
  subroutine readFeaturesBlock(this, node, nTotExtFeatures)

    !> representation of feature information
    type(TFeaturesBlock), intent(out) :: this

    !> node containing the information
    type(fnode), intent(in), pointer :: node

    !> total number of external features in the dataset
    integer, intent(in) :: nTotExtFeatures

    !> temporary storage for G-function parameters
    character(len=:), allocatable :: type
    real(dp) :: rCut, kappa, rs, eta, lambda, xi

    !> node containing the information
    type(fnode), pointer :: mappingnode, extnode, child, tmp, indChild

    !> multiple nodes of same name
    type(fnodeList), pointer :: children

    !> string buffer instance
    type(string) :: strBuffer1, strBuffer2

    !> functions emerging from an automatic generation scheme
    type(TGFunctions) :: autoFunctions, tmpFunctions

    !> auxiliary variables
    integer :: iChild, iAutoBlock, nAutoBlocks, nRadial, nAngular, atomId

    call getChild(node, 'External', child=extnode, requested=.false.)

    if (associated(extnode)) then
      call getChildValue(node, 'External', tmp)
      call getNodeName(tmp, strBuffer1)
      call getChild(tmp, 'Indices', child=indChild, requested=.false.)
      if (associated(indChild)) then
        select case (tolower(trim(char(strBuffer1))))
        case ('fromdataset')
          call getChildValue(tmp, 'Indices', strBuffer1, child=child, multiple=.true.)
          call convRangeToInt(char(strBuffer1), extnode, this%ext%indices, nTotExtFeatures)
          call setChildValue(child, '', this%ext%indices, replace=.true.)
          this%ext%nExtFeatures = size(this%ext%indices)
          if (this%ext%nExtFeatures > 0) then
            this%tExtFeatures = .true.
          else
            this%tExtFeatures = .false.
          end if
        case default
          call detailedError(extnode, 'Invalid external feature type.')
        end select
      else
        this%tExtFeatures = .false.
        this%ext%nExtFeatures = 0
      end if
    else
      this%tExtFeatures = .false.
      this%ext%nExtFeatures = 0
    end if

    call getChild(node, 'Mapping', child=mappingnode, requested=.false.)

    if (associated(mappingnode)) then
      call getChildValue(node, 'Mapping', tmp)
      call getNodeName(tmp, strBuffer1)

      call getChildValue(tmp, 'Standardization', this%mapping%tStandardize, .false.)

      select case (tolower(trim(char(strBuffer1))))
      case ('acsf')
        call getChildValue(tmp, 'Reduce', this%mapping%tReduce, .false.)
        ! settings for automatic function and parameter generation
        call getChildren(tmp, 'Function', children)
        if (getLength(children) == 0) then
          call detailedError(tmp, 'No ACSF functions specified.')
        end if
        nAutoBlocks = 0
        iAutoBlock = 0
        do iChild = 1, getLength(children)
          call getItem1(children, iChild, child)
          call getChildValue(child, '', tmp)
          call getNodeName(tmp, strBuffer2)
          type = tolower(trim(char(strBuffer2)))
          if (type == 'auto') then
            nAutoBlocks = nAutoBlocks + 1
          end if
        end do
        allocate(this%mapping%functions%func(getLength(children) - nAutoBlocks))
        do iChild = 1, getLength(children)
          call getItem1(children, iChild, child)
          call getChildValue(child, '', tmp)
          call getNodeName(tmp, strBuffer2)
          type = tolower(trim(char(strBuffer2)))
          
          ! get cutoff radius and convert Angstrom to Bohr
          call getChildValue(tmp, 'RCut', rCut)
          if (rCut <= 0.0_dp) then
            call detailedError(tmp, 'Specified cutoff less or equal zero.')
          end if
          rCut = rCut * AA__Bohr

          ! get optional atom identifier index
          call getChildValue(tmp, 'AtomID', atomId, default=0)
          if (atomId < 0) then
            call detailedError(tmp, 'Specified AtomID below zero.')
          end if
          
          select case (type)
          case ('auto')
            ! settings for automatic function and parameter generation
            call getChildValue(tmp, 'NRadial', nRadial)
            if (nRadial < 0) then
              call detailedError(tmp, 'Specified number of radial functions below zero.')
            end if
            call getChildValue(tmp, 'NAngular', nAngular)
            if (nAngular < 0) then
              call detailedError(tmp, 'Specified number of angular functions below zero.')
            end if
            call tmpFunctions%fromAutoScheme(rCut, nRadial, nAngular, atomId=atomId)
            call autoFunctions%append(tmpFunctions)
            deallocate(tmpFunctions%func)
            iAutoBlock = iAutoBlock + 1
          case ('g1')
            call TGFunction_init(this%mapping%functions%func(iChild - iAutoBlock), type, rCut)
          case ('g2')
            call getChildValue(tmp, 'eta', eta)
            call getChildValue(tmp, 'rs', rs)
            call TGFunction_init(this%mapping%functions%func(iChild - iAutoBlock), type, rCut,&
                & eta=eta, rs=rs)
          case ('g3')
            call getChildValue(tmp, 'kappa', kappa)
            call TGFunction_init(this%mapping%functions%func(iChild - iAutoBlock), type, rCut,&
                & kappa=kappa)
          case ('g4', 'g5')
            call getChildValue(tmp, 'xi', xi)
            call getChildValue(tmp, 'eta', eta)
            call getChildValue(tmp, 'lambda', lambda)
            if ((lambda > 1.0_dp) .or. ((lambda < -1.0_dp))) then
              call detailedError(tmp, 'Specified lambda parameter leaves interval [-1, 1].')
            end if
            call TGFunction_init(this%mapping%functions%func(iChild - iAutoBlock), type, rCut,&
                & xi=xi, eta=eta, lambda=lambda)
          case default
            call detailedError(tmp, 'Invalid mapping type.')
          end select
        end do

        call this%mapping%functions%append(autoFunctions)

        ! getLength(children) == 0 was already checked
        this%tMappingFeatures = .true.

      case default
        call detailedError(mappingnode, 'Invalid mapping type')
      end select

    else
      this%tMappingFeatures = .false.
    end if

    if ((.not. this%tMappingFeatures) .and. (.not. this%tExtFeatures)) then
      call detailedError(node, 'Neither mappings nor features provided.')
    end if

  end subroutine readFeaturesBlock


  !> Processes the Features HSD block.
  subroutine processFeaturesBlock(features, atomicNumbers)

    !> representation of feature information
    type(TFeaturesBlock), intent(inout) :: features

    !> atomic numbers of BPNN sub-nn species
    integer, intent(in) :: atomicNumbers(:)

    ! prevent for calculating non-existent combinations of atomic numbers
    if (size(atomicNumbers) == 1) then
      features%mapping%tReduce = .true.
    end if

    if (features%tMappingFeatures) then
      call processAcsfFunctions(features%mapping%functions, atomicNumbers,&
          & features%mapping%tReduce, features%mapping%nRadial, features%mapping%nAngular)
      features%nFeatures = features%ext%nExtFeatures + size(features%mapping%functions%func)
    else
      features%nFeatures = features%ext%nExtFeatures
    end if

  end subroutine processFeaturesBlock


  !> Processes the ACSF mapping functions.
  subroutine processAcsfFunctions(gFunctions, atomicNumbers, tReduce, nRadial, nAngular)

    !> wrapper around multiple G-functions
    type(TGFunctions), intent(inout) :: gFunctions

    !> atomic numbers of BPNN sub-nn species
    integer, intent(in) :: atomicNumbers(:)

    !> true, if species-unresolved scheme is desired
    logical, intent(in) :: tReduce

    !> obtained number of radial and angular functions
    integer, intent(out) :: nRadial, nAngular

    !> wrapper around multiple G-functions
    type(TGFunctions) :: tmpGfunctions

    !> combinations of atomic number pairs with replacement
    integer, allocatable :: comb(:,:)

    !> auxiliary variables
    integer :: iFunc, iAtNum, iComb, nSpecies, nFunctions

    ! count species-unresolved scheme
    nRadial = 0
    nAngular = 0

    do iFunc = 1, size(gFunctions%func)
      if (gFunctions%func(iFunc)%tRadial) then
        nRadial = nRadial + 1
      else
        nAngular = nAngular + 1
      end if
    end do

    if ((.not. tReduce) .and. (size(atomicNumbers) > 1)) then

      nSpecies = size(atomicNumbers)

      nFunctions = nSpecies * (nRadial + nAngular) + (factorial(nSpecies) * nAngular)&
          & / (2 * factorial(nSpecies - 1))

      ! re-count species-resolved scheme
      nRadial = 0
      nAngular = 0

      allocate(tmpGfunctions%func(nFunctions))

      ! get the possible combinations of two atomic number without double counting
      call combinationsWithReplacement(atomicNumbers, size(atomicNumbers), comb)

      do iFunc = 1, size(gFunctions%func)
        if (gFunctions%func(iFunc)%tRadial) then
          do iAtNum = 1, size(atomicNumbers)
            nRadial = nRadial + 1
            tmpGfunctions%func(nRadial+nAngular) = gFunctions%func(iFunc)
            tmpGfunctions%func(nRadial+nAngular)%atomicNumbers(1) = atomicNumbers(iAtNum)
            tmpGfunctions%func(nRadial+nAngular)%atomicNumbers(2) = 0
          end do
        else
          do iComb = 1, size(comb, dim=2)
            nAngular = nAngular + 1
            tmpGfunctions%func(nRadial+nAngular) = gFunctions%func(iFunc)
            tmpGfunctions%func(nRadial+nAngular)%atomicNumbers(:) = comb(:, iComb)
          end do
        end if
      end do

      if (allocated(gFunctions%func)) deallocate(gFunctions%func)
      call move_alloc(tmpGfunctions%func, gFunctions%func)

    end if

  end subroutine processAcsfFunctions


  !> Checks input consistency by evaluating several assertions.
  subroutine TInput_checkInputConsistency(this)

    !> instance containing parsed input
    class(TInput), intent(in), target :: this

    write(stdOut, '(A)', advance='no') 'Checking Input Consistency...'

    @:ASSERT(size(this%network%hidden) > 0)
    @:ASSERT(minval(this%network%hidden) > 0)

    @:ASSERT(size(this%network%allDims) > 2)
    @:ASSERT(minval(this%network%allDims) > 0)

    @:ASSERT(this%network%nHiddenLayer >= 0)
    select case (tolower(this%option%mode))
    case('train')
      @:ASSERT(this%network%nSubNnParams > 0)
    end select
    @:ASSERT(tolower(this%network%type) == 'bpnn')

    associate (descriptor => this%network%activation)
      @:ASSERT(descriptor == 'gaussian' .or. descriptor == 'relu' .or. descriptor == 'sigmoid'&
          & .or. descriptor == 'heaviside' .or. descriptor == 'tanh' .or. descriptor == 'linear')
    end associate

    if (tolower(this%option%mode) == 'train') then
      associate (descriptor => this%training%lossType)
        @:ASSERT(descriptor == 'mse' .or. descriptor == 'rms' .or. descriptor == 'mae'&
            & .or. descriptor == 'mape')
      end associate
    end if

    if (tolower(this%option%mode) == 'train') then
      @:ASSERT(this%training%nTrainIt >= 0)
      @:ASSERT(this%training%nPrintOut >= 1)
      @:ASSERT(this%training%nSaveNet >= 1)
    end if

    @:ASSERT(this%option%seed >= 0)

    associate (descriptor => this%option%mode)
      @:ASSERT(descriptor == 'train' .or. descriptor == 'validate' .or. descriptor == 'predict')
    end associate

    write(stdOut, '(A)') 'passed'

  end subroutine TInput_checkInputConsistency


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
