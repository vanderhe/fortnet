!--------------------------------------------------------------------------------------------------!
!  FORTNET: A Behler-Parrinello-Neural-Network Implementation                                      !
!  Copyright (C) 2020 - 2023  T. W. van der Heide                                                  !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Defines the general behavior of Fortnet.
program fortnet

  use dftbp_accuracy, only : dp
  use dftbp_message, only : error
  use dftbp_charmanip, only : toupper
  use dftbp_globalenv, only : initGlobalEnv, destructGlobalEnv, stdOut
  use dftbp_typegeometry, only : TGeometry
  use dftbp_simplealgebra, only : determinant33

  use fnet_initprogram, only : TProgramVariables, TProgramVariables_init, TTraining_initOptimizer,&
      & TNetworkBlock, TDataBlock, TAnalysisBlock
  use fnet_nestedtypes, only : TEnv, TEnv_init, TPredicts, TPredicts_init, TIntArray1D,&
      & TRealArray2D, TJacobians
  use fnet_forces, only : TForces, forceAnalysis, TGeometriesForFiniteDiff,&
      & TGeometriesForFiniteDiff_init
  use fnet_features, only : TFeatures_init, TFeatures_collect, TMappingBlock, TFeaturesBlock
  use fnet_acsf, only : TAcsf, TAcsf_init
  use fnet_bpnn, only : TBpnn, TBpnn_init
  use fnet_loss, only : minError, maxError, msLoss, rmsLoss, maLoss, TRegularizationBlock
  use fnet_netstat, only : readExtFeaturesConfig, writeBpnnHeader, writeExtFeaturesConfig,&
      & inquireExtFeatures, inquireAcsf
  use fnet_iterout, only : writeIterTrajToFile
  use fnet_fnetout, only : writeFnetout
  use fnet_fnetdata, only : TDataset, checkBpnnDatasetCompatibility, checkAcsfDatasetCompatibility,&
      & checkExtFeaturesDatasetCompatibility

#:if WITH_MPI
  use fnet_mpifx
#:endif

#:if WITH_SOCKETS
  use fnet_ipisocket, only : ipiSocketComm
  use fnet_fnetdata, only : receiveDatasetFromSocket
  use fnet_features, only : TFeatures_initForSocketComm, TFeatures_collectForSocketComm
#:endif

  implicit none


  !> instance containing program variables
  type(TProgramVariables) :: prog

  !> Behler-Parrinello-Neural-Network instance
  type(TBpnn) :: bpnn

  !> representation of ACSF mappings
  type(TAcsf) :: trainAcsf, validAcsf, forcesAcsf

  !> representation of geometries with shifted atomic coordinates
  type(TGeometriesForFiniteDiff) :: forcesGeos

  !> obtained network predictions
  type(TPredicts) :: predicts

  !> obtained atomic forces
  type(TForces) :: forces

  !> true, if current mpi process is the lead
  logical :: tLead

  !> file name of generic Fortnet fortout.xml file
  character(len=*), parameter :: fnetoutFile = 'fnetout.hdf5'

  !> file name of generic Fortnet iterout.dat file
  character(len=*), parameter :: iteroutFile = 'iterout.dat'

#:if WITH_SOCKETS

  !> representation of ACSF mappings and derivatives
  type(TAcsf) :: acsf, acsfPrime

  !> current geometry iteration (for socket communication)
  integer :: iGeoStep

  !> true, if driver should be stopped
  logical :: tStopDriver

  !> obtained global prediction (i.e. total energy contribution)
  real(dp) :: globalPrediction

  !> obtained atomic forces, shape: [nAtom]
  real(dp), allocatable :: atomicPredictions(:)

  !> obtained atomic forces, shape: [3, nAtom]
  real(dp), allocatable :: atomicForces(:,:)

  !> cell stresses
  real(dp) :: stress(3,3)

  !> cell volume
  real(dp) :: cellVol

#:endif

  !> initialise global environment
  call initGlobalEnv()

  !> initialise program variables
  call TEnv_init(prog%env)
  call TProgramVariables_init(prog)

#:if WITH_MPI
  tLead = prog%env%globalMpiComm%lead
#:else
  tLead = .true.
#:endif

  if (prog%inp%driver%tSocketDriven) then

  #:if WITH_SOCKETS
    ! stress calculation not yet supported, pass zeros instead
    stress(:,:) = 1.0_dp

    ! initializes the BPNN and ACSF instance (only call once)
    call handleInitialisationForSocketComm(prog, bpnn, acsf)

    call printSubnnDetails(prog%inp%network, bpnn%nGlobalTargets, bpnn%nAtomicTargets)
    call printMappingDetails(acsf, prog%inp%features)

    commStep: do iGeoStep = 0, prog%inp%driver%nGeoSteps
      if (iGeoStep < prog%inp%driver%nGeoSteps) then
        ! only receive geometry from socket, if there are still geometry iterations left
        call receiveDatasetFromSocket(prog%env, prog%socket, prog%inp%driver%tPeriodic,&
            & prog%inp%driver%atomicNumbers, bpnn%atomicNumbers, prog%trainDataset, tStopDriver)
        if (tStopDriver) then
          exit commStep
        end if
        if (prog%trainDataset%geos(1)%tPeriodic) then
          cellVol = abs(determinant33(prog%trainDataset%geos(1)%latvecs))
        else
          cellVol = 0.0_dp
        end if
        ! check if the dataset provides structural information needed by the ACSF
        call checkAcsfDatasetCompatibility(prog%trainDataset, acsf)
        ! bypass target number test for prediction runs
        call checkBpnnDatasetCompatibility(prog%trainDataset, bpnn%atomicNumbers,&
            & prog%trainDataset%nGlobalTargets, prog%trainDataset%nAtomicTargets,&
            & allowSpSubset=.true.)
        ! calculate ACSF mappings for the current geometry (call everytime the geometry changed)
        call calculateMappingsForSocketComm(prog%env, prog%trainDataset%geos,&
            &  prog%trainDataset%localAtToAtNum, acsf, acsfPrime)
        call predictForSocketComm(prog, bpnn, acsf, acsfPrime, prog%trainDataset%geos(1),&
            & globalPrediction, atomicPredictions, atomicForces)
        call sendEnergyAndForces(prog%env, prog%socket, globalPrediction, atomicForces, stress,&
            & cellVol)
      end if
    end do commStep

    if (prog%inp%driver%tSocketDriven .and. tLead) then
      call prog%socket%shutdown()
    end if
  #:else
    call error('Internal error: code compiled without socket support')
  #:endif

  else

    call handleInitialisation(prog, bpnn, trainAcsf)

    call printSubnnDetails(prog%inp%network, bpnn%nGlobalTargets, bpnn%nAtomicTargets)
    call printMappingDetails(trainAcsf, prog%inp%features)
    call printExternalFeatureDetails(prog%inp%features)
    call printDatasetDetails(prog%trainDataset, prog%validDataset, prog%inp%data,&
        & prog%inp%network%nSubNnParams)
    call printTrainingDetails(prog%inp%option%mode, prog%inp%training%nTrainIt,&
        & prog%inp%training%iOptimizer, prog%inp%training%lossType, prog%inp%training%threshold,&
        & prog%inp%training%tShuffle, prog%inp%training%tRegularization, prog%inp%training%regu)

    if (prog%inp%features%tMappingFeatures) then
      call calculateMappings(prog%inp%data, prog%inp%analysis, prog%trainDataset,&
          & prog%validDataset, prog%env, trainAcsf, validAcsf, forcesAcsf, forcesGeos)
    end if

    if (tLead .and. (.not. prog%inp%option%tReadNetStats)) then
      if (prog%inp%features%tMappingFeatures) then
        call trainAcsf%toFile(prog%inp%data%netstatpath)
      end if
      if (prog%inp%features%tExtFeatures) then
        call writeExtFeaturesConfig(prog%inp%data%netstatpath, prog%inp%features%ext)
      end if
    end if

    call runCore(prog, bpnn, trainAcsf, validAcsf, forcesAcsf, forcesGeos, predicts, forces)

    if (tLead .and. (prog%inp%option%mode == 'validate'&
        & .or. prog%inp%option%mode == 'predict')) then
      call writeFnetout(fnetoutFile, prog%inp%option%mode, prog%trainDataset%globalTargets,&
          & prog%trainDataset%atomicTargets, predicts, forces, prog%inp%analysis%tForces)
    end if

  end if

  call destructGlobalEnv()


contains

#:if WITH_SOCKETS
  subroutine sendEnergyAndForces(env, socket, energy, forces, stress, cellVol)

    !> Environment
    type(TEnv), intent(in) :: env

    !> Socket may be unallocated (as on follower processes)
    type(ipiSocketComm), intent(inout), allocatable :: socket

    !> total energy
    real(dp), intent(in) :: energy

    !> energy derivatives (-gradients)
    real(dp), intent(in) :: forces(:,:)

    !> stress tensor
    real(dp), intent(in) :: stress(:,:)

    !> cell volume
    real(dp), intent(in) :: cellVol

    !! true, if current process is the lead
    logical :: tLead

  #:if WITH_MPI
    tLead = prog%env%globalMpiComm%lead
  #:else
    tLead = .true.
  #:endif

    if (tLead) then
      ! stress was computed above in the force evaluation block or is 0 if aperiodic
      call socket%send(energy, forces, stress * cellVol)
    end if

  end subroutine sendEnergyAndForces
#:endif


  !> Initialises several derived types needed for running the core.
  subroutine handleInitialisation(prog, bpnn, trainAcsf)

    !> representation of program variables
    type(TProgramVariables), intent(inout) :: prog

    !> Behler-Parrinello-Neural-Network instance
    type(TBpnn), intent(out) :: bpnn

    !> representation of ACSF mappings
    type(TAcsf), intent(out) :: trainAcsf

    !! serialized initial weights and biases
    real(dp), allocatable :: weightsAndBiases(:,:)

    !! true, if current process is the lead
    logical :: tLead

    !! true, if an ACSF configuration is found in the netstat file
    logical :: tAcsf

  #:if WITH_MPI
    tLead = prog%env%globalMpiComm%lead
  #:else
    tLead = .true.
  #:endif

    write(stdOut, '(A,/)') 'Initialisation'

    select case(prog%inp%option%mode)
    case ('train')
      write(stdOut, '(A)') 'running in training mode'
    case ('validate')
      write(stdOut, '(A)') 'running in validation mode'
    case ('predict')
      write(stdOut, '(A)') 'running in prediction mode'
    end select

    write(stdout, '(A,I0)') 'random seed: ', prog%inp%option%seed
    write(stdOut, '(A,L1,/)') 'read initial netstats: ', prog%inp%option%tReadNetStats

    select case(prog%inp%option%mode)
    case ('train')
      if (prog%inp%option%tReadNetStats) then
        if (tLead) then
          prog%inp%features%nFeatures = 0
          call inquireExtFeatures(prog%inp%data%netstatpath, prog%inp%features%tExtFeatures)
          if (prog%inp%features%tExtFeatures) then
            call readExtFeaturesConfig(prog%inp%data%netstatpath, prog%inp%features%ext)
            prog%inp%features%nFeatures = prog%inp%features%nFeatures&
                & + prog%inp%features%ext%nExtFeatures
            ! check if the dataset provides suitable external features
            call checkExtFeaturesDatasetCompatibility(prog%trainDataset,&
                & prog%inp%features%ext%indices)
          end if
          call inquireAcsf(prog%inp%data%netstatpath, prog%inp%features%tMappingFeatures)
          if (prog%inp%features%tMappingFeatures) then
            call trainAcsf%fromFile(prog%inp%data%netstatpath,&
                & tReduce=prog%inp%features%mapping%tReduce,&
                & tStandardize=prog%inp%features%mapping%tStandardize,&
                & nRadial=prog%inp%features%mapping%nRadial,&
                & nAngular=prog%inp%features%mapping%nAngular)
            prog%inp%features%nFeatures = prog%inp%features%nFeatures&
                & + size(trainAcsf%gFunctions%func)
            ! check if the dataset provides structural information needed by the ACSF
            call checkAcsfDatasetCompatibility(prog%trainDataset, trainAcsf)
          end if
          call bpnn%fromFile(prog%inp%data%netstatpath)
          call checkBpnnDatasetCompatibility(prog%trainDataset, bpnn%atomicNumbers,&
              & bpnn%nGlobalTargets, bpnn%nAtomicTargets, allowSpSubset=.false.)
        end if
      #:if WITH_MPI
        call mpifx_bcast(prog%env%globalMpiComm, prog%inp%features%nFeatures)
        call mpifx_bcast(prog%env%globalMpiComm, prog%inp%features%tMappingFeatures)
        call mpifx_bcast(prog%env%globalMpiComm, prog%inp%features%tExtFeatures)
        if (prog%inp%features%tMappingFeatures) then
          call trainAcsf%syncConfig(prog%env%globalMpiComm)
        end if
        if (prog%inp%features%tExtFeatures) then
          call prog%inp%features%ext%syncConfig(prog%env%globalMpiComm)
        end if
        call bpnn%sync(prog%env%globalMpiComm)
      #:endif
      else
        if (prog%inp%features%tMappingFeatures) then
          call TAcsf_init(trainAcsf, prog%inp%features%mapping%functions,&
              & tZscore=prog%inp%features%mapping%tStandardize)
        end if
        call TBpnn_init(bpnn, prog%inp%network%allDims, prog%trainDataset%nSpecies,&
            & prog%trainDataset%nGlobalTargets, prog%trainDataset%nAtomicTargets,&
            & prog%trainDataset%atomicNumbers, rndGen=prog%rndGen,&
            & activation=prog%inp%network%activation)
        if (tLead) then
          call writeBpnnHeader(prog%inp%data%netstatpath, bpnn, prog%trainDataset%nGlobalTargets,&
              & prog%trainDataset%nAtomicTargets)
        end if
      end if
    case ('validate', 'predict')
      if (tLead) then
        prog%inp%features%nFeatures = 0
        call inquireExtFeatures(prog%inp%data%netstatpath, prog%inp%features%tExtFeatures)
        if (prog%inp%features%tExtFeatures) then
          call readExtFeaturesConfig(prog%inp%data%netstatpath, prog%inp%features%ext)
          prog%inp%features%nFeatures = prog%inp%features%nFeatures&
              & + prog%inp%features%ext%nExtFeatures
          ! check if the dataset provides suitable external features
          call checkExtFeaturesDatasetCompatibility(prog%trainDataset,&
              & prog%inp%features%ext%indices)
        end if
        call inquireAcsf(prog%inp%data%netstatpath, prog%inp%features%tMappingFeatures)
        if (prog%inp%features%tMappingFeatures) then
          call trainAcsf%fromFile(prog%inp%data%netstatpath,&
              & tReduce=prog%inp%features%mapping%tReduce,&
              & tStandardize=prog%inp%features%mapping%tStandardize,&
              & nRadial=prog%inp%features%mapping%nRadial,&
              & nAngular=prog%inp%features%mapping%nAngular)
          prog%inp%features%nFeatures = prog%inp%features%nFeatures&
              & + size(trainAcsf%gFunctions%func)
          ! check if the dataset provides structural information needed by the ACSF
          call checkAcsfDatasetCompatibility(prog%trainDataset, trainAcsf)
        end if
        call bpnn%fromFile(prog%inp%data%netstatpath)
        if (prog%inp%option%mode == 'validate') then
          call checkBpnnDatasetCompatibility(prog%trainDataset, bpnn%atomicNumbers,&
              & bpnn%nGlobalTargets, bpnn%nAtomicTargets, allowSpSubset=.true.)
        else
          ! bypass target number test for prediction runs
          call checkBpnnDatasetCompatibility(prog%trainDataset, bpnn%atomicNumbers,&
              & prog%trainDataset%nGlobalTargets, prog%trainDataset%nAtomicTargets,&
              & allowSpSubset=.true.)
        end if
      end if
    #:if WITH_MPI
      call mpifx_bcast(prog%env%globalMpiComm, prog%inp%features%nFeatures)
      call mpifx_bcast(prog%env%globalMpiComm, prog%inp%features%tMappingFeatures)
      call mpifx_bcast(prog%env%globalMpiComm, prog%inp%features%tExtFeatures)
      if (prog%inp%features%tMappingFeatures) then
        call trainAcsf%syncConfig(prog%env%globalMpiComm)
      end if
      if (prog%inp%features%tExtFeatures) then
        call prog%inp%features%ext%syncConfig(prog%env%globalMpiComm)
      end if
      call bpnn%sync(prog%env%globalMpiComm)
    #:endif
      ! check the consistency of specified analysis options
      call prog%inp%checkAnalysisConsistency(prog%trainDataset)
    end select

    call bpnn%serializedWeightsAndBiases(weightsAndBiases)
    call TTraining_initOptimizer(prog%train, prog%inp%training, weightsAndBiases)

    write(stdout, '(A,/)') repeat('-', 80)

  end subroutine handleInitialisation


#:if WITH_SOCKETS

  !> Initialises several derived types needed for running the core, when acting as socket client.
  subroutine handleInitialisationForSocketComm(prog, bpnn, acsf)

    !> representation of program variables
    type(TProgramVariables), intent(inout) :: prog

    !> Behler-Parrinello-Neural-Network instance
    type(TBpnn), intent(out) :: bpnn

    !> representation of ACSF mappings
    type(TAcsf), intent(out) :: acsf

    !! true, if current process is the lead
    logical :: tLead

  #:if WITH_MPI
    tLead = prog%env%globalMpiComm%lead
  #:else
    tLead = .true.
  #:endif

    write(stdOut, '(/A,/)') 'Initialisation:'
    write(stdOut, '(A)') 'running in prediction mode'
    write(stdOut, '(A,1X,A)') 'socket communication to host',&
        & trim(prog%inp%driver%socketInput%host)

    if (tLead) then
      prog%inp%features%nFeatures = 0
      call inquireExtFeatures(prog%inp%data%netstatpath, prog%inp%features%tExtFeatures)
      if (prog%inp%features%tExtFeatures) then
        call error('Socket operation does not support additional external input features.')
      end if
      call inquireAcsf(prog%inp%data%netstatpath, prog%inp%features%tMappingFeatures)
      if (prog%inp%features%tMappingFeatures) then
        call acsf%fromFile(prog%inp%data%netstatpath,&
            & tReduce=prog%inp%features%mapping%tReduce,&
            & tStandardize=prog%inp%features%mapping%tStandardize,&
            & nRadial=prog%inp%features%mapping%nRadial,&
            & nAngular=prog%inp%features%mapping%nAngular)
        if (prog%inp%features%mapping%nRadial + prog%inp%features%mapping%nAngular < 1) then
          call error('Parsed netstat file does not contain ACSF inputs. Aborting.')
        end if
        prog%inp%features%nFeatures = prog%inp%features%nFeatures&
            & + size(acsf%gFunctions%func)
      end if
      call bpnn%fromFile(prog%inp%data%netstatpath)
    end if
  #:if WITH_MPI
    call mpifx_bcast(prog%env%globalMpiComm, prog%inp%features%nFeatures)
    call mpifx_bcast(prog%env%globalMpiComm, prog%inp%features%tMappingFeatures)
    if (prog%inp%features%tMappingFeatures) then
      call acsf%syncConfig(prog%env%globalMpiComm)
    end if
    call bpnn%sync(prog%env%globalMpiComm)
  #:endif

    write(stdout, '(A,/)') repeat('-', 80)

  end subroutine handleInitialisationForSocketComm


  !> Carries out the calculation of structural mappings, if desired.
  subroutine calculateMappingsForSocketComm(env, geos, localAtToAtNum, acsf, acsfPrime)

    !> if compiled with mpi enabled, contains mpi communicator
    type(TEnv), intent(in) :: env

    !> geometry to calculate ACSF for
    type(TGeometry), intent(in) :: geos(:)

    !> index mapping local atom --> atomic number
    type(TIntArray1D), intent(in) :: localAtToAtNum(:)

    !> representation of acsf mapping information
    type(TAcsf), intent(inout) :: acsf, acsfPrime

    write(stdOut, '(A)', advance='no') 'Calculating ACSF...'

    ! copy ACSF configurations in advance of any calculation
    acsfPrime = acsf

    call acsf%calculate(geos, env, localAtToAtNum)
    call acsfPrime%calculatePrime(geos, env, localAtToAtNum)

    write(stdout, '(A)') 'done'

  end subroutine calculateMappingsForSocketComm


  !> Runs Fortnet's core routines, depending on the obtained configuration.
  subroutine predictForSocketComm(prog, bpnn, acsf, acsfPrime, geo, globalPrediction,&
      & atomicPredictions, atomicForces)

    !> representation of program variables
    type(TProgramVariables), intent(inout) :: prog

    !> Behler-Parrinello-Neural-Network instance
    type(TBpnn), intent(inout) :: bpnn

    !> representation of acsf mapping information
    type(TAcsf), intent(inout) :: acsf, acsfPrime

    !> geometry to calculate
    type(TGeometry), intent(in) :: geo

    !> obtained global prediction (i.e. total energy contribution)
    real(dp), intent(out) :: globalPrediction

    !> obtained atomic forces, shape: [nAtom]
    real(dp), intent(out), allocatable :: atomicPredictions(:)

    !> obtained atomic forces, shape: [3, nAtom]
    real(dp), intent(out), allocatable :: atomicForces(:,:)

    !! obtained atomic forces
    type(TForces) :: forces

    !! contains Jacobians of a single systems
    type(TJacobians) :: jacobian

    !! index mapping local atom --> global species index
    type(TIntArray1D), allocatable :: localAtToGlobalSp(:)

    !! additional external, atomic features of the dataset
    real(dp), allocatable :: extFeatures(:,:)

    !! true, if current process is the lead
    logical :: tLead

    !! number of datapoints
    integer :: nDatapoints

    !! number of ACSF mappings
    integer :: nAcsf

    !! number of additional external atomic input features
    integer :: nExtFeatures

    !! total number of atomic input features
    integer :: nFeatures

    !! number of unique elements in geometry
    integer :: nElements

  #:if WITH_MPI
    tLead = prog%env%globalMpiComm%lead
  #:else
    tLead = .true.
  #:endif

    ! we only want to calculate a single structure
    nDatapoints = 1

    ! additional external input features not supported in this case
    nExtFeatures = 0

    if (bpnn%nAtomicTargets /= 0 .or. bpnn%nGlobalTargets /= 1) then
      call error('Socket driven operation only working for exactly one global target.')
    end if

    ! count total number of features
    nAcsf = size(acsf%gFunctions%func)
    nFeatures = nAcsf + nExtFeatures

    ! since there is only one geometry, no difference between local/global species indexing
    allocate(localAtToGlobalSp(1))
    localAtToGlobalSp(1)%array = geo%species

    call TFeatures_initForSocketComm(prog%features, acsf, geo)

    if (tLead) then
      call TFeatures_collectForSocketComm(prog%features, acsf)
      if (bpnn%dims(1) /= nFeatures) then
        call error('Mismatch in number of features and BPNN input nodes.')
      end if
    end if

  #:if WITH_MPI
    call prog%features%sync(prog%env%globalMpiComm)
  #:endif

    if (allocated(acsf%vals%vals)) deallocate(acsf%vals%vals)
    if (allocated(acsfPrime%vals%vals)) deallocate(acsfPrime%vals%vals)

    write(stdOut, '(A)', advance='no') 'Feeding networks...'
    atomicPredictions = reshape(bpnn%iPredict(prog%features%trainFeatures(1)%array, geo%species),&
        & [geo%nAtom])
    write(stdOut, '(A)') 'done'
    globalPrediction = sum(atomicPredictions)

    write(stdOut, '(A)', advance='no') 'Calculating Jacobian...'
    allocate(jacobian%sys(nDatapoints))
    jacobian%sys(1) = bpnn%iJacobian(prog%features%trainFeatures(1)%array, geo%species)
    write(stdOut, '(A)') 'done'
    write(stdOut, '(A)', advance='no') 'Calculating forces...'
    forces = forceAnalysis(bpnn, prog%features%trainFeatures, jacobian, acsfPrime, prog%env,&
        & localAtToGlobalSp)
    write(stdOut, '(A,/)') 'done'
    atomicForces = forces%geos(1)%array

    write(stdout, '(A,/)') repeat('-', 80)

  end subroutine predictForSocketComm

#:endif


  !> Prints some useful details regarding the sub-networks of a BPNN.
  subroutine printSubnnDetails(network, nGlobalTargets, nAtomicTargets)

    !> representation of network information
    type(TNetworkBlock), intent(in) :: network

    !> number of system-wide training targets of BPNN
    integer, intent(in) :: nGlobalTargets

    !> number of atomic training targets of BPNN
    integer, intent(in) :: nAtomicTargets

    !! auxiliary variable
    integer :: ii

    write(stdout, '(A,/)') 'Sub-NN Details'
    write(stdout, '(A,I0)') 'inputs: ', network%allDims(1)
    write(stdout, '(A)', advance='no') 'hidden layers:'
    do ii = 1, size(network%hidden)
      write(stdout, '(A)', advance='no') ' '
      write(stdout, '(I0)', advance='no') network%hidden(ii)
    end do
    write(stdout, '(/,A,I0)') 'outputs: ', network%allDims(size(network%allDims))
    write(stdout, '(/,A,I0)') 'nr. of global outputs: ', nGlobalTargets
    write(stdout, '(A,I0,/)') 'nr. of atomic outputs: ', nAtomicTargets
    write(stdout, '(2A,/)') 'activation: ', network%activation
    write(stdout, '(A,/)') repeat('-', 80)

  end subroutine printSubnnDetails


  !> Prints some useful details regarding the structural mappings.
  subroutine printMappingDetails(acsf, features)

    !> representation of acsf mapping information
    type(TAcsf), intent(in) :: acsf

    !> representation of feature information
    type(TFeaturesBlock), intent(in) :: features

    !! auxiliary variable
    integer :: iFunc

    if (features%tMappingFeatures) then

      write(stdout, '(A,/)') 'ACSF Mappings'
      write(stdout, '(A,L1,/)') 'species-resolved: ', (.not. features%mapping%tReduce)
      write(stdout, '(A,I0)') 'nr. of radial functions: ', features%mapping%nRadial
      write(stdout, '(A,I0,/)') 'nr. of angular functions: ', features%mapping%nAngular

      do iFunc = 1, size(acsf%gFunctions%func)
        select case (acsf%gFunctions%func(iFunc)%type)
        case ('g1')
          write(stdout, '(2A,1F0.6,A)') acsf%gFunctions%func(iFunc)%type, ': rc = ',&
              & acsf%gFunctions%func(iFunc)%rCut, ','
          if (features%mapping%tReduce) then
            write(stdout, '(A,I0)') '    atomId = ', acsf%gFunctions%func(iFunc)%atomId
          else
            write(stdout, '(A,I0,A,I0)') '    atomId = ', acsf%gFunctions%func(iFunc)%atomId,&
                & ', atNum = ', acsf%gFunctions%func(iFunc)%atomicNumbers(1)
          end if
        case ('g2')
          write(stdout, '(2A,1F0.6,A,1F0.6,A,1F0.6,A)') acsf%gFunctions%func(iFunc)%type,&
              & ': rc = ', acsf%gFunctions%func(iFunc)%rCut, ', rs = ',&
              & acsf%gFunctions%func(iFunc)%rs, ', eta = ', acsf%gFunctions%func(iFunc)%eta, ','
          if (features%mapping%tReduce) then
            write(stdout, '(A,I0)') '    atomId = ', acsf%gFunctions%func(iFunc)%atomId
          else
            write(stdout, '(A,I0,A,I0)') '    atomId = ', acsf%gFunctions%func(iFunc)%atomId,&
                & ', atNum = ', acsf%gFunctions%func(iFunc)%atomicNumbers(1)
          end if
        case ('g3')
          write(stdout, '(2A,1F0.6,A,1F0.6,A)') acsf%gFunctions%func(iFunc)%type, ': rc = ',&
              & acsf%gFunctions%func(iFunc)%rCut, ', kappa = ', acsf%gFunctions%func(iFunc)%kappa,&
              & ','
          if (features%mapping%tReduce) then
            write(stdout, '(A,I0)') '    atomId = ', acsf%gFunctions%func(iFunc)%atomId
          else
            write(stdout, '(A,I0,A,I0)') '    atomId = ', acsf%gFunctions%func(iFunc)%atomId,&
                & ', atNum = ', acsf%gFunctions%func(iFunc)%atomicNumbers(1)
          end if
        case ('g4', 'g5')
          write(stdout, '(2A,1F0.6,A,1F0.6,A,1F0.6,A,1F0.6)') acsf%gFunctions%func(iFunc)%type,&
              & ': rc = ', acsf%gFunctions%func(iFunc)%rCut, ', lambda = ',&
              & acsf%gFunctions%func(iFunc)%lambda, ', eta = ', acsf%gFunctions%func(iFunc)%eta,&
              & ', xi = ', acsf%gFunctions%func(iFunc)%xi
          if (features%mapping%tReduce) then
            write(stdout, '(A,I0)') '    atomId = ', acsf%gFunctions%func(iFunc)%atomId
          else
            write(stdout, '(A,I0,A,I0,A,I0,A)') '    atomId = ',&
                & acsf%gFunctions%func(iFunc)%atomId, ', atNums = [',&
                & acsf%gFunctions%func(iFunc)%atomicNumbers(1), ', ',&
                & acsf%gFunctions%func(iFunc)%atomicNumbers(2), ']'
          end if
        case default
          call error('Invalid function type, aborting printing.')
        end select
      end do

      write(stdout, '(A,/)') repeat('-', 80)

    end if

  end subroutine printMappingDetails


  !> Prints some useful details regarding external atomic input features.
  subroutine printExternalFeatureDetails(features)

    !> representation of feature information
    type(TFeaturesBlock), intent(in) :: features

    !! auxiliary variable
    integer :: ii

    if (features%tExtFeatures) then

      write(stdout, '(A,/)') 'External Features'
      write(stdout, '(A,I0)') 'nr. of external features: ', features%ext%nExtFeatures
      write(stdout, '(A)', advance='no') 'dataset indices:'
      do ii = 1, size(features%ext%indices)
        write(stdout, '(A)', advance='no') ' '
        write(stdout, '(I0)', advance='no') features%ext%indices(ii)
      end do
      write(stdout, '(A,/)') ''
      write(stdout, '(A,/)') repeat('-', 80)

    end if

  end subroutine printExternalFeatureDetails


  !> Prints some useful details regarding the training process.
  subroutine printTrainingDetails(mode, nTrainIt, iOptimizer, lossType, threshold, tShuffle,&
      & tRegularization, regu)

    !> running mode oif current run
    character(len=*), intent(in) :: mode

    !> maximum number of training iterations
    integer, intent(in) :: nTrainIt

    !> integer ID of specified optimizert
    integer, intent(in) :: iOptimizer

    !> type of loss function to use during the training
    character(len=*), intent(in) :: lossType

    !> gradient threshold where to stop the training, if provided
    real(dp), intent(in) :: threshold

    !> wether a Knuth-shuffle should be applied to the gradient calculation of datapoints
    logical, intent(in) :: tShuffle

    !> true, if loss-based regularization is requested
    logical, intent(in) :: tRegularization

    !> data type containing variables of the Regularization block
    type(TRegularizationBlock), intent(in) :: regu

    !! string associated with the integer optimizer ID
    character(len=:), allocatable :: optimizerStr

    if (mode == 'train') then

      select case(iOptimizer)
      case(1)
        optimizerStr = 'Steepest Descent'
      case(2)
        optimizerStr = 'Conjugate Gradients'
      case(3)
        optimizerStr = 'L-BFGS'
      case(4)
        optimizerStr = 'FIRE'
      case default
        call error('Could not identify gradient optimizer by integer ID.')
      end select

      write(stdout, '(A,/)') 'Training Information'
      write(stdout, '(A,I0)') 'max. iTrain: ', nTrainIt
      write(stdout, '(A,F0.4)') 'converged at: ', threshold
      write(stdout, '(2A)') 'optimizer: ', optimizerStr
      write(stdout, '(2A)') 'loss: ', lossType
      write(stdout, '(A,L1)') 'shuffle gradients: ', tShuffle
      if (tRegularization) then
        write(stdout, '(2A,/)') 'regularization: ', regu%type
      else
        write(stdout, '(A,/)') 'regularization: /'
      end if
      write(stdout, '(A,/)') repeat('-', 80)

    end if

  end subroutine printTrainingDetails


  !> Prints some useful details regarding the dataset(s).
  subroutine printDatasetDetails(trainDataset, validDataset, data, nSubNnParams)

    !> representation of a training and validation dataset
    type(TDataset), intent(in) :: trainDataset, validDataset

    !> variables of the Data input block
    type(TDataBlock), intent(in) :: data

    !> number of paramaters (weights + biases) per sub-nn
    integer, intent(in) :: nSubNnParams

    write(stdout, '(A,/)') 'Dataset Information'
    write(stdout, '(A,I0)') 'nr. of global targets: ', trainDataset%nGlobalTargets
    write(stdout, '(A,I0)') 'nr. of atomic targets: ', trainDataset%nAtomicTargets
    write(stdout, '(A,I0,A,I0,A)') 'found: ', sum(trainDataset%weights), ' datapoints (',&
        & trainDataset%nDatapoints, ' unique ones)'
    write(stdout, '(2A)') 'in file: ', data%trainpath
    write(stdout, '(A,I0)') 'total sub-nn parameters: ', nSubNnParams
    write(stdout, '(A,F0.4,/)') 'targets per parameter: ', trainDataset%nTargetsPerParam
    if (data%tMonitorValid) then
      write(stdout, '(A,I0,3A)') 'found: ', validDataset%nDatapoints, ' external datapoints'
      write(stdout, '(2A,/)') 'in file: ', data%validpath
    end if
    write(stdout, '(A,/)') repeat('-', 80)

  end subroutine printDatasetDetails


  !> Carries out the calculation of structural mappings, if desired.
  subroutine calculateMappings(data, analysis, trainDataset, validDataset, env, trainAcsf,&
      & validAcsf, forcesAcsf, forcesGeos)

    !> variables of the Data input block
    type(TDataBlock), intent(in) :: data

    !> variables of the Analysis input block
    type(TAnalysisBlock), intent(in) :: analysis

    !> representation of a training and validation dataset
    type(TDataset), intent(in) :: trainDataset, validDataset

    !> if compiled with mpi enabled, contains mpi communicator
    type(TEnv), intent(in) :: env

    !> representation of acsf mapping information
    type(TAcsf), intent(inout) :: trainAcsf, validAcsf, forcesAcsf

    !> representation of geometries with shifted atomic coordinates
    type(TGeometriesForFiniteDiff), intent(out) :: forcesGeos

    !! serialized representation of geometries with shifted atomic coordinates
    type(TGeometry), allocatable :: forcesSerializedGeos(:)

    !! index mapping local atom --> atomic number
    type(TIntArray1D), allocatable :: localAtToAtNum(:)

    !! atom dependent scaling parameters for cutoff function
    type(TRealArray2D), allocatable :: extFeatures(:)

    !! auxiliary variables
    integer :: iGeo, ii, nTotGeometries, ind

    write(stdOut, '(A)', advance='no') 'Calculating ACSF...'

    ! copy ACSF configurations in advance of their calculation
    if (data%tMonitorValid) then
      validAcsf = trainAcsf
    end if
    if (analysis%tForces) then
      forcesAcsf = trainAcsf
    end if

    call trainAcsf%calculate(trainDataset%geos, env, trainDataset%localAtToAtNum,&
        & extFeaturesInp=trainDataset%extFeatures, weights=trainDataset%weights)

    if (data%tMonitorValid) then
      call validAcsf%calculate(validDataset%geos, env, validDataset%localAtToAtNum,&
          & extFeaturesInp=validDataset%extFeatures, zPrec=trainAcsf%zPrec)
    end if

    if (analysis%tForces) then
      if (prog%inp%analysis%forceMethod == 'finitedifferences') then
        call TGeometriesForFiniteDiff_init(forcesGeos, trainDataset%geos, analysis%delta)
        nTotGeometries = 0
        do iGeo = 1, size(forcesGeos%geos)
          nTotGeometries = nTotGeometries + 6 * size(forcesGeos%geos(iGeo)%atom)
        end do
        allocate(localAtToAtNum(nTotGeometries))
        if (trainDataset%tExtFeatures) allocate(extFeatures(nTotGeometries))
        ind = 1
        do iGeo = 1, size(forcesGeos%geos)
          do ii = 1, 6 * size(forcesGeos%geos(iGeo)%atom)
            localAtToAtNum(ind+ii-1)%array = trainDataset%localAtToAtNum(iGeo)%array
            if (trainDataset%tExtFeatures) then
              extFeatures(ind+ii-1)%array = trainDataset%extFeatures(iGeo)%array
            end if
          end do
          ind = ind + 6 * size(forcesGeos%geos(iGeo)%atom)
        end do
        call forcesGeos%serialize(forcesSerializedGeos)
        call forcesAcsf%calculate(forcesSerializedGeos, env, localAtToAtNum, zPrec=trainAcsf%zPrec,&
            & extFeaturesInp=extFeatures)
      elseif (prog%inp%analysis%forceMethod == 'analytical') then
        call forcesAcsf%calculatePrime(trainDataset%geos, env, trainDataset%localAtToAtNum,&
            & extFeaturesInp=trainDataset%extFeatures)
      end if
    end if

    write(stdout, '(A)') 'done'

  end subroutine calculateMappings


  !> Runs Fortnet's core routines, depending on the obtained configuration.
  subroutine runCore(prog, bpnn, trainAcsf, validAcsf, forcesAcsf, forcesGeos, predicts, forces)

    !> representation of program variables
    type(TProgramVariables), intent(inout) :: prog

    !> Behler-Parrinello-Neural-Network instance
    type(TBpnn), intent(inout) :: bpnn

    !> representation of ACSF mappings
    type(TAcsf), intent(inout) :: trainAcsf, validAcsf, forcesAcsf

    !> representation of geometries with shifted atomic coordinates
    type(TGeometriesForFiniteDiff), intent(in) :: forcesGeos

    !> obtained network predictions
    type(TPredicts), intent(out) :: predicts

    !> obtained atomic forces
    type(TForces), intent(out) :: forces

    !! true, if gradient got below the specified tolerance
    logical :: tConverged

    !! true, if current process is the lead
    logical :: tLead

    !! iteration of lowest training/validation loss
    integer :: iMin, iValidMin

    !! training/validation loss obtained during training
    real(dp), allocatable :: trainLoss(:), validLoss(:)

    !! training gradients obtained during training
    real(dp), allocatable :: gradients(:)

    !! summed loss of predictions, in comparison to targets
    real(dp) :: mse, rms, mae

    !! min. and max. deviation of predictions, in comparison to targets
    real(dp) :: min, max

    !! contains Jacobians of multiple systems
    type(TJacobians) :: jacobian

  #:if WITH_MPI
    tLead = prog%env%globalMpiComm%lead
  #:else
    tLead = .true.
  #:endif

    call TFeatures_init(prog%features, prog%inp%features, prog%trainDataset, prog%validDataset,&
        & trainAcsf, prog%inp%data%tMonitorValid)

    if (tLead) then
      call TFeatures_collect(prog%features, prog%inp%features, prog%trainDataset,&
          & prog%validDataset, trainAcsf, validAcsf, prog%inp%data%tMonitorValid)
      if (bpnn%dims(1) /= prog%inp%features%nFeatures) then
        call error('Mismatch in number of features and BPNN input nodes.')
      end if
    end if

  #:if WITH_MPI
    call prog%features%sync(prog%env%globalMpiComm)
  #:endif

    if (allocated(prog%trainDataset%geos)) deallocate(prog%trainDataset%geos)
    if (allocated(prog%validDataset%geos)) deallocate(prog%validDataset%geos)

    if (allocated(trainAcsf%vals%vals)) deallocate(trainAcsf%vals%vals)
    if (allocated(validAcsf%vals%vals)) deallocate(validAcsf%vals%vals)

    if (allocated(prog%trainDataset%extFeatures)) deallocate(prog%trainDataset%extFeatures)
    if (allocated(prog%validDataset%extFeatures)) deallocate(prog%validDataset%extFeatures)

    select case(prog%inp%option%mode)
    case ('train')
      write(stdOut, '(A,/)') 'Starting training...'
      if (prog%inp%data%tMonitorValid) then
        write(stdOut, '(A11,6X,A13,6X,A13,6X,A13)') 'iTrain', toupper(prog%inp%training%lossType)&
            & // '-Loss', 'Gradients', toupper(prog%inp%training%lossType) // '-V-Loss'
      else
        write(stdOut, '(A11,6X,A13,6X,A13)') 'iTrain', toupper(prog%inp%training%lossType) //&
            & '-Loss', 'Gradients'
      end if
      write(stdout, '(A)') repeat('-', 68)
      call bpnn%nTrain(prog%env, prog%rndGen, prog%train%pOptimizer, prog%trainDataset,&
          & prog%validDataset, prog%features, prog%inp%training%nTrainIt,&
          & prog%inp%training%nPrintOut, prog%inp%training%nSaveNet, prog%inp%data%netstatpath,&
          & prog%train%loss, prog%train%lossgrad, prog%train%reguLoss, prog%inp%training%regu,&
          & prog%inp%training%tShuffle, prog%inp%data%tMonitorValid, tConverged,&
          & trainLoss=trainLoss, validLoss=validLoss, gradients=gradients)

      if (tLead) then
        if (prog%inp%option%tWriteIterTraj .and. prog%inp%data%tMonitorValid) then
          call writeIterTrajToFile(iteroutFile, trainLoss=trainLoss, validLoss=validLoss,&
              & gradients=gradients)
        elseif (prog%inp%option%tWriteIterTraj .and. (.not. prog%inp%data%tMonitorValid)) then
          call writeIterTrajToFile(iteroutFile, trainLoss=trainLoss, gradients=gradients)
        end if
      end if
      write(stdout, '(A,/)') repeat('-', 68)
      if (tConverged) then
        write(stdOut, '(A,/)') 'Training finished (gradients converged)'
      else
        write(stdOut, '(A,/)') 'Training finished (max. Iterations reached)'
      end if
      write(stdout, '(A,/)') repeat('-', 68)
      write(stdOut, '(A,/)') 'Loss Analysis (global min.)'
      iMin = minloc(trainLoss, dim=1)
      write(stdOut, '(A,I0,A,(ES12.6E2))') 'iTrain: ', iMin, ', Loss: ', trainLoss(iMin)
      if (prog%inp%data%tMonitorValid) then
        iValidMin = minloc(validLoss, dim=1)
        write(stdOut, '(A,I0,A,(ES12.6E2),/)') 'iTrain: ', iValidMin, ', V-Loss: ',&
            & validLoss(iValidMin)
      else
        write(stdOut, '(A)') ''
      end if
      write(stdout, '(A,/)') repeat('-', 68)
    case ('validate', 'predict')
      write(stdOut, '(A)', advance='no') 'Start feeding...'
      call TPredicts_init(predicts, prog%trainDataset%nDatapoints, bpnn%nGlobalTargets,&
          & bpnn%nAtomicTargets, prog%trainDataset%localAtToAtNum)
      predicts%sys = bpnn%predictBatch(prog%features%trainFeatures, prog%env,&
          & prog%trainDataset%localAtToGlobalSp)
      write(stdOut, '(A)') 'done'
      if (prog%inp%analysis%tForces) then
        write(stdOut, '(A)', advance='no') 'Calculate forces...'
        if (prog%inp%analysis%forceMethod == 'finitedifferences') then
          forces = forceAnalysis(bpnn, forcesGeos, forcesAcsf, prog%env,&
              & prog%trainDataset%localAtToGlobalSp, prog%inp%analysis%delta)
        elseif (prog%inp%analysis%forceMethod == 'analytical') then
          jacobian = bpnn%nJacobian(prog%features%trainFeatures, prog%env,&
              & prog%trainDataset%localAtToGlobalSp)
          forces = forceAnalysis(bpnn, prog%features%trainFeatures, jacobian, forcesAcsf, prog%env,&
              & prog%trainDataset%localAtToGlobalSp)
        end if
        write(stdOut, '(A,/)') 'done'
      else
        write(stdOut, '(A)') ''
      end if
      write(stdout, '(A,/)') repeat('-', 80)
    end select

    select case(prog%inp%option%mode)
    case ('validate')
      min = minError(predicts, prog%trainDataset%globalTargets, prog%trainDataset%atomicTargets)
      max = maxError(predicts, prog%trainDataset%globalTargets, prog%trainDataset%atomicTargets)
      mse = msLoss(predicts, prog%trainDataset%globalTargets, prog%trainDataset%atomicTargets,&
          & prog%trainDataset%atomicWeights)
      rms = rmsLoss(predicts, prog%trainDataset%globalTargets, prog%trainDataset%atomicTargets,&
          & prog%trainDataset%atomicWeights)
      mae = maLoss(predicts, prog%trainDataset%globalTargets, prog%trainDataset%atomicTargets,&
          & prog%trainDataset%atomicWeights)
      write(stdout, '(A,/)') 'Validation'
      write(stdout, '(A,E15.6)') 'min. absolute error: ', min
      write(stdout, '(A,E15.6,/)') 'max. absolute error: ', max
      write(stdout, '(A,E15.6)') 'mean average error:      ', mae
      write(stdout, '(A,E15.6)') 'mean squared error:      ', mse
      write(stdout, '(A,E15.6,/)') 'root mean squared error: ', rms
      write(stdout, '(A,/)') repeat('-', 80)
    end select

  end subroutine runCore

end program fortnet
