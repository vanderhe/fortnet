!--------------------------------------------------------------------------------------------------!
!  FORTNET: A Behler-Parrinello-Neural-Network Implementation                                      !
!  Copyright (C) 2020 - 2021  T. W. van der Heide                                                  !
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

  use fnet_initprogram, only : TProgramVariables, TProgramVariables_init, TTraining_initOptimizer,&
      & TNetworkBlock, TDataBlock
  use fnet_nestedtypes, only : TEnv, TEnv_init, TPredicts
  use fnet_features, only : TFeatures_init, TFeatures_collect, TMappingBlock, TFeaturesBlock
  use fnet_acsf, only : TAcsf, TAcsf_init
  use fnet_bpnn, only : TBpnn, TBpnn_init
  use fnet_loss, only : minError, maxError, msLoss, rmsLoss, maLoss
  use fnet_netstat, only : readExtFeaturesConfig, writeBpnnHeader, writeExtFeaturesConfig,&
      & inquireExtFeatures, inquireAcsf
  use fnet_iterout, only : writeIterTrajToFile
  use fnet_fnetout, only : writeFnetout
  use fnet_fnetdata, only : TDataset, checkBpnnDatasetCompatibility, checkAcsfDatasetCompatibility,&
      & checkExtFeaturesDatasetCompatibility

#:if WITH_MPI
  use fnet_mpifx
#:endif

  implicit none


  !> instance containing program variables
  type(TProgramVariables) :: prog

  !> Behler-Parrinello-Neural-Network instance
  type(TBpnn) :: bpnn

  !> representation of ACSF mappings
  type(TAcsf) :: trainAcsf, validAcsf

  !> obtained network predictions
  type(TPredicts) :: predicts

  !> true, if current mpi process is the lead
  logical :: tLead

  !> file name of generic Fortnet fortout.xml file
  character(len=*), parameter :: fnetoutFile = 'fnetout.hdf5'

  !> file name of generic Fortnet iterout.dat file
  character(len=*), parameter :: iteroutFile = 'iterout.dat'

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

  call handleInitialisation(prog, bpnn, trainAcsf)

  call printSubnnDetails(prog%inp%network)
  call printMappingDetails(trainAcsf, prog%inp%features)
  call printExternalFeatureDetails(prog%inp%features)
  call printDatasetDetails(prog%trainDataset, prog%validDataset, prog%inp%data,&
      & prog%inp%network%nSubNnParams)

  if (prog%inp%features%tMappingFeatures) then
    call calculateMappings(prog%inp%data, prog%trainDataset, prog%validDataset, trainAcsf,&
        & validAcsf, prog%env)
  end if

  if (tLead .and. (.not. prog%inp%option%tReadNetStats)) then
    if (prog%inp%features%tMappingFeatures) then
      call trainAcsf%toFile(prog%inp%data%netstatpath)
    end if
    if (prog%inp%features%tExtFeatures) then
      call writeExtFeaturesConfig(prog%inp%data%netstatpath, prog%inp%features%ext)
    end if
  end if

  call runCore(prog, bpnn, trainAcsf, validAcsf, predicts)

  if (tLead .and. (prog%inp%option%mode == 'validate' .or. prog%inp%option%mode == 'predict')) then
    call writeFnetout(fnetoutFile, prog%inp%option%mode, prog%trainDataset%targets, predicts,&
        & prog%trainDataset%tAtomicTargets)
  end if

  call destructGlobalEnv()


contains

  !> Initialises several derived types needed for running the core.
  subroutine handleInitialisation(prog, bpnn, trainAcsf)

    !> representation of program variables
    type(TProgramVariables), intent(inout) :: prog

    !> Behler-Parrinello-Neural-Network instance
    type(TBpnn), intent(out) :: bpnn

    !> representation of ACSF mappings
    type(TAcsf), intent(out) :: trainAcsf

    !> serialized initial weights and biases
    real(dp), allocatable :: weightsAndBiases(:,:)

    !> true, if current process is the lead
    logical :: tLead

    !> true, if an ACSF configuration is found in the netstat file
    logical :: tAcsf

    !> true, if the BPNN is trained on atomic targets
    logical :: tAtomicTargets

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
          call bpnn%fromFile(prog%inp%data%netstatpath, tAtomicTargets)
          call checkBpnnDatasetCompatibility(prog%trainDataset, bpnn%atomicNumbers, tAtomicTargets,&
              & allowSpSubset=.false.)
        end if
      #:if WITH_MPI
        call mpifx_bcast(prog%env%globalMpiComm, prog%inp%features%nFeatures)
        call mpifx_bcast(prog%env%globalMpiComm, prog%inp%features%tMappingFeatures)
        call mpifx_bcast(prog%env%globalMpiComm, prog%inp%features%tExtFeatures)
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
            & prog%trainDataset%atomicNumbers, rndGen=prog%rndGen,&
            & activation=prog%inp%network%activation)
        if (tLead) then
          call writeBpnnHeader(prog%inp%data%netstatpath, bpnn, prog%trainDataset%tAtomicTargets)
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
        call bpnn%fromFile(prog%inp%data%netstatpath, tAtomicTargets)
        call checkBpnnDatasetCompatibility(prog%trainDataset, bpnn%atomicNumbers, tAtomicTargets,&
            & allowSpSubset=.true.)
      end if
    #:if WITH_MPI
      call mpifx_bcast(prog%env%globalMpiComm, prog%inp%features%nFeatures)
      call mpifx_bcast(prog%env%globalMpiComm, prog%inp%features%tMappingFeatures)
      call mpifx_bcast(prog%env%globalMpiComm, prog%inp%features%tExtFeatures)
      call mpifx_bcast(prog%env%globalMpiComm, prog%inp%features%tExtFeatures)
      if (prog%inp%features%tMappingFeatures) then
        call trainAcsf%syncConfig(prog%env%globalMpiComm)
      end if
      if (prog%inp%features%tExtFeatures) then
        call prog%inp%features%ext%syncConfig(prog%env%globalMpiComm)
      end if
      call bpnn%sync(prog%env%globalMpiComm)
    #:endif
    end select

    call bpnn%serializedWeightsAndBiases(weightsAndBiases)
    call TTraining_initOptimizer(prog%train, prog%inp%training, weightsAndBiases)

    write(stdout, '(A,/)') repeat('-', 80)

  end subroutine handleInitialisation


  !> Prints some useful details regarding the sub-networks of a BPNN.
  subroutine printSubnnDetails(network)

    !> representation of network information
    type(TNetworkBlock), intent(in) :: network

    !> auxiliary variable
    integer :: ii

    write(stdout, '(A,/)') 'Sub-NN Details'
    write(stdout, '(A,I0)') 'inputs: ', network%allDims(1)
    write(stdout, '(A)', advance='no') 'hidden layers:'
    do ii = 1, size(network%hidden)
      write(stdout, '(A)', advance='no') ' '
      write(stdout, '(I0)', advance='no') network%hidden(ii)
    end do
    write(stdout, '(/,A,I0)') 'outputs: ', network%allDims(size(network%allDims))
    write(stdout, '(/,2A,/)') 'activation: ', network%activation
    write(stdout, '(A,/)') repeat('-', 80)

  end subroutine printSubnnDetails


  !> Prints some useful details regarding the structural mappings.
  subroutine printMappingDetails(acsf, features)

    !> representation of acsf mapping information
    type(TAcsf), intent(in) :: acsf

    !> representation of feature information
    type(TFeaturesBlock), intent(in) :: features

    !> auxiliary variable
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

    !> auxiliary variable
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


  !> Prints some useful details regarding the dataset(s).
  subroutine printDatasetDetails(trainDataset, validDataset, data, nSubNnParams)

    !> representation of a training and validation dataset
    type(TDataset), intent(in) :: trainDataset, validDataset

    !> variables of the Data input block
    type(TDataBlock), intent(in) :: data

    !> number of paramaters (weights + biases) per sub-nn
    integer, intent(in) :: nSubNnParams

    write(stdout, '(A,/)') 'Dataset Information'
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
  subroutine calculateMappings(data, trainDataset, validDataset, trainAcsf, validAcsf, env)

    !> variables of the Data input block
    type(TDataBlock), intent(in) :: data

    !> representation of a training and validation dataset
    type(TDataset), intent(in) :: trainDataset, validDataset

    !> representation of acsf mapping information
    type(TAcsf), intent(inout) :: trainAcsf

    !> representation of acsf mapping information
    type(TAcsf), intent(inout) :: validAcsf

    !> if compiled with mpi enabled, contains mpi communicator
    type(TEnv), intent(in) :: env

    write(stdOut, '(A)', advance='no') 'Calculating ACSF...'
    if (data%tMonitorValid) then
      validAcsf = trainAcsf
    end if

    call trainAcsf%calculate(trainDataset%geos, env, trainDataset%localAtToAtNum,&
        & extFeaturesInp=trainDataset%extFeatures, weights=trainDataset%weights)

    if (data%tMonitorValid) then
      call validAcsf%calculate(validDataset%geos, env, validDataset%localAtToAtNum,&
          & extFeaturesInp=validDataset%extFeatures, zPrec=trainAcsf%zPrec)
    end if

    write(stdout, '(A)') 'done'

  end subroutine calculateMappings


  !> Runs Fortnet's core routines, depending on the obtained configuration.
  subroutine runCore(prog, bpnn, trainAcsf, validAcsf, predicts)

    !> representation of program variables
    type(TProgramVariables), intent(inout) :: prog

    !> Behler-Parrinello-Neural-Network instance
    type(TBpnn), intent(inout) :: bpnn

    !> representation of ACSF mappings
    type(TAcsf), intent(inout) :: trainAcsf, validAcsf

    !> obtained network predictions
    type(TPredicts), intent(out) :: predicts

    !> true, if gradient got below the specified tolerance
    logical :: tConverged

    !> true, if current process is the lead
    logical :: tLead

    !> iteration of lowest training/validation loss
    integer :: iMin, iValidMin

    !> training/validation loss obtained during training
    real(dp), allocatable :: trainLoss(:), validLoss(:)

    !> training gradients obtained during training
    real(dp), allocatable :: gradients(:)

    !> summed loss of predictions, in comparison to targets
    real(dp) :: mse, rms, mae

    !> min. and max. deviation of predictions, in comparison to targets
    real(dp) :: min, max

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
          & prog%train%loss, prog%train%lossgrad, prog%inp%training%tShuffle,&
          & prog%inp%data%tMonitorValid, tConverged, trainLoss=trainLoss, validLoss=validLoss,&
          & gradients=gradients)

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
      predicts = bpnn%predictBatch(prog%features%trainFeatures, prog%env,&
          & prog%trainDataset%localAtToGlobalSp, prog%trainDataset%tAtomicTargets)
      write(stdOut, '(A,/)') 'done'
      write(stdout, '(A,/)') repeat('-', 80)
    end select

    select case(prog%inp%option%mode)
    case ('validate')
      min = minError(predicts, prog%trainDataset%targets)
      max = maxError(predicts, prog%trainDataset%targets)
      mse = msLoss(predicts, prog%trainDataset%targets)
      rms = rmsLoss(predicts, prog%trainDataset%targets)
      mae = maLoss(predicts, prog%trainDataset%targets)
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
