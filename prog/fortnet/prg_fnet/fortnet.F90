!--------------------------------------------------------------------------------------------------!
!  FORTNET: A Behler-Parrinello-Neural-Network Implementation                                      !
!  Copyright (C) 2020 - 2021  T. W. van der Heide                                                  !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

program fortnet

#:if WITH_MPI
  use fnet_mpifx
  use fnet_initprogram, only : syncFeatures
#:endif

  use dftbp_accuracy, only : dp
  use dftbp_message, only : error
  use dftbp_constants, only : Bohr__AA
  use dftbp_typegeometry, only : TGeometry
  use dftbp_charmanip, only : toupper
  use dftbp_globalenv, only : initGlobalEnv, destructGlobalEnv, stdOut

  use fnet_initprogram, only : TProgramVariables, TProgramVariables_init, TArch, TData, TEnv
  use fnet_initprogram, only : initOptimizer, readAcsfFromFile
  use fnet_initprogram, only : TFeatures, TFeatures_init, TFeatures_collect
  use fnet_nestedtypes, only : TEnv_init, TPredicts
  use fnet_acsf, only : TAcsf, TAcsf_init, TAcsfParams_init
  use fnet_bpnn, only : TBpnn, TBpnn_init
  use fnet_loss, only : minError, maxError, msLoss, rmsLoss, maLoss
  use fnet_fnetout, only : writeFnetout

  implicit none


  !> instance containing program variables
  type(TProgramVariables) :: prog

  !> Behler-Parrinello-Neural-Network instance
  type(TBpnn) :: bpnn

  !> obtained network predictions
  type(TPredicts) :: predicts

  !> true, if current mpi process is the lead
  logical :: tLead

  !> acsf file name
  character(len=*), parameter :: acsfFile = 'acsf.out'

  !> file name of generic Fortnet fortout.xml file
  character(len=*), parameter :: fnetoutFile = 'fnetout.xml'

  !> initialise global environment
  call initGlobalEnv()

  !> initialise program variables
  call TProgramVariables_init(prog)
  call TEnv_init(prog%env)

#:if WITH_MPI
  tLead = prog%env%globalMpiComm%lead
#:else
  tLead = .true.
#:endif

  call handleInitialisation(prog, bpnn)

  call printSubnnDetails(prog%arch)
  call printMappingDetails(prog%acsf, prog%data%globalSpNames)
  call printExternalFeatureDetails(prog%data)
  call printDatasetDetails(prog%data, prog%arch%nSubNnParams)

  call calculateMappings(prog%acsf, prog%validAcsf, prog%data, prog%env)

  if (tLead .and. (.not. prog%option%tReadNetStats)) then
    call prog%acsf%toFile(acsfFile, prog%data%globalSpNames)
  end if

  call runCore(prog%option%mode, prog%data, bpnn, prog%features, prog%acsf, prog%validAcsf,&
      & prog%train%lossType, prog%env, predicts)

  if (tLead .and. (prog%option%mode == 'validate' .or. prog%option%mode == 'predict')) then
    call writeFnetout(fnetoutFile, prog%option%mode, prog%data%targets, predicts,&
        & prog%data%tAtomicTargets)
  end if

  call destructGlobalEnv()


contains

  subroutine handleInitialisation(prog, bpnn)

    !> representation of program variables
    type(TProgramVariables), intent(inout) :: prog

    !> Behler-Parrinello-Neural-Network instance
    type(TBpnn), intent(out) :: bpnn

    !> serialized initial weights and biases
    real(dp), allocatable :: weightsAndBiases(:,:)

    write(stdOut, '(A,/)') 'Initialisation'

    select case(prog%option%mode)
    case ('train')
      write(stdOut, '(A)') 'running in training mode'
    case ('validate')
      write(stdOut, '(A)') 'running in validation mode'
    case ('predict')
      write(stdOut, '(A)') 'running in prediction mode'
    end select

    write(stdout, '(A,I0)') 'random seed: ', prog%option%seed
    write(stdOut, '(A,L1,/)') 'read initial netstats: ', prog%option%tReadNetStats

    select case(prog%option%mode)
    case ('train')
      if (prog%option%tReadNetStats) then
        call readAcsfFromFile(prog, acsfFile)
        call bpnn%fromFile(prog)
      else
        call TAcsf_init(prog%acsf, prog%mapping%nRadial, prog%mapping%nAngular, prog%mapping%rCut,&
            & prog%mapping%speciesIds, tZscore=prog%mapping%tStandardize)
        call TAcsfParams_init(prog%acsf%param, prog%mapping%nRadial, prog%mapping%nAngular,&
            & prog%mapping%rCut)
        call TBpnn_init(bpnn, prog%arch%allDims, prog%data%nSpecies, rndGen=prog%rndGen,&
            & activation=prog%arch%activation)
      end if
    case ('validate', 'predict')
      call readAcsfFromFile(prog, acsfFile)
      call bpnn%fromFile(prog)
    end select

    call bpnn%serializedWeightsAndBiases(weightsAndBiases)
    call initOptimizer(prog, weightsAndBiases)

    write(stdout, '(A,/)') repeat('-', 80)

  end subroutine handleInitialisation


  subroutine printSubnnDetails(arch)

    !> representation of architecture informations
    type(TArch), intent(in) :: arch

    !> auxiliary variable
    integer :: ii

    write(stdout, '(A,/)') 'Sub-NN Details'
    write(stdout, '(A,I0)') 'inputs: ', arch%allDims(1)
    write(stdout, '(A)', advance='no') 'hidden layers:'
    do ii = 1, size(arch%hidden)
      write(stdout, '(A)', advance='no') ' '
      write(stdout, '(I0)', advance='no') arch%hidden(ii)
    end do
    write(stdout, '(/,A,I0)') 'outputs: ', arch%allDims(size(arch%allDims))
    write(stdout, '(/,2A,/)') 'activation: ', arch%activation
    write(stdout, '(A,/)') repeat('-', 80)

  end subroutine printSubnnDetails


  subroutine printMappingDetails(acsf, globalSpNames)

    !> representation of acsf mapping informations
    type(TAcsf), intent(in) :: acsf

    !> contains (unique) species of all dataset geometries
    character(len=*), intent(in) :: globalSpNames(:)

    !> auxiliary variable
    integer :: ii

    write(stdout, '(A,/)') 'ACSF Mappings'
    write(stdout, '(A,1F0.4,A)') 'cutoff: ', acsf%rCut * Bohr__AA, ' Angstrom'
    write(stdout, '(A,I0)') 'nr. of radial functions: ', acsf%nRadial
    write(stdout, '(A,I0,/)') 'nr. of angular functions: ', acsf%nAngular
    write(stdout, '(A)') 'species identifier: '
    do ii = 1, size(acsf%speciesIds)
      write(stdout, '(2A,F0.6)') trim(globalSpNames(ii)), ': ', acsf%speciesIds(ii)
    end do
    write(stdout, '(A)') ''
    write(stdout, '(A,L1,/)') 'Standardization: ', acsf%tZscore
    write(stdout, '(A,/)') repeat('-', 80)

  end subroutine printMappingDetails


  subroutine printExternalFeatureDetails(data)

    !> representation of dataset informations
    type(TData), intent(in) :: data

    !> auxiliary variable
    integer :: ii

    if (data%tExtFeatures) then

      write(stdout, '(A,/)') 'External Features'
      write(stdout, '(A,I0)') 'nr. of external features: ', data%nExtFeatures
      write(stdout, '(A)', advance='no') 'dataset indices:'
      do ii = 1, size(data%extFeaturesInd)
        write(stdout, '(A)', advance='no') ' '
        write(stdout, '(I0)', advance='no') data%extFeaturesInd(ii)
      end do
      write(stdout, '(A,/)') ''
      write(stdout, '(A,/)') repeat('-', 80)

    end if

  end subroutine printExternalFeatureDetails


  subroutine printDatasetDetails(data, nSubNnParams)

    !> representation of dataset informations
    type(TData), intent(in) :: data

    !> number of paramaters (weights + biases) per sub-nn
    integer, intent(in) :: nSubNnParams

    write(stdout, '(A,/)') 'Dataset Information'
    write(stdout, '(A,I0,A,I0,A)') 'found: ', sum(data%weights), ' geometries (',&
        & data%nDatapoints, ' unique ones)'
    write(stdout, '(2A)') 'in pathfile: ', data%datapath
    write(stdout, '(A,I0)') 'total sub-nn parameters: ', nSubNnParams
    write(stdout, '(A,F0.4,/)') 'targets per parameter: ', data%nTargetsPerParam
    if (data%tMonitorValid) then
      write(stdout, '(A,I0,3A)') 'found: ', data%nValidDatapoints, ' external geometries'
      write(stdout, '(2A,/)') 'in pathfile: ', data%validpath
    end if
    write(stdout, '(A,/)') repeat('-', 80)

  end subroutine printDatasetDetails


  subroutine calculateMappings(acsf, validAcsf, data, env)

    !> representation of acsf mapping informations
    type(TAcsf), intent(inout) :: acsf

    !> representation of acsf mapping informations
    type(TAcsf), intent(inout) :: validAcsf

    !> representation of dataset informations
    type(TData), intent(inout) :: data

    !> if compiled with mpi enabled, contains mpi communicator
    type(TEnv), intent(in) :: env

    write(stdOut, '(A)', advance='no') 'Calculating ACSF...'
    if (data%tMonitorValid) then
      validAcsf = acsf
    end if
    call acsf%calculate(data%geos, env, prog%data%localAtToGlobalSp, weights=prog%data%weights)

    if (data%tMonitorValid) then
      call validAcsf%calculate(data%validGeos, env, prog%data%localValidAtToGlobalSp,&
          & zPrec=acsf%zPrec)
    end if

    write(stdout, '(A)') 'done'

  end subroutine calculateMappings


  subroutine runCore(mode, data, bpnn, features, acsf, validAcsf, lossType, env, predicts)

    !> mode of current run (train, validate, run)
    character(len=*), intent(in) :: mode

    !> representation of dataset informations
    type(TData), intent(inout) :: data

    !> Behler-Parrinello-Neural-Network instance
    type(TBpnn), intent(inout) :: bpnn

    !> collected features of data and mapping block
    type(TFeatures), intent(inout) :: features

    !> representation of acsf mapping informations
    type(TAcsf), intent(inout) :: acsf

    !> representation of acsf mapping informations
    type(TAcsf), intent(inout) :: validAcsf

    !> type of loss function to use
    character(len=*), intent(in) :: lossType

    !> if compiled with mpi enabled, contains mpi communicator
    type(TEnv), intent(in) :: env

    !> obtained network predictions
    type(TPredicts), intent(out) :: predicts

    !> true, if gradient got below the specified tolerance
    logical :: tConverged

    !> true, if current process is the lead
    logical :: tLead

    !> total number of input features
    integer :: nFeatures

    !> summed loss of predictions, in comparison to targets
    real(dp) :: mse, rms, mae

    !> min. and max. deviation of predictions, in comparison to targets
    real(dp) :: min, max

  #:if WITH_MPI
    tLead = prog%env%globalMpiComm%lead
  #:else
    tLead = .true.
  #:endif

    call TFeatures_init(features, nFeatures, data, acsf, validAcsf)

    if (tLead) then
      call TFeatures_collect(features, data, acsf, validAcsf)
      if (bpnn%dims(1) /= nFeatures) then
        call error('Mismatch in number of features and BPNN input nodes.')
      end if
    end if

  #:if WITH_MPI
    call syncFeatures(features, env%globalMpiComm)
  #:endif

    deallocate(data%geos)
    deallocate(acsf%vals%vals)

    if (data%tExtFeatures) then
      deallocate(data%extFeatures)
    end if

    if (data%tMonitorValid) then
      deallocate(data%validGeos)
      deallocate(validAcsf%vals%vals)
      if (data%tExtFeatures) then
        deallocate(data%extValidFeatures)
      end if
    end if

    select case(mode)
    case ('train')
      write(stdOut, '(A,/)') 'Starting training...'
      if (data%tMonitorValid) then
        write(stdOut, '(A11,6X,A13,6X,A13,6X,A13)') 'iTrain', toupper(lossType) // '-Loss',&
            & 'Gradients', toupper(lossType) // '-V-Loss'
      else
        write(stdOut, '(A11,6X,A13,6X,A13)') 'iTrain', toupper(lossType) // '-Loss', 'Gradients'
      end if
      write(stdout, '(A)') repeat('-', 68)
      call bpnn%nTrain(prog, tConverged)
      write(stdout, '(A,/)') repeat('-', 68)
      if (tConverged) then
        write(stdOut, '(A)') 'Training finished (gradients converged)'
      else
        write(stdOut, '(A)') 'Training finished (max. Iterations reached)'
      end if
    case ('validate')
      write(stdOut, '(A)', advance='no') 'Start feeding...'
      predicts = bpnn%predictBatch(prog%features%features, prog%data%localAtToGlobalSp,&
          & prog%data%tAtomicTargets, zPrec=prog%data%zPrec)
      write(stdOut, '(A,/)') 'done'
      write(stdout, '(A,/)') repeat('-', 80)
      min = minError(predicts, prog%data%targets)
      max = maxError(predicts, prog%data%targets)
      mse = msLoss(predicts, prog%data%targets)
      rms = rmsLoss(predicts, prog%data%targets)
      mae = maLoss(predicts, prog%data%targets)
      write(stdout, '(A,/)') 'Validation'
      write(stdout, '(A,F0.6)') 'min: ', min
      write(stdout, '(A,F0.6,/)') 'max: ', max
      write(stdout, '(A,F0.6)') 'mse: ', mse
      write(stdout, '(A,F0.6)') 'rms: ', rms
      write(stdout, '(A,F0.6,/)') 'mae: ', mae
      write(stdout, '(A,/)') repeat('-', 80)
    case ('predict')
      write(stdOut, '(A)', advance='no') 'Start feeding...'
      predicts = bpnn%predictBatch(prog%features%features, prog%data%localAtToGlobalSp,&
          & prog%data%tAtomicTargets, zPrec=prog%data%zPrec)
      write(stdOut, '(A,/)') 'done'
      write(stdout, '(A,/)') repeat('-', 80)
    end select

  end subroutine runCore

end program fortnet
