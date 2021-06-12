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
#:endif

  use dftbp_accuracy, only : dp
  use dftbp_message, only : error
  use dftbp_constants, only : Bohr__AA
  use dftbp_typegeometry, only : TGeometry
  use dftbp_charmanip, only : toupper
  use dftbp_globalenv, only : initGlobalEnv, destructGlobalEnv, stdOut

  use fnet_initprogram, only : TProgramVariables, TProgramVariables_init, TArch, TData, TEnv
  use fnet_initprogram, only : initOptimizer, readAcsfFromFile, TExternal
  use fnet_initprogram, only : TFeatures, TFeatures_init, TFeatures_collect, TOption
  use fnet_nestedtypes, only : TEnv_init, TPredicts
  use fnet_acsf, only : TAcsf, TAcsf_init, TAcsfParams_init
  use fnet_bpnn, only : TBpnn, TBpnn_init
  use fnet_loss, only : minError, maxError, msLoss, rmsLoss, maLoss
  use fnet_fnetout, only : writeFnetout
  use fnet_iterout, only : writeIterTrajToFile

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

  !> file name of generic Fortnet iterout.dat file
  character(len=*), parameter :: iteroutFile = 'iterout.dat'

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
  call printMappingDetails(prog%acsf, prog%data%globalSpNames, prog%ext%atomIdIndex)
  call printExternalFeatureDetails(prog%ext)
  call printDatasetDetails(prog%data, prog%arch%nSubNnParams)

  call calculateMappings(prog%ext, prog%data, prog%acsf, prog%validAcsf, prog%env)

  if (tLead .and. (.not. prog%option%tReadNetStats)) then
    call prog%acsf%toFile(acsfFile, prog%data%globalSpNames)
  end if

  call runCore(prog, predicts)

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
        call readAcsfFromFile(prog%mapping, prog%acsf, acsfFile, prog%data%globalSpNames,&
            & prog%data%nSpecies)
        call bpnn%fromFile(prog)
      else
        call TAcsf_init(prog%acsf, prog%mapping%nRadial, prog%mapping%nAngular, prog%mapping%rCut,&
            & prog%ext%speciesIds, tZscore=prog%mapping%tStandardize)
        call TAcsfParams_init(prog%acsf%param, prog%mapping%nRadial, prog%mapping%nAngular,&
            & prog%mapping%rCut)
        call TBpnn_init(bpnn, prog%arch%allDims, prog%data%nSpecies, rndGen=prog%rndGen,&
            & activation=prog%arch%activation)
      end if
    case ('validate', 'predict')
      call readAcsfFromFile(prog%mapping, prog%acsf, acsfFile, prog%data%globalSpNames,&
          & prog%data%nSpecies)
      call bpnn%fromFile(prog)
    end select

    call bpnn%serializedWeightsAndBiases(weightsAndBiases)
    call initOptimizer(prog%train, weightsAndBiases)

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


  subroutine printMappingDetails(acsf, globalSpNames, atomIdIndex)

    !> representation of acsf mapping informations
    type(TAcsf), intent(in) :: acsf

    !> contains (unique) species of all dataset geometries
    character(len=*), intent(in) :: globalSpNames(:)

    !> dataset index of atom identifier
    integer, intent(in) :: atomIdIndex

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
    if (atomIdIndex > 0) then
      write(stdout, '(A,I0)') 'atom id index: ', atomIdIndex
    else
      write(stdout, '(A)') 'atom id index: /'
    end if
    write(stdout, '(A)') ''
    write(stdout, '(A,L1,/)') 'Standardization: ', acsf%tZscore
    write(stdout, '(A,/)') repeat('-', 80)

  end subroutine printMappingDetails


  subroutine printExternalFeatureDetails(ext)

    !> representation of external informations
    type(TExternal), intent(in) :: ext

    !> auxiliary variable
    integer :: ii

    if (ext%tExtFeatures) then

      write(stdout, '(A,/)') 'External Features'
      write(stdout, '(A,I0)') 'nr. of external features: ', ext%nExtFeatures
      write(stdout, '(A)', advance='no') 'dataset indices:'
      do ii = 1, size(ext%extFeaturesInd)
        write(stdout, '(A)', advance='no') ' '
        write(stdout, '(I0)', advance='no') ext%extFeaturesInd(ii)
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


  subroutine calculateMappings(ext, data, acsf, validAcsf, env)

    !> representation of external information
    type(TExternal), intent(in) :: ext

    !> representation of dataset information
    type(TData), intent(in) :: data

    !> representation of acsf mapping information
    type(TAcsf), intent(inout) :: acsf

    !> representation of acsf mapping information
    type(TAcsf), intent(inout) :: validAcsf

    !> if compiled with mpi enabled, contains mpi communicator
    type(TEnv), intent(in) :: env

    write(stdOut, '(A)', advance='no') 'Calculating ACSF...'
    if (data%tMonitorValid) then
      validAcsf = acsf
    end if
    call acsf%calculate(data%geos, env, prog%data%localAtToGlobalSp, weights=prog%data%weights,&
        & atomIds=ext%atomIds)

    if (data%tMonitorValid) then
      call validAcsf%calculate(data%validGeos, env, prog%data%localValidAtToGlobalSp,&
          & atomIds=ext%validAtomIds, zPrec=acsf%zPrec)
    end if

    write(stdout, '(A)') 'done'

  end subroutine calculateMappings


  subroutine runCore(prog, predicts)

    !> representation of program variables
    type(TProgramVariables), intent(inout) :: prog

    !> obtained network predictions
    type(TPredicts), intent(out) :: predicts

    !> true, if gradient got below the specified tolerance
    logical :: tConverged

    !> true, if current process is the lead
    logical :: tLead

    !> total number of input features
    integer :: nFeatures

    !> iteration of lowest training/validation loss
    integer :: iMin, iValidMin

    !> training/validation loss obtained during training
    real(dp), allocatable :: loss(:), validLoss(:)

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

    call TFeatures_init(prog%features, nFeatures, prog%data, prog%ext, prog%acsf, prog%validAcsf)

    if (tLead) then
      call TFeatures_collect(prog%features, prog%data, prog%ext, prog%acsf, prog%validAcsf)
      if (bpnn%dims(1) /= nFeatures) then
        call error('Mismatch in number of features and BPNN input nodes.')
      end if
    end if

  #:if WITH_MPI
    call prog%features%sync(prog%env%globalMpiComm)
  #:endif

    deallocate(prog%data%geos)
    deallocate(prog%acsf%vals%vals)

    if (prog%ext%tExtFeatures) then
      deallocate(prog%ext%extFeatures)
    end if

    if (prog%data%tMonitorValid) then
      deallocate(prog%data%validGeos)
      deallocate(prog%validAcsf%vals%vals)
      if (prog%ext%tExtFeatures) then
        deallocate(prog%ext%extValidFeatures)
      end if
    end if

    select case(prog%option%mode)
    case ('train')
      write(stdOut, '(A,/)') 'Starting training...'
      if (prog%data%tMonitorValid) then
        write(stdOut, '(A11,6X,A13,6X,A13,6X,A13)') 'iTrain', toupper(prog%train%lossType) //&
            & '-Loss', 'Gradients', toupper(prog%train%lossType) // '-V-Loss'
      else
        write(stdOut, '(A11,6X,A13,6X,A13)') 'iTrain', toupper(prog%train%lossType) //&
            & '-Loss', 'Gradients'
      end if
      write(stdout, '(A)') repeat('-', 68)
      call bpnn%nTrain(prog, tConverged, loss=loss, validLoss=validLoss, gradients=gradients)
      if (prog%option%tWriteIterTraj .and. prog%data%tMonitorValid) then
        call writeIterTrajToFile(iteroutFile, loss=loss, validLoss=validLoss, gradients=gradients)
      elseif (prog%option%tWriteIterTraj .and. (.not. prog%data%tMonitorValid)) then
        call writeIterTrajToFile(iteroutFile, loss=loss, gradients=gradients)
      end if
      write(stdout, '(A,/)') repeat('-', 68)
      if (tConverged) then
        write(stdOut, '(A,/)') 'Training finished (gradients converged)'
      else
        write(stdOut, '(A,/)') 'Training finished (max. Iterations reached)'
      end if
      write(stdout, '(A,/)') repeat('-', 68)
      write(stdOut, '(A,/)') 'Loss Analysis (global min.)'
      iMin = minloc(loss, dim=1)
      write(stdOut, '(A,I0,A,(ES12.6E2))') 'iTrain: ', iMin, ', Loss: ', loss(iMin)
      if (prog%data%tMonitorValid) then
        iValidMin = minloc(validLoss, dim=1)
        write(stdOut, '(A,I0,A,(ES12.6E2),/)') 'iTrain: ', iValidMin, ', V-Loss: ',&
            & validLoss(iValidMin)
      else
        write(stdOut, '(A)') ''
      end if
      write(stdout, '(A,/)') repeat('-', 68)
    case ('validate')
      write(stdOut, '(A)', advance='no') 'Start feeding...'
      predicts = bpnn%predictBatch(prog%features%features, prog%env, prog%data%localAtToGlobalSp,&
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
      predicts = bpnn%predictBatch(prog%features%features, prog%env, prog%data%localAtToGlobalSp,&
          & prog%data%tAtomicTargets, zPrec=prog%data%zPrec)
      write(stdOut, '(A,/)') 'done'
      write(stdout, '(A,/)') repeat('-', 80)
    end select

  end subroutine runCore

end program fortnet
