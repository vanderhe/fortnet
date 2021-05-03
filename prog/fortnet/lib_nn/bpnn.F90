!--------------------------------------------------------------------------------------------------!
!  FORTNET: A Behler-Parrinello-Neural-Network Implementation                                      !
!  Copyright (C) 2020 - 2021  T. W. van der Heide                                                  !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

module fnet_bpnn

  use dftbp_globalenv, only : stdOut
  use dftbp_accuracy, only: dp
  use dftbp_ranlux, only : TRanlux

  use fnet_initprogram, only : TProgramvariables
  use fnet_nestedtypes, only : TRealArray1D, TRealArray2D, TMultiLayerStruc, TMultiLayerStruc_init
  use fnet_nestedtypes, only : TBiasDerivs, TWeightDerivs, TDerivs, TDerivs_init
  use fnet_nestedtypes, only : TIntArray1D, TPredicts, TPredicts_init
  use fnet_network, only : TNetwork, TNetwork_init
  use fnet_loss, only : deviation
  use fnet_optimizers, only : TOptimizer, next

#:if WITH_MPI
  use fnet_mpifx
  use fnet_parallel, only : getStartAndEndIndex
#:endif

  implicit none

  private

  public :: TBpnn_init, TBpnn


  type :: TBpnn

    !> number of sub-nn's
    integer :: nSpecies

    !> total number of bias and weight parameters of each sub-nn
    integer :: nBiases, nWeights

    !> representation of neural networks
    type(TNetwork), allocatable :: nets(:)

    !> dimensions of all layers in the sub-nn's
    integer, allocatable :: dims(:)

  contains

    procedure :: update => TBpnn_update
    procedure :: serializedWeightsAndBiases => TBpnn_serializedWeightsAndBiases
    procedure :: serialWeightsAndBiasesFillup => TBpnn_serialWeightsAndBiasesFillup
    procedure :: resetActivations => TBpnn_resetActivations
    procedure :: collectOutput => TBpnn_collectOutput
    procedure :: iTrain => TBpnn_iTrain
    procedure :: nTrain => TBpnn_nTrain
    procedure :: iPredict => TBpnn_iPredict
    procedure :: predictBatch => TBpnn_predictBatch
    procedure :: fromFile => TBpnn_fromFile
    procedure :: toFile => TBpnn_toFile

  end type TBpnn


contains

  subroutine TBpnn_init(this, dims, nSpecies, rndGen, activation)

    !> representation of a Behler-Parrinello neural network
    type(TBpnn), intent(out) :: this

    !> dimensions of all layers in the sub-nn's
    integer, intent(in) :: dims(:)

    !> number of different species in training data
    integer, intent(in) :: nSpecies

    !> luxury pseudorandom generator instance
    type(TRanlux), intent(inout), optional :: rndGen

    !> type of activation function to use for all layers, except output layer (linear)
    character(len=*), intent(in), optional :: activation

    !> Species identifier
    integer :: iSpecies

    this%dims = dims
    this%nSpecies = nSpecies

    allocate(this%nets(nSpecies))

    do iSpecies = 1, nSpecies
      call TNetwork_init(this%nets(iSpecies), dims, rndGen=rndGen, descriptor=activation)
    end do

    ! assume that every sub-nn has the same architecture
    this%nBiases = this%nets(1)%nBiases
    this%nWeights = this%nets(1)%nWeights

  end subroutine TBpnn_init


  subroutine TBpnn_nTrain(this, prog, tConverged)

    !> representation of a Behler-Parrinello neural network
    class(TBpnn), intent(inout) :: this

    !> representation of program variables
    type(TProgramVariables), intent(inout) :: prog

    !> true, if gradient got below the specified tolerance
    logical, intent(out) :: tConverged

    !> predictions of a single datapoint
    real(dp), allocatable :: iPredict(:,:)

    !> network predictions during the training
    type(TPredicts) :: predicts, resPredicts, validPredicts

    !> total weight and bias gradients of the current system
    type(TDerivs) :: dd, ddRes

    !> temporary weight and bias gradient storage of the current atom
    type(TDerivs) :: ddTmp

    !> total (validation) rms loss, by comparison of predictions and targets
    real(dp) :: loss, validLoss

    !> euclidean norm of total gradient
    real(dp) :: totGradNorm

    !> true, if current process is the lead
    logical :: tLead

    !> auxiliary variables
    logical :: tPrintOut, tSaveNet

    !> auxiliary variables
    integer :: iIter, iSys, iGlobalSp, iLayer, iStart, iEnd

    call TDerivs_init(this%dims, size(this%nets), dd)
    call TDerivs_init(this%dims, size(this%nets), ddRes)

    call TPredicts_init(predicts, prog%data%targets)
    call TPredicts_init(resPredicts, prog%data%targets)
    call TPredicts_init(validPredicts, prog%data%targets)

  #:if WITH_MPI
    tLead = prog%env%globalMpiComm%lead
    call getStartAndEndIndex(prog%data%nDatapoints, prog%env%globalMpiComm%size,&
        & prog%env%globalMpiComm%rank, iStart, iEnd)
    call syncWeightsAndBiases(this, prog%env%globalMpiComm)
  #:else
    tLead = .true.
    iStart = 1
    iEnd = prog%data%nDatapoints
  #:endif

    lpIter: do iIter = 1, prog%train%nTrainIt

      lpSystem: do iSys = iStart, iEnd
        call this%iTrain(prog%features%features(iSys)%array, prog%data%zTargets(iSys)%array,&
            & prog%data%localAtToGlobalSp(iSys)%array, prog%data%tAtomicTargets, ddTmp, iPredict)
        ! collect outputs and gradients of the systems in range
        predicts%sys(iSys)%array(:,:) = iPredict
        do iGlobalSp = 1, size(this%nets)
          do iLayer = 1, size(this%dims)
            dd%dw(iGlobalSp)%dw(iLayer)%array = dd%dw(iGlobalSp)%dw(iLayer)%array +&
                & ddTmp%dw(iGlobalSp)%dw(iLayer)%array * real(prog%data%weights(iSys), dp)
            dd%db(iGlobalSp)%db(iLayer)%array = dd%db(iGlobalSp)%db(iLayer)%array +&
                & ddTmp%db(iGlobalSp)%db(iLayer)%array * real(prog%data%weights(iSys), dp)
          end do
        end do
      end do lpSystem

    #:if WITH_MPI
      ! collect outputs of all nodes
      do iSys = 1, prog%data%nDatapoints
        call mpifx_allreduce(prog%env%globalMpiComm, predicts%sys(iSys)%array,&
            & resPredicts%sys(iSys)%array, MPI_SUM)
      end do
      ! sync gradients between MPI nodes
      do iGlobalSp = 1, size(this%nets)
        do iLayer = 1, size(this%dims)
          call mpifx_allreduce(prog%env%globalMpiComm, dd%dw(iGlobalSp)%dw(iLayer)%array,&
              & ddRes%dw(iGlobalSp)%dw(iLayer)%array, MPI_SUM)
          call mpifx_allreduce(prog%env%globalMpiComm, dd%db(iGlobalSp)%db(iLayer)%array,&
              & ddRes%db(iGlobalSp)%db(iLayer)%array, MPI_SUM)
        end do
      end do
    #:else
      resPredicts = predicts
      ddRes = dd
    #:endif

      if (tLead) then
        if (prog%data%tZscore) then
          call removeZscore(resPredicts, prog%data%zPrec)
        end if
        loss = prog%train%loss(resPredicts, prog%data%targets)
        if (prog%data%tMonitorValid) then
          validPredicts = this%predictBatch(prog%features%validFeatures,&
              & prog%data%localValidAtToGlobalSp, prog%data%tAtomicTargets, zPrec=prog%data%zPrec)
          validLoss = prog%train%loss(validPredicts, prog%data%validTargets)
        end if
        call this%update(prog%train%pOptimizer, ddRes, loss, sum(prog%data%weights), totGradNorm,&
            & tConverged)
      end if

    #:if WITH_MPI
      call mpifx_bcast(prog%env%globalMpiComm, tConverged)
      call syncWeightsAndBiases(this, prog%env%globalMpiComm)
    #:endif

      tPrintOut = modulo(iIter, prog%train%nPrintOut) == 0 .or. iIter == prog%train%nTrainIt

      if (tPrintOut) then
        if (prog%data%tMonitorValid) then
          write(stdout, '(I10,5X,E15.6,4X,E15.6,4X,E15.6)') iIter, loss, totGradNorm, validLoss
        else
          write(stdout, '(I10,5X,E15.6,4X,E15.6)') iIter, loss, totGradNorm
        end if
      end if

      tSaveNet = modulo(iIter, prog%train%nSaveNet) == 0 .or. iIter == prog%train%nTrainIt

      if (tLead .and. tSaveNet) then
        call this%toFile(prog)
      end if

      if (tConverged) then
        if (prog%data%tMonitorValid) then
          write(stdout, '(I10,5X,E15.6,4X,E15.6,4X,E15.6)') iIter, loss, totGradNorm, validLoss
        else
          write(stdout, '(I10,5X,E15.6,4X,E15.6)') iIter, loss, totGradNorm
        end if
        if (tLead) then
          call this%toFile(prog)
        end if
        exit lpIter
      end if

      call dd%reset()
      call ddRes%reset()

    end do lpIter

  end subroutine TBpnn_nTrain


  subroutine removeZscore(output, zPrec)

    !> neural network output
    type(TPredicts), intent(inout) :: output

    !> storage container of means and variances to calculate z-score, shape: [nTargets, 2]
    real(dp), intent(in) :: zPrec(:,:)

    !> auxiliary variables
    integer :: iSys, iAtom

    do iSys = 1, size(output%sys)
      do iAtom = 1, size(output%sys(iSys)%array, dim=2)
        output%sys(iSys)%array(:, iAtom) = output%sys(iSys)%array(:, iAtom) * zPrec(:, 2) +&
            & zPrec(:, 1)
      end do
    end do

  end subroutine removeZscore


#:if WITH_MPI
  subroutine syncWeightsAndBiases(this, comm)

    !> representation of a Behler-Parrinello neural network
    type(TBpnn), intent(in) :: this

    !> mpi communicator with some additional information
    type(mpifx_comm), intent(in) :: comm

    !> auxiliary variables
    integer :: iGlobalSp, iLayer

    do iGlobalSp = 1, size(this%nets)

      ! broadcast all weights
      do iLayer = 1, size(this%dims)
        call mpifx_bcast(comm, this%nets(iGlobalSp)%layers(iLayer)%ww)
      end do

      ! broadcast all biases
      do iLayer = 1, size(this%dims)
        call mpifx_bcast(comm, this%nets(iGlobalSp)%layers(iLayer)%bb)
      end do

    end do

  end subroutine syncWeightsAndBiases
#:endif


  subroutine TBpnn_iTrain(this, input, targets, localAtToGlobalSp, tAtomic, dd, predicts)

    !> representation of a Behler-Parrinello neural network
    class(TBpnn), intent(inout) :: this

    !> features of all atoms of the current system/geometry, shape: [nFeatures, nAtoms]
    real(dp), intent(in) :: input(:,:)

    !> target data
    !> shape: [nTargets, nAtoms] for atomic targets or [nTargets, 1] for system targets
    real(dp), intent(in) :: targets(:,:)

    !> index mapping local atom --> global species index
    integer, intent(in) :: localAtToGlobalSp(:)

    !> true, if network is trained on atomic properties
    logical, intent(in) :: tAtomic

    !> total weight and bias gradients of the current system
    type(TDerivs), intent(out) :: dd

    !> summed predictions of all sub-nn's
    !> shape: [nTargets, nAtoms] for atomic predictions or [nTargets, 1] for system predictions
    real(dp), intent(out), allocatable :: predicts(:,:)

    !> representation of temporary layer storage container
    type(TMultiLayerStruc) :: tmpLayer

    !> temporary weight gradient storage of the current atom
    type(TWeightDerivs) :: dwTmp

    !> temporary bias gradient storage of the current atom
    type(TBiasDerivs) :: dbTmp

    !> calculated deviation, by comparison of predictions and targets
    real(dp), allocatable :: deviation(:,:)

    !> temporary output storage of sub-nn's
    real(dp), allocatable :: tmpOut(:)

    !> temporary real valued storage
    real(dp), allocatable :: globalPredicts(:)

    !> auxiliary variable
    integer :: iAtom, iGlobalSp, iLayer

    call TMultiLayerStruc_init(this%dims, size(input, dim=2), tmpLayer)
    call TDerivs_init(this%dims, size(this%nets), dd)

    allocate(predicts(size(targets, dim=1), size(input, dim=2)))

    do iAtom = 1, size(input, dim=2)
      iGlobalSp = localAtToGlobalSp(iAtom)
      call this%nets(iGlobalSp)%fprop(input(:, iAtom), state=tmpLayer%struc(iAtom)%layers,&
          & out=tmpOut)
      predicts(:, iAtom) = tmpOut
    end do

    ! calculate absolute deviation between predictions and targets
    if (.not. tAtomic) then
      globalPredicts = sum(predicts, dim=2)
      allocate(deviation(size(predicts, dim=1), size(input, dim=2)))
      do iAtom = 1, size(input, dim=2)
        deviation(:, iAtom) = globalPredicts - targets(:, 1)
      end do
      deallocate(predicts)
      allocate(predicts(size(targets, dim=1), 1))
      predicts(:, 1) = globalPredicts
    else
      deviation = predicts - targets
    end if

    ! collect gradients and normalize them by atom number
    do iAtom = 1, size(input, dim=2)
      iGlobalSp = localAtToGlobalSp(iAtom)
      this%nets(iGlobalSp)%layers = tmpLayer%struc(iAtom)%layers
      call this%nets(iGlobalSp)%bprop(deviation(:, iAtom), dwTmp, dbTmp)
      do iLayer = 1, size(this%dims)
        dd%dw(iGlobalSp)%dw(iLayer)%array = dd%dw(iGlobalSp)%dw(iLayer)%array +&
            & dwTmp%dw(iLayer)%array / real(size(input, dim=2), dp)
        dd%db(iGlobalSp)%db(iLayer)%array = dd%db(iGlobalSp)%db(iLayer)%array +&
            & dbTmp%db(iLayer)%array / real(size(input, dim=2), dp)
      end do
    end do

  end subroutine TBpnn_iTrain


  subroutine TBpnn_update(this, pOptimizer, dd, loss, nDatapoints, totGradNorm, tConverged)

    !> representation of a Behler-Parrinello neural network
    class(TBpnn), intent(inout) :: this

    !> general function optimizer
    type(TOptimizer), intent(inout) :: pOptimizer(:)

    !> weight and bias gradients
    type(TDerivs), intent(in) :: dd

    !> loss function value of current training iteration
    real(dp), intent(in) :: loss

    !> number of systems, used for gradient normalization
    integer, intent(in) :: nDatapoints

    !> euclidean norm of total gradient
    real(dp), intent(out) :: totGradNorm

    !> true, if gradient got below the specified tolerance
    logical, intent(out) :: tConverged

    !> serialized weights and biases
    real(dp), allocatable :: weightsAndBiases(:,:), newWeightsAndBiases(:,:)

    !> serialized weight and bias gradients
    real(dp), allocatable :: ddSerial(:,:)

    !> auxiliary variable
    integer :: iGlobalSp

    call dd%serialized(this%nBiases + this%nWeights, ddSerial)
    call this%serializedWeightsAndBiases(weightsAndBiases)

    ddSerial = ddSerial / real(nDatapoints, dp)
    totGradNorm = norm2(ddSerial)

    allocate(newWeightsAndBiases(size(weightsAndBiases, dim=1), size(weightsAndBiases, dim=2)))

    do iGlobalSp = 1, size(this%nets)
      call next(pOptimizer(iGlobalSp), loss, ddSerial(:, iGlobalSp),&
          & newWeightsAndBiases(:, iGlobalSp), tConverged)
    end do

    call this%serialWeightsAndBiasesFillup(newWeightsAndBiases)

  end subroutine TBpnn_update


  subroutine TBpnn_serialWeightsAndBiasesFillup(this, weightsAndBiases)

    !> representation of a Behler-Parrinello neural network
    class(TBpnn), intent(inout) :: this

    !> serialized weights and biases, shape: [nTotParams, nSpecies]
    real(dp), intent(in) :: weightsAndBiases(:,:)

    !> Species identifier
    integer :: iGlobalSp

    do iGlobalSp = 1, this%nSpecies
      call this%nets(iGlobalSp)%serialWeightsAndBiasesFillup(weightsAndBiases(:, iGlobalSp))
    end do

  end subroutine TBpnn_serialWeightsAndBiasesFillup


  subroutine TBpnn_serializedWeightsAndBiases(this, weightsAndBiases)

    !> representation of a Behler-Parrinello neural network
    class(TBpnn), intent(in) :: this

    !> serialized weights and biases
    real(dp), intent(out), allocatable :: weightsAndBiases(:,:)

    !> temporary buffer
    real(dp), allocatable :: buffer(:)

    !> Species identifier
    integer :: iGlobalSp

    allocate(weightsAndBiases(this%nBiases + this%nWeights, this%nSpecies))

    do iGlobalSp = 1, this%nSpecies
      call this%nets(iGlobalSp)%serializedWeightsAndBiases(buffer)
      weightsAndBiases(:, iGlobalSp) = buffer
    end do

  end subroutine TBpnn_serializedWeightsAndBiases


  subroutine TBpnn_resetActivations(this)

    !> representation of a Behler-Parrinello neural network
    class(TBpnn), intent(inout) :: this

    !> Species identifier
    integer :: iSpecies

    do iSpecies = 1, this%nSpecies
      call this%nets(iSpecies)%resetActivations()
    end do

  end subroutine TBpnn_resetActivations


  subroutine TBpnn_collectOutput(this, output)

    !> representation of a Behler-Parrinello neural network
    class(TBpnn), intent(in) :: this

    !> summed output of all sub-nn's
    real(dp), intent(out), allocatable :: output(:)

    !> temporary store for output parts
    real(dp), allocatable :: tmp(:)

    !> Species identifier
    integer :: iSpecies

    allocate(output(this%dims(size(this%dims))))
    output(:) = 0.0_dp

    do iSpecies = 1, this%nSpecies
      call this%nets(iSpecies)%getOutput(tmp)
      output = output + tmp
    end do

  end subroutine TBpnn_collectOutput


  function TBpnn_iPredict(this, input, localAtToGlobalSp, tAtomic) result(output)

    !> representation of a Behler-Parrinello neural network
    class(TBpnn), intent(in) :: this

    !> single sample of input data to calculate output for
    real(dp), intent(in) :: input(:,:)

    !> index mapping local atom --> global species index
    integer, intent(in) :: localAtToGlobalSp(:)

    !> true, if network is trained on atomic properties
    logical, intent(in) :: tAtomic

    !> atomic contributions to the network prediction
    real(dp), allocatable :: atomicOut(:,:)

    !> final prediction for current input sample
    real(dp), allocatable :: output(:,:)

    !> number of predictions for each sub-nn
    integer :: nOutput

    !> number of atoms in the current system
    integer :: nAtoms

    !> auxiliary variable
    integer :: iAtom, iGlobalSp

    nAtoms = size(input, dim=2)
    nOutput = this%dims(size(this%dims))

    allocate(atomicOut(nOutput, nAtoms))

    do iAtom = 1, nAtoms
      iGlobalSp = localAtToGlobalSp(iAtom)
      atomicOut(:, iAtom) = this%nets(iGlobalSp)%iPredict(input(:, iAtom))
    end do

    if (.not. tAtomic) then
      allocate(output(nOutput, 1))
      output(:, 1) = sum(atomicOut, dim=2)
    else
      output = atomicOut
    end if

  end function TBpnn_iPredict


  function TBpnn_predictBatch(this, features, localAtToGlobalSp, tAtomic, zPrec) result(predicts)

    !> representation of a Behler-Parrinello neural network
    class(TBpnn), intent(in) :: this

    !> atomic features as network input
    type(TRealArray2D), intent(in) :: features(:)

    !> index mapping local atom --> global species index
    type(TIntArray1D), intent(in) :: localAtToGlobalSp(:)

    !> true, if network is trained on atomic properties
    logical, intent(in) :: tAtomic

    !> means and variances to calculate z-score, shape: [nTargets, 2]
    real(dp), intent(in), optional :: zPrec(:,:)

    !> network predictions during the training
    type(TPredicts) :: predicts

    !> auxiliary variable
    integer :: iSys

    allocate(predicts%sys(size(features)))

    do iSys = 1, size(features)
      predicts%sys(iSys)%array = this%iPredict(features(iSys)%array, localAtToGlobalSp(iSys)%array,&
          & tAtomic)
    end do

    if (present(zPrec)) then
      call removeZscore(predicts, zPrec)
    end if

  end function TBpnn_predictBatch


  subroutine TBpnn_toFile(this, prog)

    !> representation of a Behler-Parrinello neural network
    class(TBpnn), intent(in) :: this

    !> representation of program variables
    type(TProgramVariables), intent(in) :: prog

    !> auxiliary variable
    integer :: iSpecies

    do iSpecies = 1, size(this%nets)
      call this%nets(iSpecies)%toFile(prog, iSpecies)
    end do

  end subroutine TBpnn_toFile


  subroutine TBpnn_fromFile(this, prog)

    !> representation of a Behler-Parrinello neural network
    class(TBpnn), intent(out) :: this

    !> representation of program variables
    type(TProgramVariables), intent(inout) :: prog

    !> auxiliary variable
    integer :: iSpecies

    call TBpnn_init(this, prog%arch%allDims, prog%data%nSpecies, activation=prog%arch%activation)

    do iSpecies = 1, prog%data%nSpecies
      call this%nets(iSpecies)%fromFile(prog, iSpecies)
    end do

  end subroutine TBpnn_fromFile

end module fnet_bpnn
