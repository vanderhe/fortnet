!--------------------------------------------------------------------------------------------------!
!  FORTNET: A Behler-Parrinello-Neural-Network Implementation                                      !
!  Copyright (C) 2020 - 2022  T. W. van der Heide                                                  !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Implements the Behler-Parrinello-Neural-Network topology.
module fnet_bpnn

  use h5lt
  use hdf5

  use dftbp_message, only : error
  use dftbp_globalenv, only : stdOut
  use dftbp_accuracy, only: dp
  use dftbp_ranlux, only : TRanlux
  use dftbp_charmanip, only : tolower, i2c
  use dftbp_constants, only : elementSymbol

  use fnet_nestedtypes, only : TRealArray2D, TMultiLayerStruc, TMultiLayerStruc_init, TBiasDerivs,&
      & TWeightDerivs, TDerivs, TDerivs_init, TIntArray1D, TPredicts, TPredicts_init, TEnv,&
      & TJacobian, TJacobians
  use fnet_hdf5fx, only : h5ltfx_read_dataset_int_f, h5ltfx_read_dataset_double_f,&
      & h5ltfxmake_dataset_int_f, h5ltfxmake_dataset_double_f
  use fnet_features, only : TFeatures
  use fnet_network, only : TNetwork, TNetwork_init
  use fnet_loss, only : lossFunc, lossGradientFunc, reguFunc, TRegularizationBlock
  use fnet_optimizers, only : TOptimizer, next
  use fnet_fnetdata, only : TDataset
  use fnet_random, only : knuthShuffle

#:if WITH_MPI
  use fnet_mpifx
  use fnet_parallel, only : getStartAndEndIndex
#:endif

  implicit none

  private

  public :: TBpnn, TBpnn_init


  type :: TBpnn

    !> number of sub-nn's
    integer :: nSpecies

    !> atomic numbers of sub-nn species
    integer, allocatable :: atomicNumbers(:)

    !> total number of bias and weight parameters of each sub-nn
    integer :: nBiases, nWeights

    !> number of system-wide training targets of BPNN
    integer :: nGlobalTargets

    !> number of atomic training targets of BPNN
    integer :: nAtomicTargets

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
    procedure :: sysTrain => TBpnn_sysTrain
    procedure :: updateGradients => TBpnn_updateGradients
    procedure :: nTrain => TBpnn_nTrain
    procedure :: iPredict => TBpnn_iPredict
    procedure :: predictBatch => TBpnn_predictBatch
    procedure :: iJacobian => TBpnn_iJacobian
    procedure :: nJacobian => TBpnn_nJacobian
    procedure :: fromFile => TBpnn_fromFile
    procedure :: toFile => TBpnn_toFile

  #:if WITH_MPI
    procedure :: sync => TBpnn_sync
  #:endif

  end type TBpnn


contains

  !> Initialises a BPNN instance.
  subroutine TBpnn_init(this, dims, nSpecies, nGlobalTargets, nAtomicTargets, atomicNumbers,&
      & rndGen, activation)

    !> representation of a Behler-Parrinello neural network
    type(TBpnn), intent(out) :: this

    !> dimensions of all layers in the sub-nn's
    integer, intent(in) :: dims(:)

    !> number of different species in training data
    integer, intent(in) :: nSpecies

    !> number of system-wide training targets of BPNN
    integer, intent(in) :: nGlobalTargets

    !> number of atomic training targets of BPNN
    integer, intent(in) :: nAtomicTargets

    !> atomic numbers of sub-nn species
    integer, intent(in) :: atomicNumbers(:)

    !> luxury pseudorandom generator instance
    type(TRanlux), intent(inout), optional :: rndGen

    !> type of activation function to use for all layers, except output layer (linear)
    character(len=*), intent(in), optional :: activation

    !> Species identifier
    integer :: iSpecies

    this%dims = dims
    this%nSpecies = nSpecies
    this%atomicNumbers = atomicNumbers

    this%nGlobalTargets = nGlobalTargets
    this%nAtomicTargets = nAtomicTargets

    allocate(this%nets(nSpecies))

    do iSpecies = 1, nSpecies
      call TNetwork_init(this%nets(iSpecies), dims, rndGen=rndGen, descriptor=activation)
    end do

    ! assume that every sub-nn has the same architecture
    this%nBiases = this%nets(1)%nBiases
    this%nWeights = this%nets(1)%nWeights

  end subroutine TBpnn_init


  !> Performs multiple training iterations of a BPNN.
  subroutine TBpnn_nTrain(this, env, rndGen, pOptimizer, trainDataset, validDataset, features,&
      & nTrainIt, nPrintOut, nSaveNet, netstatpath, loss, lossgrad, reguLoss, regu, tShuffle,&
      & tMonitorValid, tConverged, trainLoss, validLoss, gradients)

    !> representation of a Behler-Parrinello neural network
    class(TBpnn), intent(inout) :: this

    !> contains mpi communicator, if compiled with mpi enabled
    type(TEnv), intent(in) :: env

    !> luxury pseudorandom generator instance
    type(TRanlux), intent(inout) :: rndGen

    !> general function optimizer
    type(TOptimizer), intent(inout) :: pOptimizer(:)

    !> representation of training and validation dataset
    type(TDataset), intent(in) :: trainDataset, validDataset

    !> collected features, including external and ACSF features
    type(TFeatures), intent(in) :: features

    !> maximum number of training iterations
    integer, intent(in) :: nTrainIt

    !> printout loss/gradient information every nPrintOut steps
    integer, intent(in) :: nPrintOut

    !> save network status every nSaveNet steps
    integer, intent(in) :: nSaveNet

    !> filename of netstat file
    character(len=*), intent(in) :: netstatpath

    !> loss function procedure
    procedure(lossFunc), intent(in), pointer :: loss

    !> procedure, pointing to the choosen loss function gradient
    procedure(lossGradientFunc), intent(in), pointer :: lossgrad

    !> regularization loss function procedure
    procedure(reguFunc), intent(in), pointer :: reguLoss

    !> contains variables of the regularization
    type(TRegularizationBlock), intent(in) :: regu

    !> true, if a Knuth-shuffle should be applied to the gradients calculations
    logical, intent(in) :: tShuffle

    !> true, if validation monitoring is desired
    logical, intent(in) :: tMonitorValid

    !> true, if gradient got below the specified tolerance
    logical, intent(out) :: tConverged

    !> optional, iteration-resolved total loss
    real(dp), intent(out), allocatable, optional :: trainLoss(:)

    !> optional, iteration-resolved total validation loss
    real(dp), intent(out), allocatable, optional :: validLoss(:)

    !> optional, iteration-resolved total gradients
    real(dp), intent(out), allocatable, optional :: gradients(:)

    !> network predictions during the training
    type(TPredicts) :: predicts, resPredicts, validPredicts

    !> total weight and bias gradients of the current system
    type(TDerivs) :: dd, ddRes

    !> temporary (validation) loss function container
    real(dp), allocatable :: tmpLoss(:), tmpValidLoss(:)

    !> temporary iteration-resolved euclidean norm of total gradient
    real(dp), allocatable :: tmpGradients(:)

    !> shuffle array to randomize the order of gradient calculations
    integer, allocatable :: shuffle(:)

    !> true, if current process is the lead
    logical :: tLead

    !> auxiliary variables
    logical :: tPrintOut, tSaveNet

    !> auxiliary variables
    integer :: ii, iIter, iStart, iEnd, iTmpIter, iLastIter

    call TDerivs_init(this%dims, size(this%nets), dd)
    call TDerivs_init(this%dims, size(this%nets), ddRes)

    call TPredicts_init(predicts, trainDataset%nDatapoints, this%nGlobalTargets,&
        & this%nAtomicTargets, trainDataset%localAtToAtNum)
    call TPredicts_init(resPredicts, trainDataset%nDatapoints, this%nGlobalTargets,&
        & this%nAtomicTargets, trainDataset%localAtToAtNum)
    if (tMonitorValid) then
      call TPredicts_init(validPredicts, validDataset%nDatapoints, this%nGlobalTargets,&
          & this%nAtomicTargets, validDataset%localAtToAtNum)
    end if

    allocate(tmpLoss(nTrainIt))
    allocate(tmpGradients(nTrainIt))

    if (tMonitorValid) then
      allocate(tmpValidLoss(nTrainIt))
      tmpValidLoss(:) = 0.0_dp
    end if

  #:if WITH_MPI
    tLead = env%globalMpiComm%lead
    call getStartAndEndIndex(trainDataset%nDatapoints, env%globalMpiComm%size,&
        & env%globalMpiComm%rank, iStart, iEnd)
    call syncWeightsAndBiases(this, env%globalMpiComm)
  #:else
    tLead = .true.
    iStart = 1
    iEnd = trainDataset%nDatapoints
  #:endif

    ! generate shuffle array, if desired
    allocate(shuffle(trainDataset%nDatapoints))
    shuffle = [(ii, ii = 1, trainDataset%nDatapoints)]
    if (tLead .and. tShuffle) then
      call knuthShuffle(rndGen, shuffle)
    end if

  #:if WITH_MPI
    call mpifx_bcast(env%globalMpiComm, shuffle)
  #:endif

    call this%updateGradients(env, iStart, iEnd, shuffle, trainDataset, features%trainFeatures,&
        & lossgrad, dd, ddRes, predicts, resPredicts)

    if (tLead) then
      tmpLoss(1) = loss(resPredicts, trainDataset%globalTargets, trainDataset%atomicTargets,&
          & trainDataset%atomicWeights, weights=trainDataset%weights)
    end if
    if (tMonitorValid) then
      validPredicts%sys = this%predictBatch(features%validFeatures, env,&
          & validDataset%localAtToGlobalSp)
      if (tLead) then
        tmpValidLoss(1) = loss(validPredicts, validDataset%globalTargets,&
            & validDataset%atomicTargets, validDataset%atomicWeights)
      end if
    end if

    lpIter: do iIter = 1, nTrainIt

      iTmpIter = iIter

      if (tLead) then
        call this%update(pOptimizer, ddRes, regu, reguLoss, tmpLoss(iIter),&
            & sum(trainDataset%weights), tmpGradients(iIter), tConverged)
      end if

    #:if WITH_MPI
      call mpifx_bcast(env%globalMpiComm, tConverged)
      call syncWeightsAndBiases(this, env%globalMpiComm)
    #:endif

      call ddRes%reset()

      if (tLead .and. tShuffle) then
        call knuthShuffle(rndGen, shuffle)
      end if

    #:if WITH_MPI
      call mpifx_bcast(env%globalMpiComm, shuffle)
    #:endif

      call this%updateGradients(env, iStart, iEnd, shuffle, trainDataset, features%trainFeatures,&
          & lossgrad, dd, ddRes, predicts, resPredicts)

      if (tLead) then
        tmpLoss(iIter) = loss(resPredicts, trainDataset%globalTargets, trainDataset%atomicTargets,&
            & trainDataset%atomicWeights, weights=trainDataset%weights)
      end if
      if (tMonitorValid) then
        validPredicts%sys = this%predictBatch(features%validFeatures, env,&
            & validDataset%localAtToGlobalSp)
        if (tLead) then
          tmpValidLoss(iIter) = loss(validPredicts, validDataset%globalTargets,&
              & validDataset%atomicTargets, validDataset%atomicWeights)
        end if
      end if

      tPrintOut = (modulo(iIter, nPrintOut) == 0) .or. (iIter == nTrainIt)

      if (tPrintOut) then
        if (tMonitorValid) then
          write(stdout, '(I10,5X,E15.6,4X,E15.6,4X,E15.6)') iIter, tmpLoss(iIter),&
              & tmpGradients(iIter), tmpValidLoss(iIter)
        else
          write(stdout, '(I10,5X,E15.6,4X,E15.6)') iIter, tmpLoss(iIter), tmpGradients(iIter)
        end if
      end if

      tSaveNet = (modulo(iIter, nSaveNet) == 0) .or. (iIter == nTrainIt)

      if (tLead .and. tSaveNet) then
        call this%toFile(netstatpath)
      end if

      if (tConverged) then
        if (tMonitorValid) then
          write(stdout, '(I10,5X,E15.6,4X,E15.6,4X,E15.6)') iIter, tmpLoss(iIter),&
              & tmpGradients(iIter), tmpValidLoss(iIter)
        else
          write(stdout, '(I10,5X,E15.6,4X,E15.6)') iIter, tmpLoss(iIter), tmpGradients(iIter)
        end if
        if (tLead) then
          call this%toFile(netstatpath)
        end if
        exit lpIter
      end if

    end do lpIter

  #:if WITH_MPI
    call mpifx_bcast(env%globalMpiComm, tmpLoss)
    if (tMonitorValid) then
      call mpifx_bcast(env%globalMpiComm, tmpValidLoss)
    end if
    call mpifx_bcast(env%globalMpiComm, tmpGradients)
    call mpifx_allreduce(env%globalMpiComm, iTmpIter, iLastIter, MPI_MAX)
  #:else
    iLastIter = iTmpIter
  #:endif

    ! crop loss to actual range of past iterations and, if desired, pass through
    if (present(trainLoss)) then
      allocate(trainLoss(iLastIter))
      trainLoss(:) = tmpLoss(1:iLastIter)
    end if
    if (present(gradients)) then
      allocate(gradients(iLastIter))
      gradients(:) = tmpGradients(1:iLastIter)
    end if
    if (tMonitorValid .and. present(validLoss)) then
      allocate(validLoss(iLastIter))
      validLoss(:) = tmpValidLoss(1:iLastIter)
    end if

  end subroutine TBpnn_nTrain


  !> Calculates new gradients in weight-bias space based on the current BPNN state.
  subroutine TBpnn_updateGradients(this, env, iStart, iEnd, shuffle, trainDataset, trainFeatures,&
      & lossgrad, dd, resDd, predicts, resPredicts)

    !> representation of a Behler-Parrinello neural network
    class(TBpnn), intent(inout) :: this

    !> contains mpi communicator, if compiled with mpi enabled
    type(TEnv), intent(in) :: env

    !> system index to start and end at
    integer, intent(in) :: iStart, iEnd

    !> shuffle indices for permuting the gradient calculation
    integer, intent(in) :: shuffle(:)

    !> training dataset representation
    type(TDataset), intent(in) :: trainDataset

    !> collection of input features for training
    type(TRealArray2D), intent(in) :: trainFeatures(:)

    !> procedure, pointing to the choosen loss function gradient
    procedure(lossGradientFunc), intent(in), pointer :: lossgrad

    !> total weight and bias gradients of the current system
    type(TDerivs), intent(inout) :: dd

    !> resulting total weight and bias gradients
    type(TDerivs), intent(inout) :: resDd

    !> network predictions during the training
    type(TPredicts), intent(inout) :: predicts

    !> resulting network predictions
    type(TPredicts), intent(inout) :: resPredicts

    !> temporary weight and bias gradient storage of the current atom
    type(TDerivs) :: ddTmp

    !> auxiliary variables
    integer :: ii, iSys, iGlobalSp, iLayer

    lpSystem: do ii = iStart, iEnd
      iSys = shuffle(ii)
      call this%sysTrain(trainFeatures(iSys)%array, trainDataset%globalTargets(iSys)%array,&
          & trainDataset%atomicTargets(iSys)%array, trainDataset%atomicWeights(iSys)%array,&
          & trainDataset%localAtToGlobalSp(iSys)%array, lossgrad, ddTmp,&
          & predicts%sys(iSys)%array)
      ! collect gradients of the system
      do iGlobalSp = 1, size(this%nets)
        do iLayer = 1, size(this%dims)
          dd%dw(iGlobalSp)%dw(iLayer)%array = dd%dw(iGlobalSp)%dw(iLayer)%array +&
              & ddTmp%dw(iGlobalSp)%dw(iLayer)%array * real(trainDataset%weights(iSys), dp)
          dd%db(iGlobalSp)%db(iLayer)%array = dd%db(iGlobalSp)%db(iLayer)%array +&
              & ddTmp%db(iGlobalSp)%db(iLayer)%array * real(trainDataset%weights(iSys), dp)
        end do
      end do
    end do lpSystem

  #:if WITH_MPI
    ! collect outputs of all nodes
    do iSys = 1, trainDataset%nDatapoints
      call mpifx_allreduce(env%globalMpiComm, predicts%sys(iSys)%array,&
          & resPredicts%sys(iSys)%array, MPI_SUM)
    end do
    ! sync gradients between MPI nodes
    do iGlobalSp = 1, size(this%nets)
      do iLayer = 1, size(this%dims)
        call mpifx_allreduce(env%globalMpiComm, dd%dw(iGlobalSp)%dw(iLayer)%array,&
            & resDd%dw(iGlobalSp)%dw(iLayer)%array, MPI_SUM)
        call mpifx_allreduce(env%globalMpiComm, dd%db(iGlobalSp)%db(iLayer)%array,&
            & resDd%db(iGlobalSp)%db(iLayer)%array, MPI_SUM)
      end do
    end do
  #:else
    resPredicts = predicts
    resDd = dd
  #:endif

    ! reset the temporary predictions, important if random shuffling is applied
    do iSys = 1, size(predicts%sys)
      predicts%sys(iSys)%array(:,:) = 0.0_dp
    end do

    ! reset the temporary gradient storage type
    call dd%reset()

  end subroutine TBpnn_updateGradients


#:if WITH_MPI

  !> Synchronizes the entire BPNN derived type.
  subroutine TBpnn_sync(this, comm)

    !> representation of a Behler-Parrinello neural network
    class(TBpnn), intent(inout) :: this

    !> mpi communicator with some additional information
    type(mpifx_comm), intent(in) :: comm

    !> auxiliary variables
    integer :: iNet, iLayer, dims0d
    integer, allocatable :: tmpDims(:)
    character(len=1000) :: tmpStr
    character(len=:), allocatable :: tmpTransferType

    call mpifx_bcast(comm, this%nSpecies)
    call mpifx_bcast(comm, this%nBiases)
    call mpifx_bcast(comm, this%nWeights)

    if (comm%lead) then
      dims0d = size(this%atomicNumbers)
    end if
    call mpifx_bcast(comm, dims0d)
    if (.not. comm%lead) then
      if (allocated(this%atomicNumbers)) deallocate(this%atomicNumbers)
      allocate(this%atomicNumbers(dims0d))
    end if
    call mpifx_bcast(comm, this%atomicNumbers)

    if (comm%lead) then
      dims0d = size(this%dims)
    end if
    call mpifx_bcast(comm, dims0d)
    if (.not. comm%lead) then
      if (allocated(this%dims)) deallocate(this%dims)
      allocate(this%dims(dims0d))
    end if
    call mpifx_bcast(comm, this%dims)

    if (comm%lead) then
      dims0d = size(this%nets)
    end if
    call mpifx_bcast(comm, dims0d)
    if (.not. comm%lead) then
      if (allocated(this%nets)) deallocate(this%nets)
      allocate(this%nets(dims0d))
    end if

    do iNet = 1, size(this%nets)
      call mpifx_bcast(comm, this%nets(iNet)%nBiases)
      call mpifx_bcast(comm, this%nets(iNet)%nWeights)
      ! provide network dimensions to non-leads
      if (comm%lead) then
        dims0d = size(this%nets(iNet)%dims)
        tmpDims = this%nets(iNet)%dims
      end if
      call mpifx_bcast(comm, dims0d)
      if (.not. comm%lead) then
        if (allocated(tmpDims)) deallocate(tmpDims)
        allocate(tmpDims(dims0d))
      end if
      call mpifx_bcast(comm, tmpDims)
      ! provide network transfer function
      if (comm%lead) then
        dims0d = len(this%nets(iNet)%transferType)
        tmpStr = this%nets(iNet)%transferType
      end if
      call mpifx_bcast(comm, dims0d)
      if (.not. comm%lead) then
        if (allocated(tmpTransferType)) deallocate(tmpTransferType)
      end if
      call mpifx_bcast(comm, tmpStr)
      tmpTransferType = tolower(trim(tmpStr))
      ! synchronize network layers
      if (.not. comm%lead) then
        if (allocated(this%nets(iNet)%dims)) deallocate(this%nets(iNet)%dims)
        if (allocated(this%nets(iNet)%layers)) deallocate(this%nets(iNet)%layers)
        if (allocated(this%nets(iNet)%transferType)) deallocate(this%nets(iNet)%transferType)
        call TNetwork_init(this%nets(iNet), tmpDims, descriptor=tmpTransferType)
      end if
    end do

    ! broadcast actual layer-resolved values
    do iNet = 1, size(this%nets)
      do iLayer = 1, size(this%nets(iNet)%dims)
        call mpifx_bcast(comm, this%nets(iNet)%layers(iLayer)%aa)
        call mpifx_bcast(comm, this%nets(iNet)%layers(iLayer)%aarg)
      end do
      do iLayer = 1, size(this%nets(iNet)%dims) - 1
        call mpifx_bcast(comm, this%nets(iNet)%layers(iLayer + 1)%bb)
        call mpifx_bcast(comm, this%nets(iNet)%layers(iLayer)%ww)
      end do
    end do

  end subroutine TBpnn_sync


  !> Synchronizes the BPNN network parameters between the MPI nodes.
  subroutine syncWeightsAndBiases(this, comm)

    !> representation of a Behler-Parrinello neural network
    type(TBpnn), intent(inout) :: this

    !> mpi communicator with some additional information
    type(mpifx_comm), intent(in) :: comm

    !> auxiliary variables
    integer :: iGlobalSp, iLayer

    do iGlobalSp = 1, size(this%nets)
      ! broadcast all weights and biases of a sub-network
      do iLayer = 1, size(this%dims)
        call mpifx_bcast(comm, this%nets(iGlobalSp)%layers(iLayer)%ww)
        call mpifx_bcast(comm, this%nets(iGlobalSp)%layers(iLayer)%bb)
      end do
    end do

  end subroutine syncWeightsAndBiases
#:endif


  !> Performs a training iteration for a single system, i.e. set of input features.
  subroutine TBpnn_sysTrain(this, input, globalTargets,&
      & atomicTargets, atomicWeights, localAtToGlobalSp, lossgrad, dd, predicts)

    !> representation of a Behler-Parrinello neural network
    class(TBpnn), intent(inout) :: this

    !> features of all atoms of the current system/geometry, shape: [nFeatures, nAtoms]
    real(dp), intent(in) :: input(:,:)

    !> system-wide target data, shape: [nGlobalTargets]
    real(dp), intent(in) :: globalTargets(:)

    !> atomic target data, shape: [nAtomicTargets, nAtoms]
    real(dp), intent(in) :: atomicTargets(:,:)

    !> contains atomic gradient weights, exp. shape: [nAtoms]
    real(dp), intent(in) :: atomicWeights(:)

    !> maps local atom index --> global species index
    integer, intent(in) :: localAtToGlobalSp(:)

    !> procedure, pointing to the choosen loss function gradient
    procedure(lossGradientFunc), intent(in), pointer :: lossgrad

    !> total weight and bias gradients of the current system
    type(TDerivs), intent(out) :: dd

    !> summed (atom-resolved) predictions of all sub-nn's
    !> shape: [nGlobalTargets + nAtomicTargets, nAtoms]
    real(dp), intent(out) :: predicts(:,:)

    !> representation of temporary layer storage container
    type(TMultiLayerStruc) :: tmpLayer

    !> temporary weight gradient storage of the current atom
    type(TWeightDerivs) :: dwTmp

    !> temporary bias gradient storage of the current atom
    type(TBiasDerivs) :: dbTmp

    !> loss gradients w.r.t predictions and targets (only system-wide targets)
    real(dp), allocatable :: globalLossGrad(:)

    !> loss gradients w.r.t predictions and targets
    real(dp), allocatable :: lossgrads(:,:)

    !> temporary output storage of sub-nn's
    real(dp), allocatable :: tmpOut(:)

    !> temporary real valued storage for summed up system-wide predictions
    real(dp), allocatable :: globalPredicts(:)

    !> auxiliary variable
    integer :: iAtom, iGlobalSp, iLayer

    call TMultiLayerStruc_init(this%dims, size(input, dim=2), tmpLayer)
    call TDerivs_init(this%dims, size(this%nets), dd)

    allocate(lossgrads, mold=predicts)

    do iAtom = 1, size(input, dim=2)
      iGlobalSp = localAtToGlobalSp(iAtom)
      call this%nets(iGlobalSp)%fprop(input(:, iAtom), state=tmpLayer%struc(iAtom)%layers,&
          & out=tmpOut)
      predicts(:, iAtom) = tmpOut
    end do

    if (this%nGlobalTargets > 0) then
      globalPredicts = sum(predicts(1:this%nGlobalTargets, :), dim=2)
      globalLossGrad = reshape(lossgrad(reshape(globalPredicts, [this%nGlobalTargets, 1]),&
          & reshape(globalTargets, [this%nGlobalTargets, 1])), [this%nGlobalTargets])
      do iAtom = 1, size(input, dim=2)
        lossgrads(1:this%nGlobalTargets, iAtom) = globalLossGrad
      end do
    end if

    if (this%nAtomicTargets > 0) then
      lossgrads(this%nGlobalTargets + 1:, :) = lossgrad(predicts(this%nGlobalTargets + 1:, :),&
          & atomicTargets)
    end if

    ! collect gradients and normalize them by atom number
    do iAtom = 1, size(input, dim=2)
      iGlobalSp = localAtToGlobalSp(iAtom)
      this%nets(iGlobalSp)%layers = tmpLayer%struc(iAtom)%layers
      call this%nets(iGlobalSp)%bprop(lossgrads(:, iAtom), dwTmp, dbTmp)
      do iLayer = 1, size(this%dims)
        dd%dw(iGlobalSp)%dw(iLayer)%array = dd%dw(iGlobalSp)%dw(iLayer)%array +&
            & dwTmp%dw(iLayer)%array * atomicWeights(iAtom) / real(size(input, dim=2), dp)
        dd%db(iGlobalSp)%db(iLayer)%array = dd%db(iGlobalSp)%db(iLayer)%array +&
            & dbTmp%db(iLayer)%array * atomicWeights(iAtom) / real(size(input, dim=2), dp)
      end do
    end do

  end subroutine TBpnn_sysTrain


  !> Updates the BPNN parameters based on the obtained gradients in weight-bias space.
  subroutine TBpnn_update(this, pOptimizer, dd, regu, reguLoss, loss, nDatapoints, totGradNorm,&
      & tConverged)

    !> representation of a Behler-Parrinello neural network
    class(TBpnn), intent(inout) :: this

    !> general function optimizer
    type(TOptimizer), intent(inout) :: pOptimizer(:)

    !> weight and bias gradients
    type(TDerivs), intent(inout) :: dd

    !> contains variables of the regularization
    type(TRegularizationBlock), intent(in) :: regu

    !> regularization loss function procedure
    procedure(reguFunc), intent(in), pointer :: reguLoss

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

    !> sub-network resolved loss values
    real(dp), allocatable :: subnetLoss(:)

    !> auxiliary variable
    integer :: iGlobalSp

    ! add loss based regularization, if desired
    do iGlobalSp = 1, size(dd%dw)
      call dd%dw(iGlobalSp)%elasticNetRegularization(this%nets(iGlobalSp)%layers, this%nWeights,&
          & regu%strength, regu%alpha)
    end do

    call dd%serialized(this%nBiases + this%nWeights, ddSerial)
    call this%serializedWeightsAndBiases(weightsAndBiases)

    ! add loss based regularization, if desired
    allocate(subnetLoss(size(this%nets)))
    do iGlobalSp = 1, size(this%nets)
      subnetLoss(iGlobalSp) = loss + reguLoss(weightsAndBiases(1:this%nWeights, iGlobalSp),&
          & regu%strength, regu%alpha)
    end do

    ddSerial = ddSerial / real(nDatapoints, dp)
    totGradNorm = norm2(ddSerial)

    allocate(newWeightsAndBiases(size(weightsAndBiases, dim=1), size(weightsAndBiases, dim=2)))

    do iGlobalSp = 1, size(this%nets)
      call next(pOptimizer(iGlobalSp), subnetLoss(iGlobalSp), ddSerial(:, iGlobalSp),&
          & newWeightsAndBiases(:, iGlobalSp), tConverged)
    end do

    call this%serialWeightsAndBiasesFillup(newWeightsAndBiases)

  end subroutine TBpnn_update


  !> Inserts serialized network-resolved weights and biases into the BPNN structure.
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


  !> Maps the BPNN parameters to a serialized array, suitable for the gradient-based optimizer.
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


  !> Resets the activation values of a BPNN instance.
  subroutine TBpnn_resetActivations(this)

    !> representation of a Behler-Parrinello neural network
    class(TBpnn), intent(inout) :: this

    !> Species identifier
    integer :: iSpecies

    do iSpecies = 1, this%nSpecies
      call this%nets(iSpecies)%resetActivations()
    end do

  end subroutine TBpnn_resetActivations


  !> Collects the outputs of all sub-networks (for predicting system-wide properties).
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


  function TBpnn_iPredict(this, input, localAtToGlobalSp) result(atomicOutput)

    !> representation of a Behler-Parrinello neural network
    class(TBpnn), intent(in) :: this

    !> single sample of input data to calculate output for
    real(dp), intent(in) :: input(:,:)

    !> index mapping local atom --> global species index
    integer, intent(in) :: localAtToGlobalSp(:)

    !> atomic contributions to the network prediction
    real(dp), allocatable :: atomicOutput(:,:)

    !> number of predictions for each sub-nn
    integer :: nOutput

    !> number of atoms in the current system
    integer :: nAtoms

    !> auxiliary variables
    integer :: iAtom, iGlobalSp

    nAtoms = size(input, dim=2)
    nOutput = this%dims(size(this%dims))

    allocate(atomicOutput(nOutput, nAtoms))

    do iAtom = 1, nAtoms
      iGlobalSp = localAtToGlobalSp(iAtom)
      atomicOutput(:, iAtom) = this%nets(iGlobalSp)%iPredict(input(:, iAtom))
    end do

  end function TBpnn_iPredict


  !> Calculates the Jacobian matrix element of the atomic outputs w.r.t. the input features.
  function TBpnn_iJacobian(this, input, localAtToGlobalSp) result(jacobian)

    !> representation of a Behler-Parrinello neural network
    class(TBpnn), intent(in) :: this

    !> single sample of input data to calculate Jacobian matrices for
    real(dp), intent(in) :: input(:,:)

    !> index mapping local atom --> global species index
    integer, intent(in) :: localAtToGlobalSp(:)

    !> derivatives of outputs w.r.t. atomic input features
    type(TJacobian) :: jacobian

    !> number of atoms in the current system
    integer :: nAtoms

    !> auxiliary variables
    integer :: iAtom, iGlobalSp

    nAtoms = size(input, dim=2)

    allocate(jacobian%atom(nAtoms))

    do iAtom = 1, nAtoms
      iGlobalSp = localAtToGlobalSp(iAtom)
      jacobian%atom(iAtom)%array = this%nets(iGlobalSp)%fdevi(input(:, iAtom))
    end do

  end function TBpnn_iJacobian


  !> Calculates the network's Jacobian for a batch of input features.
  function TBpnn_nJacobian(this, features, env, localAtToGlobalSp) result(resJacobian)

    !> representation of a Behler-Parrinello neural network
    class(TBpnn), intent(in) :: this

    !> atomic features as network input
    type(TRealArray2D), intent(in) :: features(:)

    !> if compiled with mpi enabled, contains mpi communicator
    type(TEnv), intent(in) :: env

    !> index mapping local atom --> global species index
    type(TIntArray1D), intent(in) :: localAtToGlobalSp(:)

    !> network predictions during the training
    type(TJacobians) :: jacobian, resJacobian

    !> true, if current process is the lead
    logical :: tLead

    !> auxiliary variables
    integer :: iSys, iAtom, iStart, iEnd

  #:if WITH_MPI
    tLead = env%globalMpiComm%lead
    call getStartAndEndIndex(size(features), env%globalMpiComm%size, env%globalMpiComm%rank,&
        & iStart, iEnd)
  #:else
    tLead = .true.
    iStart = 1
    iEnd = size(features)
  #:endif

    allocate(jacobian%sys(size(features)))

    do iSys = 1, size(features)
      allocate(jacobian%sys(iSys)%atom(size(features(iSys)%array, dim=2)))
      do iAtom = 1, size(features(iSys)%array, dim=2)
        allocate(jacobian%sys(iSys)%atom(iAtom)%array(this%dims(size(this%dims)), this%dims(1)))
        jacobian%sys(iSys)%atom(iAtom)%array(:,:) = 0.0_dp
      end do
    end do

    resJacobian = jacobian

    do iSys = iStart, iEnd
      jacobian%sys(iSys) = this%iJacobian(features(iSys)%array, localAtToGlobalSp(iSys)%array)
    end do

  #:if WITH_MPI
    do iSys = 1, size(features)
      do iAtom = 1, size(features(iSys)%array, dim=2)
        call mpifx_allreduce(env%globalMpiComm, jacobian%sys(iSys)%atom(iAtom)%array,&
            & resJacobian%sys(iSys)%atom(iAtom)%array, MPI_SUM)
      end do
    end do
  #:else
    resJacobian = jacobian
  #:endif

  end function TBpnn_nJacobian


  !> Calculates the network output for a batch of input features.
  function TBpnn_predictBatch(this, features, env, localAtToGlobalSp) result(resPredicts)

    !> representation of a Behler-Parrinello neural network
    class(TBpnn), intent(in) :: this

    !> atomic features as network input
    type(TRealArray2D), intent(in) :: features(:)

    !> if compiled with mpi enabled, contains mpi communicator
    type(TEnv), intent(in) :: env

    !> index mapping local atom --> global species index
    type(TIntArray1D), intent(in) :: localAtToGlobalSp(:)

    !> network predictions during the training
    type(TRealArray2D), allocatable :: predicts(:), resPredicts(:)

    !> true, if current process is the lead
    logical :: tLead

    !> auxiliary variables
    integer :: iSys, iStart, iEnd

  #:if WITH_MPI
    tLead = env%globalMpiComm%lead
    call getStartAndEndIndex(size(features), env%globalMpiComm%size, env%globalMpiComm%rank,&
        & iStart, iEnd)
  #:else
    tLead = .true.
    iStart = 1
    iEnd = size(features)
  #:endif

    allocate(predicts(size(features)))

    do iSys = 1, size(features)
      allocate(predicts(iSys)%array(this%dims(size(this%dims)), size(features(iSys)%array,&
          & dim=2)))
      predicts(iSys)%array(:,:) = 0.0_dp
    end do

    resPredicts = predicts

    do iSys = iStart, iEnd
      predicts(iSys)%array(:,:) = this%iPredict(features(iSys)%array,&
          & localAtToGlobalSp(iSys)%array)
    end do

  #:if WITH_MPI
    do iSys = 1, size(features)
      call mpifx_allreduce(env%globalMpiComm, predicts(iSys)%array,&
          & resPredicts(iSys)%array, MPI_SUM)
    end do
  #:else
    resPredicts = predicts
  #:endif

  end function TBpnn_predictBatch


  !> Writes weight and bias parameters to netstat file.
  subroutine TBpnn_toFile(this, fname)

    !> representation of a Behler-Parrinello neural network
    class(TBpnn), intent(in) :: this

    !> filename (will be fortnet.hdf5)
    character(len=*), intent(in) :: fname

    !> various specifier flags
    integer(hid_t) :: file_id, netstat_id, bpnn_id, subnet_id, layer_id

    !> name of current subnetwork and layer
    character(len=:), allocatable :: netname, layername

    !> auxiliary variables
    integer :: iErr, iNet, iLayer, tExist

    ! open the hdf5 interface
    call h5open_f(iErr)

    ! open the netstat file
    call h5fopen_f(fname, H5F_ACC_RDWR_F, file_id, iErr)

    ! open the netstat group
    call h5gopen_f(file_id, 'netstat', netstat_id, iErr)

    ! open the bpnn group
    call h5gopen_f(netstat_id, 'bpnn', bpnn_id, iErr)

    do iNet = 1, this%nSpecies
      netname = trim(elementSymbol(this%atomicNumbers(iNet))) // '-subnetwork'

      ! open the subnetwork group
      call h5gopen_f(bpnn_id, netname, subnet_id, iErr)

      do iLayer = 1, size(this%nets(iNet)%dims) - 1
        layername = 'layer' // i2c(iLayer)

        ! open the layer group
        call h5gopen_f(subnet_id, layername, layer_id, iErr)

        tExist = h5ltfind_dataset_f(layer_id, 'weights')

        if (tExist == 1) then
          call h5ldelete_f(layer_id, 'weights', iErr)
        end if

        tExist = h5ltfind_dataset_f(layer_id, 'bias')

        if (tExist == 1) then
          call h5ldelete_f(layer_id, 'bias', iErr)
        end if

        ! write weight matrices
        call h5ltfxmake_dataset_double_f(layer_id, 'weights', this%nets(iNet)%layers(iLayer)%ww)

        ! write bias vectors
        call h5ltfxmake_dataset_double_f(layer_id, 'bias', this%nets(iNet)%layers(iLayer + 1)%bb)

        ! close the layer group
        call h5gclose_f(layer_id, iErr)

      end do

      ! close the subnetwork group
      call h5gclose_f(subnet_id, iErr)

    end do

    ! close the bpnn group
    call h5gclose_f(bpnn_id, iErr)

    ! close the netstat group
    call h5gclose_f(netstat_id, iErr)

    ! close the netstat file
    call h5fclose_f(file_id, iErr)

    ! close the hdf5 interface
    call h5close_f(iErr)

  end subroutine TBpnn_toFile


  !> Reads a BPNN from a netstat file.
  subroutine TBpnn_fromFile(this, fname)

    !> representation of a Behler-Parrinello neural network
    class(TBpnn), intent(out) :: this

    !> filename (will be fortnet.hdf5)
    character(len=*), intent(in) :: fname

    !> various specifier flags
    integer(hid_t) :: file_id, netstat_id, bpnn_id, subnet_id, layer_id

    !> name of current subnetwork and layer
    character(len=:), allocatable :: netname, layername

    !> temporary activation function type
    character(len=:), allocatable :: activation, tmpActivation

    !> dimensions of all layers in the sub-nn's
    integer, allocatable :: dims(:), tmpDims(:)

    !> atomic numbers of sub-nn species
    integer, allocatable :: atomicNumbers(:)

    !> number of system-wide training targets of BPNN
    integer :: nGlobalTargets

    !> number of atomic training targets of BPNN
    integer :: nAtomicTargets

    !> temporary storage container
    integer :: tmpInt(1)
    character(len=100) :: tmpStr

    !> auxiliary variables
    integer :: iErr, iNet, iLayer

    ! open the hdf5 interface
    call h5open_f(iErr)

    ! open the netstat file
    call h5fopen_f(fname, H5F_ACC_RDONLY_F, file_id, iErr)

    ! open the netstat group
    call h5gopen_f(file_id, 'netstat', netstat_id, iErr)

    ! open the bpnn group
    call h5gopen_f(netstat_id, 'bpnn', bpnn_id, iErr)

    ! read atomic numbers of sub-nn
    call h5ltfx_read_dataset_int_f(bpnn_id, 'atomicnumbers', atomicNumbers)

    ! read the output type, i.e. number of atomic and global targets during training
    call h5ltget_attribute_int_f(bpnn_id, './', 'nglobaltargets', tmpInt, iErr)
    nGlobalTargets = tmpInt(1)
    call h5ltget_attribute_int_f(bpnn_id, './', 'natomictargets', tmpInt, iErr)
    nAtomicTargets = tmpInt(1)

    if (nGlobalTargets < 0) then
      call error('Error while reading BPNN from netstat file.&
          & Negative number of system-wide targets.')
    end if

    if (nAtomicTargets < 0) then
      call error('Error while reading BPNN from netstat file.&
          & Negative number of atomic targets.')
    end if

    if ((nGlobalTargets + nAtomicTargets) == 0) then
      call error('Error while reading BPNN from netstat file. Total number of targets is zero')
    end if

    do iNet = 1, size(atomicNumbers)
      netname = trim(elementSymbol(atomicNumbers(iNet))) // '-subnetwork'

      ! open the subnetwork group
      call h5gopen_f(bpnn_id, netname, subnet_id, iErr)

      ! read layer topology
      call h5ltfx_read_dataset_int_f(subnet_id, 'topology', tmpDims)
      if (.not. allocated(dims)) then
        dims = tmpDims
      end if
      if (.not. all(shape(dims) == shape(tmpDims))) then
        call error('Error while reading sub-network topology. Currently only equally shaped'&
            & //' networks are supported.')
      end if
      if (.not. all(dims == tmpDims)) then
        call error('Error while reading sub-network topology. Currently only equal topologies'&
            & //' are supported.')
      end if
      dims = tmpDims

      ! read transfer function type
      call h5ltget_attribute_string_f(subnet_id, './', 'activation', tmpStr, iErr)
      tmpActivation = tolower(trim(tmpStr))
      if (.not. allocated(activation)) then
        activation = tmpActivation
      end if
      if (.not. (tmpActivation == activation)) then
        call error('Error while reading sub-network transfer functions. Currently only networks'&
            & //' with equal transfer are supported.')
      end if
      activation = tmpActivation

      ! close the subnetwork group
      call h5gclose_f(subnet_id, iErr)

    end do

    ! allocate sub-nn's and their layer structures
    call TBpnn_init(this, dims, size(atomicNumbers), nGlobalTargets, nAtomicTargets, atomicNumbers,&
        & activation=activation)

    do iNet = 1, size(this%atomicNumbers)
      netname = trim(elementSymbol(this%atomicNumbers(iNet))) // '-subnetwork'

      ! open the subnetwork group
      call h5gopen_f(bpnn_id, netname, subnet_id, iErr)

      do iLayer = 1, size(this%nets(iNet)%dims) - 1
        layername = 'layer' // i2c(iLayer)

        ! open the layer group
        call h5gopen_f(subnet_id, layername, layer_id, iErr)

        ! read weight matrices
        call h5ltfx_read_dataset_double_f(layer_id, 'weights', this%nets(iNet)%layers(iLayer)%ww)

        ! read bias vectors
        call h5ltfx_read_dataset_double_f(layer_id, 'bias', this%nets(iNet)%layers(iLayer + 1)%bb)

        ! close the layer group
        call h5gclose_f(layer_id, iErr)

      end do

      ! close the subnetwork group
      call h5gclose_f(subnet_id, iErr)

    end do

    ! close the bpnn group
    call h5gclose_f(bpnn_id, iErr)

    ! close the netstat group
    call h5gclose_f(netstat_id, iErr)

    ! close the netstat file
    call h5fclose_f(file_id, iErr)

    ! close the hdf5 interface
    call h5close_f(iErr)

  end subroutine TBpnn_fromFile

end module fnet_bpnn
