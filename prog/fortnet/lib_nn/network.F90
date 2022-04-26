!--------------------------------------------------------------------------------------------------!
!  FORTNET: A Behler-Parrinello-Neural-Network Implementation                                      !
!  Copyright (C) 2020 - 2022  T. W. van der Heide                                                  !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Implements a vanilla, dense, feed-forward multilayer-perceptron.
module fnet_network

  use dftbp_assert
  use dftbp_accuracy, only: dp
  use dftbp_ranlux, only : TRanlux
  use dftbp_math_blasroutines, only : gemv, gemm

  use fnet_nestedtypes, only : TBiasDerivs, TBiasDerivs_init
  use fnet_nestedtypes, only : TWeightDerivs, TWeightDerivs_init
  use fnet_layer, only : TLayer, TLayer_init

  implicit none

  private

  public :: TNetwork, TNetwork_init


  type :: TNetwork

    !> total number of bias and weight parameters
    integer :: nBiases, nWeights

    !> dimensions of all layers in the network
    integer, allocatable :: dims(:)

    !> representation of the network layers
    type(TLayer), allocatable :: layers(:)

    !> type of transfer function
    character(len=:), allocatable :: transferType

  contains

    procedure :: initLayers => TNetwork_initLayers
    procedure :: countParams => TNetwork_countParams
    procedure :: fprop => TNetwork_fprop
    procedure :: fdevi => TNetwork_fdevi
    procedure :: bprop => TNetwork_bprop
    procedure :: getOutput => TNetwork_getOutput
    procedure :: iPredict => TNetwork_iPredict
    procedure :: nPredict => TNetwork_nPredict
    procedure :: resetActivations => TNetwork_resetActivations
    procedure :: setTransferFunc => TNetwork_setTransferFunc
    procedure :: serializedWeightsAndBiases => TNetwork_serializedWeightsAndBiases
    procedure :: serialWeightsAndBiasesFillup => TNetwork_serialWeightsAndBiasesFillup

  end type TNetwork

contains

  !> Initialises feed-forward multilayer-perceptron instance.
  subroutine TNetwork_init(this, dims, rndGen, descriptor)

    !> representation of a neural network
    type(TNetwork), intent(inout) :: this

    !> dimensions of all layers in the network
    integer, intent(in) :: dims(:)

    !> luxury pseudorandom generator instance
    type(TRanlux), intent(inout), optional :: rndGen

    !> type of transfer function to use
    character(len=*), intent(in), optional :: descriptor

    this%dims = dims

    call this%initLayers(dims, rndGen=rndGen)

    if (present(descriptor)) then
      call this%setTransferFunc(descriptor)
      this%transferType = descriptor
    else
      call this%setTransferFunc('tanh')
      this%transferType = 'tanh'
    end if

    call this%countParams()

  end subroutine TNetwork_init


  !> Initialises the inherent layer structure.
  subroutine TNetwork_initLayers(this, dims, rndGen)

    !> representation of a neural network
    class(TNetwork), intent(inout) :: this

    !> dimensions of all layers in the network
    integer, intent(in) :: dims(:)

    !> luxury pseudorandom generator instance
    type(TRanlux), intent(inout), optional :: rndGen

    !! auxiliary variable
    integer :: ii

    this%dims = dims

    if (.not. allocated(this%layers)) allocate(this%layers(size(dims)))

    do ii = 1, size(dims) - 1
      call TLayer_init(this%layers(ii), dims(ii), dims(ii + 1), rndGen=rndGen)
    end do

    call TLayer_init(this%layers(size(dims)), dims(size(dims)), 1, rndGen=rndGen)

    this%layers(1)%bb = 0.0_dp
    this%layers(size(dims))%ww = 0.0_dp

  end subroutine TNetwork_initLayers


  !> Counts the weight and bias parameters of the network instance.
  subroutine TNetwork_countParams(this)

    !> representation of a neural network
    class(TNetwork), intent(inout) :: this

    !! auxiliary variable
    integer :: iLayer

    this%nBiases = sum(this%dims)

    this%nWeights = 0
    do iLayer = 1, size(this%dims) - 1
      this%nWeights = this%nWeights + this%dims(iLayer) * this%dims(iLayer + 1)
    end do
    this%nWeights = this%nWeights + this%dims(size(this%dims))

  end subroutine TNetwork_countParams


  !> Propagates input features through the network and stores the activations (+ arguments).
  subroutine TNetwork_fprop(this, xx, state, out)

    !> representation of a neural network
    class(TNetwork), intent(inout) :: this

    !> input to feed network with
    real(dp), intent(in) :: xx(:)

    !> copy of representation of the network layers
    type(TLayer), intent(out), allocatable, optional :: state(:)

    !> outputs of current neural network instance
    real(dp), intent(out), allocatable, optional :: out(:)

    !! temporary matrix * vector storage
    real(dp), allocatable :: tmpVec(:)

    !! auxiliary variable
    integer :: ii

    this%layers(1)%aa = xx

    do ii = 2, size(this%layers)
      allocate(tmpVec(size(transpose(this%layers(ii - 1)%ww), dim=1)))
      tmpVec(:) = 0.0_dp
      call gemv(tmpVec, transpose(this%layers(ii - 1)%ww), this%layers(ii - 1)%aa)
      this%layers(ii)%aarg = tmpVec + this%layers(ii)%bb
      this%layers(ii)%aa = this%layers(ii)%transfer(this%layers(ii)%aarg)
      deallocate(tmpVec)
    end do

    if (present(state)) state = this%layers
    if (present(out)) call this%getOutput(out)

  end subroutine TNetwork_fprop


  function TNetwork_fdevi(this, xx) result(jacobi)

    !> representation of a neural network
    class(TNetwork), intent(in) :: this

    !> input to feed network with
    real(dp), intent(in) :: xx(:)

    !! forward derivatives, i.e. jacobian matrix w.r.t. input features
    real(dp), allocatable :: jacobi(:,:)

    !! activation and network input of neurons of current layer
    real(dp), allocatable :: aa(:), aarg(:)

    !! transfer derivatives
    real(dp), allocatable :: deriv(:), diagDeriv(:,:)

    !! temporary matrix * vector and matrix * matrix storage
    real(dp), allocatable :: tmpVec(:), tmpMat(:,:), tmpMat2(:,:)

    !! auxiliary variables
    integer :: ii, jj

    aa = xx

    allocate(jacobi(size(aa), size(aa)))
    jacobi(:,:) = 0.0_dp

    do ii = 1, size(jacobi, dim=1)
      jacobi(ii, ii) = 1.0_dp
    end do

    do ii = 2, size(this%layers)

      allocate(tmpVec(size(transpose(this%layers(ii - 1)%ww), dim=1)))
      tmpVec(:) = 0.0_dp
      call gemv(tmpVec, transpose(this%layers(ii - 1)%ww), aa)
      aarg = tmpVec + this%layers(ii)%bb
      deallocate(tmpVec)

      aa = this%layers(ii)%transfer(aarg)
      deriv = this%layers(ii)%transferDeriv(aarg)

      if (allocated(diagDeriv)) deallocate(diagDeriv)
      allocate(diagDeriv(size(deriv), size(deriv)))
      diagDeriv(:,:) = 0.0_dp
      do jj = 1, size(this%layers(ii)%transferDeriv(aarg))
        diagDeriv(jj, jj) = deriv(jj)
      end do

      allocate(tmpMat(size(transpose(this%layers(ii - 1)%ww), dim=1), size(jacobi, dim=2)))
      tmpMat(:,:) = 0.0_dp
      call gemm(tmpMat, transpose(this%layers(ii - 1)%ww), jacobi)
      allocate(tmpMat2(size(diagDeriv, dim=1), size(tmpMat, dim=2)))
      tmpMat2(:,:) = 0.0_dp
      call gemm(tmpMat2, diagDeriv, tmpMat)
      jacobi = tmpMat2
      deallocate(tmpMat, tmpMat2)

    end do

  end function TNetwork_fdevi


  !> Determines the gradients in weight-bias space by a simple back-propagation algorithm.
  subroutine TNetwork_bprop(this, lossgrad, dw, db)

    !> representation of a neural network
    class(TNetwork), intent(inout) :: this

    !> loss gradient, with respect to predictions and targets
    real(dp), intent(in) :: lossgrad(:)

    !> weight gradients
    type(TWeightDerivs), intent(out) :: dw

    !> bias gradients
    type(TBiasDerivs), intent(out) :: db

    !! temporary matrix * vector storage
    real(dp), allocatable :: tmpVec(:)

    !! number of arrays in the bias structure
    integer :: nArrays

    !! auxiliary variable
    integer :: ii

    call TWeightDerivs_init(this%dims, dw)
    call TBiasDerivs_init(this%dims, db)

    nArrays = size(this%dims)

    db%db(nArrays)%array = lossgrad * this%layers(nArrays)%transferDeriv(this%layers(nArrays)%aarg)
    dw%dw(nArrays - 1)%array(:,:) = 0.0_dp
    call gemm(dw%dw(nArrays - 1)%array,&
        & reshape(this%layers(nArrays - 1)%aa, [this%dims(nArrays - 1), 1]),&
        & reshape(db%db(nArrays)%array, [1, this%dims(nArrays)]))

    do ii = size(this%dims) - 1, 2, -1

      allocate(tmpVec(size(this%layers(ii)%ww, dim=1)))
      tmpVec(:) = 0.0_dp
      call gemv(tmpVec, this%layers(ii)%ww, db%db(ii + 1)%array)
      db%db(ii)%array = tmpVec * this%layers(ii)%transferDeriv(this%layers(ii)%aarg)
      deallocate(tmpVec)

      dw%dw(ii - 1)%array(:,:) = 0.0_dp
      call gemm(dw%dw(ii - 1)%array, reshape(this%layers(ii - 1)%aa, [this%dims(ii - 1), 1]),&
          & reshape(db%db(ii)%array, [1, this%dims(ii)]))

    end do

  end subroutine TNetwork_bprop


  !> Extracts the activations of the last layer of the network, i.e. its output.
  subroutine TNetwork_getOutput(this, output)

    !> representation of a neural network
    class(TNetwork), intent(in) :: this

    !> output of current network
    real(dp), intent(out), allocatable :: output(:)

    output = this%layers(size(this%dims))%aa

  end subroutine TNetwork_getOutput


  !> Calculates the network output for a single set of input features.
  function TNetwork_iPredict(this, xx) result(aa)

    !> representation of a neural network
    class(TNetwork), intent(in) :: this

    !> single sample of input data to calculate output for
    real(dp), intent(in) :: xx(:)

    !! activation values
    real(dp), allocatable :: aa(:)

    !! temporary matrix * vector storage
    real(dp), allocatable :: tmpVec(:)

    !! auxiliary variable
    integer :: ii

    allocate(tmpVec(size(transpose(this%layers(1)%ww), dim=1)))
    tmpVec(:) = 0.0_dp
    call gemv(tmpVec, transpose(this%layers(1)%ww), xx)
    aa = this%layers(2)%transfer(tmpVec + this%layers(2)%bb)
    deallocate(tmpVec)

    do ii = 3, size(this%layers)

      allocate(tmpVec(size(transpose(this%layers(ii - 1)%ww), dim=1)))
      tmpVec(:) = 0.0_dp
      call gemv(tmpVec, transpose(this%layers(ii - 1)%ww), aa)
      aa = this%layers(ii)%transfer(tmpVec + this%layers(ii)%bb)
      deallocate(tmpVec)

    end do

  end function TNetwork_iPredict


  !> Calculates the network output for a batch of input features.
  function TNetwork_nPredict(this, xx) result(aa)

    !> representation of a neural network
    class(TNetwork), intent(in) :: this

    !> batch of input data to calculate output for
    real(dp), intent(in) :: xx(:,:)

    !! activation values
    real(dp), allocatable :: aa(:,:)

    !! auxiliary variable
    integer :: ii

    allocate(aa(this%dims(size(this%dims)), size(xx, dim=2)))

    do ii = 1, size(xx, dim=2)
      aa(:, ii) = this%iPredict(xx(:, ii))
    end do

  end function TNetwork_nPredict


  !> Sets the activation/transfer function for all layers of the network.
  subroutine TNetwork_setTransferFunc(this, descriptor)

    !> representation of a neural network
    class(TNetwork), intent(inout) :: this

    !> type of transfer function to use
    character(len=*), intent(in) :: descriptor

    !! auxiliary variable
    integer :: iLayer

    do iLayer = 1, size(this%dims)
      call this%layers(iLayer)%setTransferFunc(descriptor)
    end do

    !! to be able to predict arbitrary real numbers, always choose a linear output function
    call this%layers(size(this%dims))%setTransferFunc('linear')

  end subroutine TNetwork_setTransferFunc


  !> Maps the network parameters to a serialized array, suitable for the gradient-based optimizer.
  subroutine TNetwork_serializedWeightsAndBiases(this, weightsAndBiases)

    !> representation of a neural network
    class(TNetwork), intent(in) :: this

    !> serialized weights and biases
    real(dp), intent(out), allocatable :: weightsAndBiases(:)

    !! auxiliary variables
    integer :: ii, iLayer, ind

    allocate(weightsAndBiases(this%nBiases + this%nWeights))

    ind = 1

    ! serialize all weights
    do iLayer = 1, size(this%dims)
      do ii = 1, size(this%layers(iLayer)%ww, dim=2)
        weightsAndBiases(ind:ind + size(this%layers(iLayer)%ww, dim=1) - 1) =&
            & this%layers(iLayer)%ww(:, ii)
        ind = ind + size(this%layers(iLayer)%ww, dim=1)
      end do
    end do

    ! serialize all biases
    do iLayer = 1, size(this%dims)
      weightsAndBiases(ind:ind + size(this%layers(iLayer)%bb) - 1) = this%layers(iLayer)%bb
      ind = ind + size(this%layers(iLayer)%bb)
    end do

  end subroutine TNetwork_SerializedWeightsAndBiases


  !> Inserts serialized weights and biases into the network structure.
  subroutine TNetwork_serialWeightsAndBiasesFillup(this, weightsAndBiases)

    !> representation of a neural network
    class(TNetwork), intent(inout) :: this

    !> serialized weights and biases
    real(dp), intent(in) :: weightsAndBiases(:)

    !! auxiliary variables
    integer :: ii, iLayer, ind

    ind = 1

    ! fillup all weights
    do iLayer = 1, size(this%dims)
      do ii = 1, size(this%layers(iLayer)%ww, dim=2)
        this%layers(iLayer)%ww(:, ii) =&
            & weightsAndBiases(ind:ind + size(this%layers(iLayer)%ww, dim=1) - 1)
        ind = ind + size(this%layers(iLayer)%ww, dim=1)
      end do
    end do

    ! fillup all biases
    do iLayer = 1, size(this%dims)
      this%layers(iLayer)%bb(:) = weightsAndBiases(ind:ind + size(this%layers(iLayer)%bb) - 1)
      ind = ind + size(this%layers(iLayer)%bb)
    end do

  end subroutine TNetwork_serialWeightsAndBiasesFillup


  !> Resets the activations of all layers and all neurons of the network.
  subroutine TNetwork_resetActivations(this)

    !> representation of a neural network
    class(TNetwork), intent(inout) :: this

    !! auxiliary variable
    integer :: ii

    do ii = 1, size(this%layers)
      this%layers(ii)%aarg(:) = 0.0_dp
      this%layers(ii)%aa(:) = 0.0_dp
    end do

  end subroutine TNetwork_resetActivations

end module fnet_network
