!--------------------------------------------------------------------------------------------------!
!  FORTNET: A Behler-Parrinello-Neural-Network Implementation                                      !
!  Copyright (C) 2020 - 2022  T. W. van der Heide                                                  !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Defines a layer of a dense feed-forward multilayer-perceptron.
module fnet_layer

  use dftbp_accuracy, only: dp
  use dftbp_ranlux, only : TRanlux

  use fnet_transfer, only : transferFunc, gaussian, gaussianDeriv, relu, reluDeriv, lrelu,&
      & lreluDeriv, softPlus, softPlusDeriv, bentIdentity, bentIdentityDeriv, arctan, arctanDeriv,&
      & sigmoid, sigmoidDeriv, heaviside, heavisideDeriv, tanhf, tanhDeriv, linear, linearDeriv
  use fnet_random, only : normalXavier

  implicit none

  private

  public :: TLayer, TLayer_init


  type :: TLayer

    !> activation status of all neurons
    real(dp), allocatable :: aa(:)

    !> applied biases of all neurons
    real(dp), allocatable :: bb(:)

    !> weights of all neurons
    real(dp), allocatable :: ww(:,:)

    !> input to calculate transfer function, gets storen during fprop
    real(dp), allocatable :: aarg(:)

    !> descriptor of transfer function to use
    character(len=:), allocatable :: descriptor

    procedure(transferFunc), pointer, nopass :: transfer => null()
    procedure(transferFunc), pointer, nopass :: transferDeriv => null()

  contains

    procedure :: setTransferFunc => TLayer_setTransferFunc

  end type TLayer


contains

  !> Initialises a single layer instance (+ initial activations if desired).
  subroutine TLayer_init(this, nCurrentNeurons, nNextNeurons, rndGen)

    !> representation of a neural network layer
    type(TLayer), intent(inout) :: this

    !> number of neurons in the adjacent layer
    integer, intent(in) :: nNextNeurons

    !> number of neurons in the current layer
    integer, intent(in) :: nCurrentNeurons

    !> luxury pseudorandom generator instance
    type(TRanlux), intent(inout), optional :: rndGen

    allocate(this%aa(nCurrentNeurons))
    allocate(this%aarg(nCurrentNeurons))
    allocate(this%ww(nCurrentNeurons, nNextNeurons))
    allocate(this%bb(nCurrentNeurons))

    this%aa(:) = 0.0_dp
    this%aarg(:) = 0.0_dp

    if (present(rndGen)) then
      call normalXavier(rndGen, this%ww, nNextNeurons, nCurrentNeurons)
    else
      this%ww(:,:) = 0.0_dp
    end if

    this%bb(:) = 0.0_dp

  end subroutine TLayer_init


  !> Sets the activation/transfer function as well as its derivative.
  subroutine TLayer_setTransferFunc(this, descriptor)

    !> representation of a neural network layer
    class(TLayer), intent(inout) :: this

    !> descriptor of transfer function to use
    character(len=*), intent(in) :: descriptor

    select case(trim(descriptor))

      case('gaussian')
        this%transfer => gaussian
        this%transferDeriv => gaussianDeriv
        this%descriptor = 'gaussian'

      case('relu')
        this%transfer => relu
        this%transferDeriv => reluDeriv
        this%descriptor = 'relu'

      case('lrelu')
        this%transfer => lrelu
        this%transferDeriv => lreluDeriv
        this%descriptor = 'lrelu'

      case('softplus')
        this%transfer => softPlus
        this%transferDeriv => softPlusDeriv
        this%descriptor = 'softplus'

      case('bent')
        this%transfer => bentIdentity
        this%transferDeriv => bentIdentityDeriv
        this%descriptor = 'bent'

      case('atan')
        this%transfer => arctan
        this%transferDeriv => arctanDeriv
        this%descriptor = 'atan'

      case('sigmoid')
        this%transfer => sigmoid
        this%transferDeriv => sigmoidDeriv
        this%descriptor = 'sigmoid'

      case('heaviside')
        this%transfer => heaviside
        this%transferDeriv => heavisideDeriv
        this%descriptor = 'heaviside'

      case('tanh')
        this%transfer => tanhf
        this%transferDeriv => tanhDeriv
        this%descriptor = 'tanh'

      case('linear')
        this%transfer => linear
        this%transferDeriv => linearDeriv
        this%descriptor = 'linear'

      case default
        this%transfer => tanhf
        this%transferDeriv => tanhDeriv
        this%descriptor = 'tanh'

    end select

  end subroutine TLayer_setTransferFunc

end module fnet_layer
