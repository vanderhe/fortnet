!--------------------------------------------------------------------------------------------------!
!  FORTNET: A Behler-Parrinello-Neural-Network Implementation                                      !
!  Copyright (C) 2020 - 2021  T. W. van der Heide                                                  !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Defines several nested derived types used by other modules.
module fnet_nestedtypes

  use dftbp_accuracy, only: dp
  use dftbp_steepdesc, only : TSteepdesc
  use dftbp_conjgrad, only : TConjGrad
  use dftbp_lbfgs, only : TLbfgs
  use dftbp_fire, only : TFire

  use fnet_layer, only : TLayer, TLayer_init

#:if WITH_MPI
  use fnet_mpifx
#:endif

  implicit none

  private

  public :: TJacobian, TJacobians
  public :: TIntArray1D, TIntArray2D, TRealArray1D, TRealArray2D, TRealArray3D
  public :: TBiasDerivs, TBiasDerivs_init, TWeightDerivs, TWeightDerivs_init
  public :: TWrapSteepDesc, TWrapConjGrad, TWrapLbfgs, TWrapFire
  public :: TSingleLayerStruc, TMultiLayerStruc, TMultiLayerStruc_init
  public :: TDerivs, TDerivs_init
  public :: TPredicts, TPredicts_init
  public :: TEnv, TEnv_init


  type :: TIntArray1D

    !> onedimensional integer array
    integer, allocatable :: array(:)

  end type TIntArray1D


  type :: TIntArray2D

    !> twodimensional integer array
    integer, allocatable :: array(:)

  end type TIntArray2D


  type :: TRealArray1D

    !> onedimensional real array
    real(dp), allocatable :: array(:)

  end type TRealArray1D


  type :: TRealArray2D

    !> twodimensional real array
    real(dp), allocatable :: array(:,:)

  end type TRealArray2D


  type :: TRealArray3D

    !> threedimensional real array
    real(dp), allocatable :: array(:,:,:)

  end type TRealArray3D


  type :: TBiasDerivs

    !> contains derivatives of biases
    type(TRealArray1D), allocatable :: db(:)

  end type TBiasDerivs


  type :: TWeightDerivs

    !> contains derivatives of weights
    type(TRealArray2D), allocatable :: dw(:)

  contains

    procedure :: elasticNetRegularization => TWeightDerivs_elasticNetRegularization

  end type TWeightDerivs


  type :: TPredicts

    !> contains predictions of networks
    type(TRealArray2D), allocatable :: sys(:)

  end type TPredicts


  type :: TJacobian

    !> contains Jacobians of a single system
    type(TRealArray2D), allocatable :: atom(:)

  end type TJacobian


  type :: TJacobians

    !> contains Jacobians of multiple systems
    type(TJacobian), allocatable :: sys(:)

  end type TJacobians


  type :: TWrapSteepDesc

    !> steepest descent optimizer
    type(TSteepDesc), allocatable :: pSteepDesc

  end type TWrapSteepDesc


  type :: TWrapConjGrad

    !> conjugate gradients optimizer
    type(TConjGrad), allocatable :: pConjGrad

  end type TWrapConjGrad


  type :: TWrapLbfgs

    !> limited memory bfgs optimizer
    type(TLbfgs), allocatable :: pLbfgs

  end type TWrapLbfgs


  type :: TWrapFire

    !> fire optimizer
    type(TFire), allocatable :: pFire

  end type TWrapFire


  type :: TSingleLayerStruc

    !> representation of the network layers
    type(TLayer), allocatable :: layers(:)

  end type TSingleLayerStruc


  type :: TMultiLayerStruc

    !> representation of multiple network layers
    type(TSingleLayerStruc), allocatable :: struc(:)

  end type TMultiLayerStruc


  type :: TDerivs

    !> contains derivatives of biases
    type(TBiasDerivs), allocatable :: db(:)

    !> contains derivatives of weights
    type(TWeightDerivs), allocatable :: dw(:)

  contains

    procedure :: serialized => TDerivs_serialized
    procedure :: reset => TDerivs_reset

  end type TDerivs


  !> If compiled with mpi enabled, contains mpi communicator
  type TEnv

  #:if WITH_MPI
    !> mpi communicator with some additional information
    type(mpifx_comm) :: globalMpiComm
  #:endif

    !> true, if z-score standardization should be applied
    logical :: tWithMpi

  end type TEnv


contains

  !> Initialises a prediction structure based on a reference.
  pure subroutine TPredicts_init(this, reference)

    !> representation of neural network predictions
    type(TPredicts), intent(out) :: this

    !> reference array to get allocation structure from
    type(TRealArray2D), intent(in) :: reference(:)

    !> auxiliary variable
    integer :: iSys

    allocate(this%sys(size(reference)))

    do iSys = 1, size(reference)
      allocate(this%sys(iSys)%array(size(reference(iSys)%array, dim=1),&
          & size(reference(iSys)%array, dim=2)))
      this%sys(iSys)%array(:,:) = 0.0_dp
    end do

  end subroutine TPredicts_init


  !> Initialises an MPI environment.
  subroutine TEnv_init(this)

    !> mpi environment instance
    type(TEnv), intent(out) :: this

  #:if WITH_MPI
    ! initialise mpi environment
    call this%globalMpiComm%init()
    this%tWithMpi = .true.
  #:else
    this%tWithMpi = .false.
  #:endif

  end subroutine TEnv_init


  !> Initialises a structure that holds the derivatives of bias parameters.
  pure subroutine TBiasDerivs_init(dims, db)

    !> dimensions of all layers in the network
    integer, intent(in) :: dims(:)

    !> representation of bias derivatives
    type(TBiasDerivs), intent(out) :: db

    !> number of arrays in the bias structure
    integer :: nArrays

    !> auxiliary variable
    integer :: ii

    nArrays = size(dims)

    allocate(db%db(nArrays))

    do ii = 1, nArrays - 1
      allocate(db%db(ii)%array(dims(ii)))
      db%db(ii)%array(:) = 0.0_dp
    end do

    allocate(db%db(nArrays)%array(dims(nArrays)))
    db%db(nArrays)%array(:) = 0.0_dp

  end subroutine TBiasDerivs_init


  !> Initialises a structure that holds the derivatives of weight parameters.
  pure subroutine TWeightDerivs_init(dims, dw)

    !> dimensions of all layers in the network
    integer, intent(in) :: dims(:)

    !> weight derivatives
    type(TWeightDerivs), intent(out) :: dw

    !> number of arrays in the weight structure
    integer :: nArrays

    !> auxiliary variable
    integer :: ii

    nArrays = size(dims)

    allocate(dw%dw(nArrays))

    do ii = 1, nArrays - 1
      allocate(dw%dw(ii)%array(dims(ii), dims(ii + 1)))
      dw%dw(ii)%array(:,:) = 0.0_dp
    end do

    allocate(dw%dw(nArrays)%array(dims(ii), 1))
    dw%dw(nArrays)%array(:,:) = 0.0_dp

  end subroutine TWeightDerivs_init


  !> Applies an elastic net regularization ontop of existing gradients.
  pure subroutine TWeightDerivs_elasticNetRegularization(this, layers, nWeights, lambda, alpha)

    !> weight derivatives
    class(TWeightDerivs), intent(inout) :: this

    !> layers to extract weight structure from
    type(TLayer), intent(in) :: layers(:)

    !> number of weights to regularize
    integer, intent(in) :: nWeights

    !> strength of the regularization
    real(dp), intent(in) :: lambda

    !> elastic net weighting between lasso and ridge (0 = ridge, 1 = lasso)
    real(dp), intent(in) :: alpha

    !> temporary array to build up signum function values
    real(dp), allocatable :: sgn(:,:)

    !> auxiliary variable
    integer :: iLayer

    do iLayer = 1, size(this%dw)
      sgn = sign(1.0_dp, layers(iLayer)%ww)
      where (layers(iLayer)%ww == 0.0_dp)
        sgn(:,:) = 0.0_dp
      elsewhere
        sgn(:,:) = sgn
      end where
      this%dw(iLayer)%array(:,:) = this%dw(iLayer)%array + lambda / real(nWeights, dp)&
          & * ((1.0_dp - alpha) * layers(iLayer)%ww + alpha * sgn)
    end do

  end subroutine TWeightDerivs_elasticNetRegularization


  !> Initialises a structure that holds multiple network layer structures.
  subroutine TMultiLayerStruc_init(dims, nStrucs, this)

    !> dimensions of all layers in the sub-nn's
    integer, intent(in) :: dims(:)

    !> number of layer structures
    integer, intent(in) :: nStrucs

    !> representation of temporary layer storage container
    type(TMultiLayerStruc), intent(out) :: this

    !> auxiliary variables
    integer :: ii, jj

    allocate(this%struc(nStrucs))

    do ii = 1, nStrucs
      allocate(this%struc(ii)%layers(size(dims)))
      do jj = 1, size(dims) - 1
        call TLayer_init(this%struc(ii)%layers(jj), dims(jj), dims(jj + 1))
      end do
      call TLayer_init(this%struc(ii)%layers(jj), dims(jj), 1)
    end do

  end subroutine TMultiLayerStruc_init


  !> Initialises a wrapper around bias and weight derivatives.
  subroutine TDerivs_init(dims, nStrucs, derivs)

    !> dimensions of all layers in the sub-nn's
    integer, intent(in) :: dims(:)

    !> number of sub-structures to initialise
    integer, intent(in) :: nStrucs

    !> representation weight and bias derivatives
    class(TDerivs), intent(out) :: derivs

    !> auxiliary variables
    integer :: iStruc

    allocate(derivs%dw(nStrucs))
    allocate(derivs%db(nStrucs))

    do iStruc = 1, nStrucs
      call TWeightDerivs_init(dims, derivs%dw(iStruc))
      call TBiasDerivs_init(dims, derivs%db(iStruc))
    end do

  end subroutine TDerivs_init


  !> Serialises the gradients of a derivative type.
  pure subroutine TDerivs_serialized(this, nTotParams, ddSerial)

    !> representation of weight and bias derivatives
    class(TDerivs), intent(in) :: this

    !> number of total parameters in each sub-nn
    integer, intent(in) :: nTotParams

    !> all gradients collected in 1d-array, shape: [nTotGrads, nSpecies]
    real(dp), intent(out), allocatable :: ddSerial(:,:)

    !> auxiliary variables
    integer :: iStruc, iLayer, nStrucs, nArrays, ii, jj, ind

    nStrucs = size(this%db)
    nArrays = size(this%db(1)%db)
    allocate(ddSerial(nTotParams, nStrucs))

    do iStruc = 1, nStrucs

      ind = 1

      ! serialize all weight gradients
      do iLayer = 1, nArrays
        do ii = 1, size(this%dw(iStruc)%dw(iLayer)%array, dim=2)
          do jj = 1, size(this%dw(iStruc)%dw(iLayer)%array, dim=1)
            ddSerial(ind, iStruc) = this%dw(iStruc)%dw(iLayer)%array(jj, ii)
            ind = ind + 1
          end do
        end do
      end do

      ! serialize all bias gradients
      do iLayer = 1, nArrays
        do ii = 1, size(this%db(iStruc)%db(iLayer)%array)
          ddSerial(ind, iStruc) = this%db(iStruc)%db(iLayer)%array(ii)
          ind = ind + 1
        end do
      end do

    end do

  end subroutine TDerivs_serialized


  !> Resets a derivative structure to hold zero-valued entries.
  subroutine TDerivs_reset(this)

    !> representation of weight and bias derivatives
    class(TDerivs), intent(inout) :: this

    !> auxiliary variables
    integer :: iStruc, iArray, nStrucs, nArrays

    nStrucs = size(this%db)
    nArrays = size(this%db(1)%db)

    do iStruc = 1, nStrucs
      do iArray = 1, nArrays - 1
        this%dw(iStruc)%dw(iArray)%array(:,:) = 0.0_dp
      end do
      this%dw(iStruc)%dw(nArrays)%array(:,:) = 0.0_dp
    end do

    do iStruc = 1, nStrucs
      do iArray = 1, nArrays - 1
        this%db(iStruc)%db(iArray)%array(:) = 0.0_dp
      end do
      this%db(iStruc)%db(nArrays)%array(:) = 0.0_dp
    end do

  end subroutine TDerivs_reset

end module fnet_nestedtypes
