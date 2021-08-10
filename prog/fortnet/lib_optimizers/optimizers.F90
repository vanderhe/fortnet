!--------------------------------------------------------------------------------------------------!
!  FORTNET: A Behler-Parrinello-Neural-Network Implementation                                      !
!  Copyright (C) 2020 - 2021  T. W. van der Heide                                                  !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Defines a general interface to the optimization algorithms.
module fnet_optimizers

  use dftbp_accuracy, only : dp

  use dftbp_steepdesc, only : TSteepDesc, init, reset, next
  use dftbp_conjgrad, only : TConjGrad, init, reset, next
  use dftbp_lbfgs, only : TLbfgs, TLbfgs_init
  use dftbp_fire, only : TFire, TFire_init

  implicit none
  private


  public :: TOptimizer
  public :: init, reset, next
  public :: optimizerTypes

  !> Interface type for the various optimization algorithms.
  type TOptimizer

    private

    !> integer ID of specified optimizer
    integer :: iOptimizer

    !> steepest descent optimizer
    type(TSteepDesc), allocatable :: pSteepDesc

    !> conjugate gradients optimizer
    type(TConjGrad), allocatable :: pConjGrad

    !> limited memory bfgs optimizer
    type(TLbfgs), allocatable :: pLbfgs

    !> fire optimizer
    type(TFire), allocatable :: pFire

  end type TOptimizer


  !> Initialises an optimizer instance.
  interface init

    module procedure TOptimizer_initTSteepDesc
    module procedure TOptimizer_initTConjGrad
    module procedure TOptimizer_initTLbfgs
    module procedure TOptimizer_initTFire

  end interface


  !> Resets the optimizer.
  interface reset
    module procedure TOptimizer_reset
  end interface


  !> Delivers the next point in the minimization.
  interface next
    module procedure TOptimizer_next
  end interface


  !> Assigns an integer ID to each optimizer.
  type :: TOptimizerTypesEnum

    integer :: none = 0
    integer :: steepDesc = 1
    integer :: conjGrad = 2
    integer :: lbfgs = 3
    integer :: fire = 4

  end type TOptimizerTypesEnum

  type(TOptimizerTypesEnum), parameter :: optimizerTypes = TOptimizerTypesEnum()


contains

  !> Creates a general optimizer with a steepest descent instance.
  subroutine TOptimizer_initTSteepDesc(pSteepDesc, this)

    !> an already initialized steepest descent instance
    type(TSteepDesc), intent(inout), allocatable :: pSteepDesc

    !> general optimizer instance
    type(TOptimizer), intent(out) :: this

    this%iOptimizer = optimizerTypes%steepDesc
    call move_alloc(pSteepDesc, this%pSteepDesc)

  end subroutine TOptimizer_initTSteepDesc


  !> Creates a general optimizer with a conjugate gradient instance.
  subroutine TOptimizer_initTConjGrad(pConjGrad, this)

    !> an already initialized conjugate gradient instance
    type(TConjGrad), intent(inout), allocatable :: pConjGrad

    !> general optimizer instance
    type(TOptimizer), intent(out) :: this

    this%iOptimizer = optimizerTypes%conjGrad
    call move_alloc(pConjGrad, this%pConjGrad)

  end subroutine TOptimizer_initTConjGrad


  !> Creates a general optimizer with a limited memory bfgs instance.
  subroutine TOptimizer_initTLbfgs(pLbfgs, this)

    !> an already initialized lbfgs instance
    type(TLbfgs), intent(inout), allocatable :: pLbfgs

    !> general optimizer instance
    type(TOptimizer), intent(out) :: this

    this%iOptimizer = optimizerTypes%lbfgs
    call move_alloc(pLbfgs, this%pLbfgs)

  end subroutine TOptimizer_initTLbfgs


  !> Creates a general optimizer with a fire instance.
  subroutine TOptimizer_initTFire(pFire, this)

    !> an already initialized fire instance
    type(TFire), intent(inout), allocatable :: pFire

    !> general optimizer instance
    type(TOptimizer), intent(out) :: this

    this%iOptimizer = optimizerTypes%fire
    call move_alloc(pFire, this%pFire)

  end subroutine TOptimizer_initTFire


  !> Resets a general optimizer.
  subroutine TOptimizer_reset(this, initialVals)

    !> general optimizer instance
    type(TOptimizer), intent(inout) :: this

    !> initial values to optimize
    real(dp), intent(in) :: initialVals(:)

    select case (this%iOptimizer)

    case(optimizerTypes%conjGrad)

      call reset(this%pConjGrad, initialVals)

    case(optimizerTypes%steepDesc)

      call reset(this%pSteepDesc, initialVals)

    case (optimizerTypes%lbfgs)

      call this%pLbfgs%reset(initialVals)

    case (optimizerTypes%fire)

      call this%pFire%reset(initialVals)

    end select

  end subroutine TOptimizer_reset


  !> Delivers the next optimized values.
  !> Function value and gradients of the initial values must get passed at first invocation.
  subroutine TOptimizer_next(this, funcVal, grads, outputVals, tConverged)

    !> general optimizer instance
    type(TOptimizer), intent(inout) :: this

    !> function value for last step returned by this routine
    real(dp), intent(in) :: funcVal

    !> Gradient in the last point
    real(dp), intent(in) :: grads(:)

    !> updated values of the optimization
    real(dp), intent(out) :: outputVals(:)

    !> true, if gradient got below the specified tolerance
    logical, intent(out) :: tConverged

    select case (this%iOptimizer)

    case(optimizerTypes%conjGrad)

      call next(this%pConjGrad, funcVal, grads, outputVals, tConverged)

    case (optimizerTypes%steepDesc)

      call next(this%pSteepDesc, grads, outputVals, tConverged)

    case (optimizerTypes%lbfgs)

      call this%pLbfgs%next(funcVal, grads, outputVals, tConverged)

    case (optimizerTypes%fire)

      call this%pFire%next(grads, outputVals, tConverged)

    end select

  end subroutine TOptimizer_next

end module fnet_optimizers
