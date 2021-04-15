!--------------------------------------------------------------------------------------------------!
!  FORTNET: A Behler-Parrinello-Neural-Network Implementation                                      !
!  Copyright (C) 2020 - 2021  T. W. van der Heide                                                  !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

module fnet_network

  use dftbp_assert
  use dftbp_message, only : error
  use dftbp_charmanip, only : tolower
  use dftbp_accuracy, only: dp
  use dftbp_ranlux, only : TRanlux

  use fnet_initprogram
  use fnet_nestedtypes, only : TRealArray1D, TRealArray2D
  use fnet_nestedtypes, only : TBiasDerivs, TBiasDerivs_init
  use fnet_nestedtypes, only : TWeightDerivs, TWeightDerivs_init
  use fnet_layer, only : TLayer, TLayer_init
  use fnet_loss, only : deviation

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

  contains

    procedure :: initLayers => TNetwork_initLayers
    procedure :: countParams => TNetwork_countParams
    procedure :: fprop => TNetwork_fprop
    procedure :: bprop => TNetwork_bprop
    procedure :: getOutput => TNetwork_getOutput
    procedure :: iPredict => TNetwork_iPredict
    procedure :: nPredict => TNetwork_nPredict
    procedure :: resetActivations => TNetwork_resetActivations
    procedure :: setTransferFunc => TNetwork_setTransferFunc
    procedure :: serializedWeightsAndBiases => TNetwork_serializedWeightsAndBiases
    procedure :: serialWeightsAndBiasesFillup => TNetwork_serialWeightsAndBiasesFillup
    procedure :: toFile => TNetwork_toFile
    procedure :: fromFile => TNetwork_fromFile

  end type TNetwork

contains

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
    else
      call this%setTransferFunc('tanh')
    end if

    call this%countParams()

  end subroutine TNetwork_init


  subroutine TNetwork_initLayers(this, dims, rndGen)

    !> representation of a neural network
    class(TNetwork), intent(inout) :: this

    !> dimensions of all layers in the network
    integer, intent(in) :: dims(:)

    !> luxury pseudorandom generator instance
    type(TRanlux), intent(inout), optional :: rndGen

    !> auxiliary variable
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


  subroutine TNetwork_countParams(this)

    !> representation of a neural network
    class(TNetwork), intent(inout) :: this

    !> auxiliary variable
    integer :: iLayer

    this%nBiases = sum(this%dims)

    this%nWeights = 0
    do iLayer = 1, size(this%dims) - 1
      this%nWeights = this%nWeights + this%dims(iLayer) * this%dims(iLayer + 1)
    end do
    this%nWeights = this%nWeights + this%dims(size(this%dims))

  end subroutine TNetwork_countParams


  subroutine TNetwork_fprop(this, xx, state, out)

    !> representation of a neural network
    class(TNetwork), intent(inout) :: this

    !> input to feed network with
    real(dp), intent(in) :: xx(:)

    !> copy of representation of the network layers
    type(TLayer), intent(out), allocatable, optional :: state(:)

    !> outputs of current neural network instance
    real(dp), intent(out), allocatable, optional :: out(:)

    !> auxiliary variable
    integer :: ii

    this%layers(1)%aa = xx

    do ii = 2, size(this%layers)
      this%layers(ii)%aarg = matmul(transpose(this%layers(ii - 1)%ww), this%layers(ii - 1)%aa) +&
          & this%layers(ii)%bb
      this%layers(ii)%aa = this%layers(ii)%transfer(this%layers(ii)%aarg)
    end do

    state = this%layers

    call this%getOutput(out)

  end subroutine TNetwork_fprop


  subroutine TNetwork_bprop(this, loss, dw, db)

    !> representation of a neural network
    class(TNetwork), intent(inout) :: this

    !> loss, by comparison of predictions and targets
    real(dp), intent(in) :: loss(:)

    !> weight gradients
    type(TWeightDerivs), intent(out) :: dw

    !> bias gradients
    type(TBiasDerivs), intent(out) :: db

    !> number of arrays in the bias structure
    integer :: nArrays

    !> auxiliary variable
    integer :: ii

    call TWeightDerivs_init(this%dims, dw)
    call TBiasDerivs_init(this%dims, db)

    nArrays = size(this%dims)

    db%db(nArrays)%array = loss * this%layers(nArrays)%transferDeriv(this%layers(nArrays)%aarg)
    dw%dw(nArrays - 1)%array = matmul(reshape(this%layers(nArrays - 1)%aa,&
        & [this%dims(nArrays - 1), 1]), reshape(db%db(nArrays)%array, [1, this%dims(nArrays)]))

    do ii = size(this%dims) - 1, 2, -1
      db%db(ii)%array = matmul(this%layers(ii)%ww, db%db(ii + 1)%array) *&
          & this%layers(ii)%transferDeriv(this%layers(ii)%aarg)
      dw%dw(ii - 1)%array = matmul(reshape(this%layers(ii - 1)%aa, [this%dims(ii - 1), 1]),&
          & reshape(db%db(ii)%array, [1, this%dims(ii)]))
    end do

  end subroutine TNetwork_bprop


  subroutine TNetwork_getOutput(this, output)

    !> representation of a neural network
    class(TNetwork), intent(in) :: this

    !> output of current network
    real(dp), intent(out), allocatable :: output(:)

    output = this%layers(size(this%dims))%aa

  end subroutine TNetwork_getOutput


  function TNetwork_iPredict(this, xx) result(aa)

    !> representation of a neural network
    class(TNetwork), intent(in) :: this

    !> single sample of input data to calculate output for
    real(dp), intent(in) :: xx(:)

    !> activation values
    real(dp), allocatable :: aa(:)

    !> auxiliary variable
    integer :: ii

    aa = this%layers(2)%transfer(matmul(transpose(this%layers(1)%ww), xx) + this%layers(2)%bb)

    do ii = 3, size(this%layers)
      aa = this%layers(ii)%transfer(matmul(transpose(this%layers(ii - 1)%ww), aa) +&
          & this%layers(ii)%bb)
    end do

  end function TNetwork_iPredict


  function TNetwork_nPredict(this, xx) result(aa)

    !> representation of a neural network
    class(TNetwork), intent(in) :: this

    !> batch of input data to calculate output for
    real(dp), intent(in) :: xx(:,:)

    !> activation values
    real(dp), allocatable :: aa(:,:)

    !> auxiliary variable
    integer :: ii

    allocate(aa(this%dims(size(this%dims)), size(xx, dim=2)))

    do ii = 1, size(xx, dim=2)
      aa(:, ii) = this%iPredict(xx(:, ii))
    end do

  end function TNetwork_nPredict


  subroutine TNetwork_setTransferFunc(this, descriptor)

    !> representation of a neural network
    class(TNetwork), intent(inout) :: this

    !> type of transfer function to use
    character(len=*), intent(in) :: descriptor

    !> Layer identifier
    integer :: iLayer

    do iLayer = 1, size(this%dims)
      call this%layers(iLayer)%setTransferFunc(descriptor)
    end do

    !! to be able to predict arbitrary real numbers, always choose a linear output function
    call this%layers(size(this%dims))%setTransferFunc('linear')

  end subroutine TNetwork_setTransferFunc


  subroutine TNetwork_serializedWeightsAndBiases(this, weightsAndBiases)

    !> representation of a neural network
    class(TNetwork), intent(in) :: this

    !> serialized weights and biases
    real(dp), intent(out), allocatable :: weightsAndBiases(:)

    !> Layer identifier
    integer :: iLayer

    !> auxiliary variable
    integer :: ii, ind

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
      weightsAndBiases(ind:ind + size(this%layers(iLayer)%bb) - 1) = this%layers(iLayer)%bb(:)
      ind = ind + size(this%layers(iLayer)%bb)
    end do

  end subroutine TNetwork_SerializedWeightsAndBiases


  subroutine TNetwork_serialWeightsAndBiasesFillup(this, weightsAndBiases)

    !> representation of a neural network
    class(TNetwork), intent(inout) :: this

    !> serialized weights and biases
    real(dp), intent(in) :: weightsAndBiases(:)

    !> Layer identifier
    integer :: iLayer

    !> auxiliary variable
    integer :: ii, ind

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


  subroutine TNetwork_resetActivations(this)

    !> representation of a neural network
    class(TNetwork), intent(inout) :: this

    !> auxiliary variable
    integer :: ii

    do ii = 1, size(this%layers)
      this%layers(ii)%aarg(:) = 0.0_dp
      this%layers(ii)%aa(:) = 0.0_dp
    end do

  end subroutine TNetwork_resetActivations


  ! subroutine TNetwork_update(this, dw, db, eta)

  !   !> representation of a neural network
  !   class(TNetwork), intent(inout) :: this

  !   !> weight gradients
  !   type(TRealArray2D), intent(in) :: dw(:)

  !   !> bias gradients
  !   type(TRealArray1D), intent(in) :: db(:)

  !   !> learning rate
  !   real(dp), intent(in) :: eta

  !   !> auxiliary variables
  !   integer :: ii, jj

  !   ! update biases
  !   do ii = 2, size(this%dims)
  !     this%layers(ii)%bb = this%layers(ii)%bb - eta * db(ii)%array
  !   end do

  !   ! update weights
  !   do jj = 1, size(this%dims) - 1
  !     this%layers(jj)%ww = this%layers(jj)%ww - eta * dw(jj)%array
  !   end do

  !   end subroutine TNetwork_update


  subroutine TNetwork_toFile(this, prog, iSpecies)

    !> representation of a neural network
    class(TNetwork), intent(in) :: this

    !> representation of program variables
    type(TProgramVariables), intent(in) :: prog

    !> global species index of sub-nn
    integer, intent(in) :: iSpecies

    !> unique fileunit
    integer :: fd

    !> auxiliary variables
    integer :: ii, jj

    open(newunit=fd, file=prog%data%netstatNames(iSpecies), form='formatted', status='replace',&
        & action='write')

    if (prog%data%tAtomicTargets) then
      write(fd, '(2A10)') prog%arch%type, 'atomic'
    else
      write(fd, '(2A10)') prog%arch%type, 'global'
    end if

    write(fd, '(A)') prog%data%globalSpNames(iSpecies)

    write(fd, *) size(this%dims)
    write(fd, *) this%dims

    write(fd, '(A)') prog%arch%activation

    do ii = 2, size(this%dims)
      write(fd, '(ES26.16E3)') this%layers(ii)%bb
    end do
    do jj = 1, size(this%dims) - 1
      write(fd, '(ES26.16E3)') this%layers(jj)%ww
    end do

    close(fd)

  end subroutine TNetwork_toFile


  subroutine TNetwork_fromFile(this, prog, iSpecies)

    !> representation of a neural network
    class(TNetwork), intent(out) :: this

    !> representation of program variables
    type(TProgramVariables), intent(inout) :: prog

    !> global species index of sub-nn
    integer, intent(in) :: iSpecies

    !> unique fileunit
    integer :: fd

    !> auxiliary variables
    integer :: ii, jj, nLayers
    character(len=50) :: spName, archType, targetType

    open(newunit=fd, file=prog%data%netstatNames(iSpecies), form='formatted', status='old',&
        & action='read')

    read(fd, *) archType, targetType
    prog%arch%type = tolower(trim(archType))
    if (tolower(trim(targetType)) == 'atomic') then
      prog%data%tAtomicTargets = .true.
    elseif (tolower(trim(targetType)) == 'global') then
      prog%data%tAtomicTargets = .false.
    else
      call error("Unrecognized target type in file '" // prog%data%netstatNames(iSpecies) // "'.")
    end if

    read(fd, *) spName
    @:ASSERT(tolower(prog%data%globalSpNames(iSpecies)) == tolower(spName))

    read(fd, *) nLayers
    @:ASSERT(nLayers > 2)
    if (allocated(prog%arch%allDims)) deallocate(prog%arch%allDims)
    allocate(prog%arch%allDims(nLayers))
    read(fd, *) prog%arch%allDims

    prog%arch%hidden = prog%arch%allDims(2:size(prog%arch%allDims) - 1)
    prog%arch%nHiddenLayer = size(prog%arch%hidden)

    read(fd, *) prog%arch%activation

    call TNetwork_init(this, prog%arch%allDims, descriptor=prog%arch%activation)

    do ii = 2, size(this%dims)
      read(fd, *) this%layers(ii)%bb
    end do
    do jj = 1, size(this%dims) - 1
      read(fd, *) this%layers(jj)%ww
    end do

    close(fd)

  end subroutine TNetwork_fromFile

end module fnet_network
