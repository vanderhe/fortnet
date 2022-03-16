!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2021  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

#! suffix and kinds for real types
#:set REAL_KIND_PARAMS = [('real', 's'), ('dble', 'd')]

!> Contains F90 wrapper functions for some commonly used blas calls needed in the code. The
!> interface of all BLAS calls must be defined in the module blas.
module dftbp_math_blasroutines

  use dftbp_assert
  use dftbp_accuracy, only : dp, rsp, rdp

  implicit none
  private

  public ::  gemv, gemm


  !> General matrix vector multiply y := alpha*a*x + beta*y
  !> Wrapper for the level 2 blas routine
  interface gemv
    module procedure gemv_real
    module procedure gemv_dble
    #:for suffix, _ in REAL_KIND_PARAMS
      module procedure gemv231_${suffix}$
    #:endfor
    #:for suffix, _ in REAL_KIND_PARAMS
      module procedure gemv242_${suffix}$
    #:endfor
  end interface gemv


  !> Interface to GEMM routines evaluates C := alpha*op( A )*op( B ) + beta*C, where op( X ) is one
  !> of op( X ) = X or op( X ) = X'

  !> Wrapper for the level 3 blas routine
  interface gemm
    module procedure gemm_real
    module procedure gemm_dble
    #:for SUFFIX, _ in REAL_KIND_PARAMS
      module procedure gemm332_${SUFFIX}$
    #:endfor
  end interface gemm


contains

  !> real matrix*vector product
  subroutine gemv_real(y,a,x,alpha,beta,trans)

    !> vector
    real(rsp), intent(inout) :: y(:)

    !> matrix
    real(rsp), intent(in) :: a(:,:)

    !> vector
    real(rsp), intent(in) :: x(:)

    !> optional scaling factor (defaults to 1)
    real(rsp), intent(in), optional :: alpha

    !> optional scaling factor (defaults to 0)
    real(rsp), intent(in), optional :: beta

    !> optional transpose (defaults to 'n'), allowed choices are 'n', 'N', 't', 'T', 'c' and 'C'
    character, intent(in), optional :: trans

    integer :: n, m
    character :: iTrans
    real(rsp) :: iAlpha, iBeta

    if (present(trans)) then
      iTrans = trans
    else
      iTrans = 'n'
    end if
    if (present(alpha)) then
      iAlpha = alpha
    else
      iAlpha = 1.0_rsp
    end if
    if (present(beta)) then
      iBeta = beta
    else
      iBeta = 0.0_rsp
    end if

    @:ASSERT(iTrans == 'n' .or. iTrans == 'N' .or. iTrans == 't' .or. iTrans == 'T' .or.&
        & iTrans == 'c' .or. iTrans == 'C')
    @:ASSERT(((size(a,dim=1) == size(y)) .and. (iTrans == 'n' .or. iTrans == 'N')) .or.&
        & (size(a,dim=1) == size(x)))
    @:ASSERT(((size(a,dim=2) == size(x)) .and. (iTrans == 'n' .or. iTrans == 'N')) .or.&
        & (size(a,dim=2) == size(y)))

    m = size(a,dim=1)
    n = size(a,dim=2)

    call sgemv( iTrans, m, n, iAlpha, a, m, x, 1, iBeta, y, 1 )

  end subroutine gemv_real


  !> double precision matrix*vector product
  subroutine gemv_dble(y,a,x,alpha,beta,trans)

    !> vector
    real(rdp), intent(inout) :: y(:)

    !> matrix
    real(rdp), intent(in) :: a(:,:)

    !> vector
    real(rdp), intent(in) :: x(:)

    !> optional scaling factor (defaults to 1)
    real(rdp), intent(in), optional :: alpha

    !> optional scaling factor (defaults to 0)
    real(rdp), intent(in), optional :: beta

    !> optional transpose (defaults to 'n'), allowed choices are 'n', 'N', 't', 'T', 'c' and 'C'
    character, intent(in), optional :: trans

    integer :: n, m
    character :: iTrans
    real(rdp) :: iAlpha, iBeta

    if (present(trans)) then
      iTrans = trans
    else
      iTrans = 'n'
    end if
    if (present(alpha)) then
      iAlpha = alpha
    else
      iAlpha = 1.0_rdp
    end if
    if (present(beta)) then
      iBeta = beta
    else
      iBeta = 0.0_rdp
    end if

    @:ASSERT(iTrans == 'n' .or. iTrans == 'N' .or. iTrans == 't' .or. iTrans == 'T' .or.&
        & iTrans == 'c' .or. iTrans == 'C')
    @:ASSERT(((size(a,dim=1) == size(y)) .and. (iTrans == 'n' .or. iTrans == 'N')) .or.&
        & (size(a,dim=1) == size(x)))
    @:ASSERT(((size(a,dim=2) == size(x)) .and. (iTrans == 'n' .or. iTrans == 'N')) .or.&
        & (size(a,dim=2) == size(y)))

    m = size(a,dim=1)
    n = size(a,dim=2)

    call dgemv( iTrans, m, n, iAlpha, a, m, x, 1, iBeta, y, 1 )

  end subroutine gemv_dble


  #:for suffix, kind in REAL_KIND_PARAMS

    !> Generalized matrix vector contraction Cij = Aijk * Bk
    subroutine gemv231_${suffix}$(y, a, x, alpha, beta, trans)

      !> matrix
      real(r${kind}$p), intent(inout), contiguous, target :: y(:,:)

      !> matrix
      real(r${kind}$p), intent(in), contiguous, target :: a(:,:,:)

      !> vector
      real(r${kind}$p), intent(in) :: x(:)

      !> optional scaling factor (defaults to 1)
      real(r${kind}$p), intent(in), optional :: alpha

      !> optional scaling factor (defaults to 0)
      real(r${kind}$p), intent(in), optional :: beta

      !> optional transpose (defaults to 'n'), allowed choices are 'n', 'N', 't', 'T', 'c' and 'C'
      character, intent(in), optional :: trans

      real(r${kind}$p), pointer :: pY(:)
      real(r${kind}$p), pointer :: pA(:,:)

      pY(1 : size(y)) => y
      pA(1 : size(a, dim=1) * size(a, dim=2), 1 : size(a, dim=3)) => a
      call gemv(pY, pA, x, alpha, beta, trans)

    end subroutine gemv231_${suffix}$

  #:endfor

  #:for suffix, kind in REAL_KIND_PARAMS

    !> Generalized matrix vector contraction Cij = Aijk * Bk
    subroutine gemv242_${suffix}$(y, a, x, alpha, beta, trans)

      !> matrix
      real(r${kind}$p), intent(inout), contiguous, target :: y(:,:)

      !> matrix
      real(r${kind}$p), intent(in), contiguous, target :: a(:,:,:,:)

      !> matrix
      real(r${kind}$p), intent(in), contiguous, target :: x(:,:)

      !> optional scaling factor (defaults to 1)
      real(r${kind}$p), intent(in), optional :: alpha

      !> optional scaling factor (defaults to 0)
      real(r${kind}$p), intent(in), optional :: beta

      !> optional transpose (defaults to 'n'), allowed choices are 'n', 'N', 't', 'T', 'c' and 'C'
      character, intent(in), optional :: trans

      real(r${kind}$p), pointer :: pX(:)
      real(r${kind}$p), pointer :: pY(:)
      real(r${kind}$p), pointer :: pA(:,:)

      pY(1:size(y)) => y
      pA(1:size(a, dim=1)*size(a, dim=2), 1:size(a, dim=3)*size(a, dim=4)) => a
      pX(1:size(x)) => x
      call gemv(pY, pA, pX, alpha, beta, trans)

    end subroutine gemv242_${suffix}$

  #:endfor


#:for suffix, kind in REAL_KIND_PARAMS

  !> ${suffix}$ matrix*matrix product
  subroutine gemm_${suffix}$(C, A, B, alpha, beta, transA, transB, n, m, k, lda, ldb, ldc)

    !> general matrix output
    real(r${kind}$p), intent(inout) :: C(:,:)

    !> general matrix
    real(r${kind}$p), intent(in) :: A(:,:)

    !> general matrix
    real(r${kind}$p), intent(in) :: B(:,:)

    !> defaults to 1 if not set
    real(r${kind}$p), intent(in), optional :: alpha

    !> defaults to 0 if not set
    real(r${kind}$p), intent(in), optional :: beta

    !> optional transpose of A matrix (defaults to 'n'), allowed choices are 'n', 'N', 't', 'T', 'c'
    !> and 'C'
    character, intent(in), optional :: transA

    !> optional transpose of B matrix (defaults to 'n'), allowed choices are 'n', 'N', 't', 'T', 'c'
    !> and 'C'
    character, intent(in), optional :: transB

    !> specifies the number of columns of the matrix C
    integer, intent(in), optional :: n

    !> specifies the number of rows of the matrix C
    integer, intent(in), optional :: m

    !> specifies the internal number of elements in Op(A)_ik Op(B)_kj
    integer, intent(in), optional :: k

    !> leading dimensions
    integer, intent(in), optional :: lda, ldb, ldc

    integer :: ilda, ildb, ildc
    integer :: in, im, ik
    character :: iTransA, iTransB
    real(r${kind}$p) :: iAlpha, iBeta

    if (present(transA)) then
      iTransA = transA
    else
      iTransA = 'n'
    end if
    if (present(transB)) then
      iTransB = transB
    else
      iTransB = 'n'
    end if

    @:ASSERT(iTransA == 'n' .or. iTransA == 'N' .or. iTransA == 't'&
        & .or. iTransA == 'T' .or. iTransA == 'c' .or. iTransA == 'C')
    @:ASSERT(iTransB == 'n' .or. iTransB == 'N' .or. iTransB == 't'&
        & .or. iTransB == 'T' .or. iTransB == 'c' .or. iTransB == 'C')

    if (present(alpha)) then
      iAlpha = alpha
    else
      iAlpha = 1.0_r${kind}$p
    end if
    if (present(beta)) then
      iBeta = beta
    else
      iBeta = 0.0_r${kind}$p
    end if

  #:for CASES in [('a'), ('b'), ('c')]
    if (present(ld${CASES}$)) then
      ild${CASES}$ = ld${CASES}$
    else
      ild${CASES}$ = size(${CASES}$,dim=1)
    end if
  #:endfor

    if (present(m)) then
      im = m
    else
      if (iTransA == 'n' .or. iTransA == 'N') then
        im = size(A,dim=1)
      else
        im = size(A,dim=2)
      end if
    end if
    if (present(n)) then
      in = n
    else
      in = size(c,dim=2)
    end if
    if (present(k)) then
      ik = k
    else
      if (iTransA == 'n' .or. iTransA == 'N') then
        ik = size(A,dim=2)
      else
        ik = size(A,dim=1)
      end if
    end if

    @:ASSERT(im>0)
    @:ASSERT(in>0)
    @:ASSERT(ik>0)
    @:ASSERT(((ilda>=im).and.(iTransA == 'n' .or. iTransA == 'N'))&
        & .or. (size(a,dim=2)>=im))
    @:ASSERT(ildc>=im)
    @:ASSERT(((size(b,dim=2)>=in).and.(iTransB == 'n' .or. iTransB == 'N'))&
        & .or. (ildb>=in))
    @:ASSERT(size(c,dim=2)>=in)
    @:ASSERT(((size(a,dim=2)>=ik).and.(iTransA == 'n' .or. iTransA == 'N'))&
        & .or. (ilda>=ik))
    @:ASSERT(((ildb>=ik).and.(iTransB == 'n' .or. iTransB == 'N'))&
        & .or. (size(b,dim=2)>=ik))

    call ${kind}$gemm(iTransA,iTransB,im,in,ik,iAlpha,A,ilda,B,ildb,iBeta,C,ildc)

  end subroutine gemm_${suffix}$

#:endfor


#:for suffix, kind in REAL_KIND_PARAMS

  !> Generalized real matrix matrix contraction (Cijl = Aijk * Bkl)
  subroutine gemm332_${suffix}$(C, A, B, alpha, beta, transA, transB)

    !> general matrix output
    real(r${kind}$p), intent(inout), target, contiguous :: C(:,:,:)

    !> general matrix
    real(r${kind}$p), intent(in), target, contiguous :: A(:,:,:)

    !> general matrix
    real(r${kind}$p), intent(in) :: B(:,:)

    !> defaults to 1 if not set
    real(r${kind}$p), intent(in), optional :: alpha

    !> defaults to 0 if not set
    real(r${kind}$p), intent(in), optional :: beta

    !> optional transpose of A matrix (defaults to 'n'), allowed choices are 'n', 'N', 't', 'T',
    !> 'c' and 'C'. Note this acts on the compound index ij
    character, intent(in), optional :: transA

    !> optional transpose of B matrix (defaults to 'n'), allowed choices are 'n', 'N', 't', 'T',
    !> 'c' and 'C'
    character, intent(in), optional :: transB

    real(r${kind}$p), pointer :: pA(:,:), pC(:,:)

    pA(1 : size(A, dim=1) * size(A, dim=2), 1 : size(A, dim=3)) => A
    pC(1 : size(C, dim=1) * size(C, dim=2), 1 : size(C, dim=3)) => C
    call gemm(pC, pA, B, alpha, beta, transA, transB)

  end subroutine gemm332_${suffix}$

#:endfor

end module dftbp_math_blasroutines
