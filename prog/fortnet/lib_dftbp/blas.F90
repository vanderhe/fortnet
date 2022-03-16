!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2021  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#! Note: This module contains preprocessor variable substitutions in subroutine names (${NAME}$)
#! which may break the documentation system. Make sure you preprocess this file before passing it
#! to a source code documentation tool.

#:include 'common.fypp'

!> Interface wrapper for the blas routines.
!>
!> ALL BLAS routines which are called from the main code must be included here.
module dftbp_extlibs_blas

  use dftbp_accuracy, only : rsp, rdp

  public


  interface

    !> performs one of the matrix-vector operations
    !> y := alpha*A*x + beta*y, or y := alpha*A**T*x + beta*y,
    subroutine sgemv(trans, mm, nn, alpha, aa, lda, xx, incx, beta, yy, incy)
      import rsp

      !> should transposition be used
      character, intent(in) :: trans

      !> matrix sizing
      integer, intent(in) :: mm

      !> matrix  size
      integer, intent(in) :: nn

      !> scaling factor
      real(rsp), intent(in) :: alpha

      !> leading matrix dimension
      integer, intent(in) :: lda

      !> matrix A
      real(rsp), intent(in) :: aa(lda, *)

      !> vector
      real(rsp), intent(in) :: xx(*)

      !> stride
      integer, intent(in) :: incx

      !> scale factor
      real(rsp), intent(in) :: beta

      !> vector
      real(rsp), intent(inout) :: yy(*)

      !> stride
      integer, intent(in) :: incy
    end subroutine sgemv


    !> performs one of the matrix-vector operations
    !> y := alpha*A*x + beta*y, or y := alpha*A**T*x + beta*y,
    subroutine dgemv(trans, mm, nn, alpha, aa, lda, xx, incx, beta, yy, incy)
      import rdp

      !> should transposition be used
      character, intent(in) :: trans

      !> matrix sizing
      integer, intent(in) :: mm

      !> matrix  size
      integer, intent(in) :: nn

      !> scale factor
      real(rdp), intent(in) :: alpha

      !> leading matrix dimension
      integer, intent(in) :: lda

      !> matrix A
      real(rdp), intent(in) :: aa(lda, *)

      !> vector
      real(rdp), intent(in) :: xx(*)

      !> stride
      integer, intent(in) :: incx

      !> scaling factor
      real(rdp), intent(in) :: beta

      !> vector
      real(rdp), intent(inout) :: yy(*)

      !> stride
      integer, intent(in) :: incy
    end subroutine dgemv


    !> performs one of the matrix-matrix operations
    !> C := alpha*op( A )*op( B ) + beta*C, where  op( X ) is one of
    !> op( X ) = X, or op( X ) = X**T
    subroutine sgemm(transa, transb, mm, nn, kk, alpha, aa, lda, bb, ldb, beta,&
        & cc, ldc)
      import rsp

      !> On entry, TRANSA specifies the form of op( A ) to be used
      character, intent(in) :: transa

      !> on entry specifies op(B)
      !> should transposition be used
      character, intent(in) :: transb

      !> matrix sizing
      integer, intent(in) :: mm

      !> matrix  size
      integer, intent(in) :: nn

      !> shared index size
      integer, intent(in) :: kk

      !> scaling factor
      real(rsp), intent(in) :: alpha

      !> leading matrix dimension
      integer, intent(in) :: lda

      !> matrix A
      real(rsp), intent(in) :: aa(lda, *)

      !> leading matrix dimension
      integer, intent(in) :: ldb

      !> matrix B
      real(rsp), intent(in) :: bb(ldb, *)

      !> scale factor
      real(rsp), intent(in) :: beta

      !> leading matrix dimension
      integer, intent(in) :: ldc

      !> matrix C
      real(rsp), intent(inout) :: cc(ldc, *)
    end subroutine sgemm


    !> performs one of the matrix-matrix operations
    !> C := alpha*op( A )*op( B ) + beta*C, where  op( X ) is one of
    !> op( X ) = X, or op( X ) = X**T
    subroutine dgemm(transa, transb, mm, nn, kk, alpha, aa, lda, bb, ldb, beta,&
        & cc, ldc)
      import rdp

      !> On entry, TRANSA specifies the form of op( A ) to be used
      character, intent(in) :: transa

      !> on entry specifies op(B)
      !> should transposition be used
      character, intent(in) :: transb

      !> matrix sizing
      integer, intent(in) :: mm

      !> matrix  size
      integer, intent(in) :: nn

      !> shared index size
      integer, intent(in) :: kk

      !> scale factor
      real(rdp), intent(in) :: alpha

      !> leading matrix dimension
      integer, intent(in) :: lda

      !> matrix A
      real(rdp), intent(in) :: aa(lda, *)

      !> leading matrix dimension
      integer, intent(in) :: ldb

      !> matrix B
      real(rdp), intent(in) :: bb(ldb, *)

      !> scaling factor
      real(rdp), intent(in) :: beta

      !> leading matrix dimension
      integer, intent(in) :: ldc

      !> matrix C
      real(rdp), intent(inout) :: cc(ldc, *)
    end subroutine dgemm

  end interface

end module dftbp_extlibs_blas
