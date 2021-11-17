module lapack_interface
  use, intrinsic :: iso_fortran_env, only : real64
  implicit none

  interface
    subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
      import
      integer, intent(in) :: info, lda, ldb, n, nrhs
      integer, intent(out):: ipiv(*)
      real(real64), intent(inout) :: A(LDA, *), B(LDB, *)
    end subroutine dgesv
  end interface

end module lapack_interface
