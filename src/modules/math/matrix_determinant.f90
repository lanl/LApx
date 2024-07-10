module matrix_determinant_mod
  use, intrinsic :: IEEE_ARITHMETIC
  use kinds
  ! Returns the inverse of a matrix calculated by finding the LU
  ! decomposition.  Depends on LAPACK.
  implicit none
  private :: matdet2, matdet3

contains

  subroutine matrixDeterminat(A,det)
    real(k_real), dimension(:,:), intent(in) :: A
    real(k_real), intent(out) :: det
    integer :: n, m
    n = size(A,1)
    m = size(A,2)
    if (n/=m) then
      write(*,*) "matrix inverse only works with square matrices !!"
      stop
    end if

    select case (n)
    case (2)
      call matdet2(A,det)
    case (3)
      call matdet3(A,det)
    case default
      error stop "matric determinant implemented only for rank 2 and 3 matrices"
      ! call lapackDeterminant(A, Ainv)
    end select

  end subroutine

  subroutine matdet2(A, det)
    use tensor_math_mod, only : rotateTensor4, getSymmetricPart
    use print_utils_mod, only : printToScreen
    real(k_real), dimension(2,2), intent(in) :: A
    real(k_real), intent(inout) :: det
      det = A(1,1)*A(2,2) - A(1,2)*A(2,1)
  end subroutine

  subroutine matdet3(A, det)
    use tensor_math_mod, only : rotateTensor4, getSymmetricPart
    use print_utils_mod, only : printToScreen
    real(k_real), dimension(3,3), intent(in) :: A
    real(k_real), intent(inout) :: det
      det = A(1,1)*A(2,2)*A(3,3) + A(1,2)*A(2,3)*A(3,1) + A(1,3)*A(2,1)*A(3,2) &
            - A(1,3)*A(2,2)*A(3,1) - A(1,2)*A(2,1)*A(3,3) - A(1,1)*A(2,3)*A(3,2)
  end subroutine

end module matrix_determinant_mod
