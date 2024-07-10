module matrix_inversion
  use, intrinsic :: IEEE_ARITHMETIC
  use kinds
  ! Returns the inverse of a matrix calculated by finding the LU
  ! decomposition.  Depends on LAPACK.
  implicit none
  private :: matrixInverseInPlace, matrixInverseOutOfPlace

  interface matrixInverseSymmetric
    procedure matrixInverseInPlace
    procedure matrixInverseOutOfPlace
  end interface

  interface matrixInverseNonSymmetric
    procedure lapackInverseNonSymmetric
  end interface
contains

  subroutine matrixInverseInPlace(A)
    real(k_real), dimension(:,:), intent(inout) :: A
    real(k_real), dimension(size(A,1),size(A,2)) :: A_copy

    A_copy = A
    call matrixInverseOutOfPlace(A_copy, A)
  end subroutine

  subroutine matrixInverseOutOfPlace(A,Ainv)
    real(k_real), dimension(:,:), intent(in) :: A
    real(k_real), dimension(size(A,1),size(A,2)), intent(out) :: Ainv
    integer :: n, m
    n = size(A,1)
    m = size(A,2)
    if (n/=m) then
      write(*,*) "matrix inverse only works with square matrices !!"
      stop
    end if

    select case (n)
    case (2)
      call matinv2(A,Ainv)
    case (3)
      call matinv3(A,Ainv)
    case (4)
      call matinv4(A,Ainv)
    case default
      call lapackInverseSymmetric(A, Ainv)
    end select

  end subroutine

  subroutine lapackInverseSymmetric(A, Ainv, info_out)
    use tensor_math_mod, only : rotateTensor4, getSymmetricPart
    use print_utils_mod, only : printToScreen
    real(k_real), dimension(:,:), intent(in) :: A
    real(k_real), dimension(size(A,1),size(A,2)), intent(inout) :: Ainv
    integer, optional :: info_out
    integer :: n, info, i, j

    ! External procedures defined in LAPACK
    external DPOTRI
    external DPOTRF

#ifdef __DEBUG__
      ! check if matrix is symmetric
      real(k_real), dimension(size(A,1),size(A,2)) :: A_skew
      A_skew = (A-transpose(A))*0.5_k_real

      if (.not.all(abs(A_skew).lt.1e-6_k_real)) then
        call printToScreen(A, "Matrix to invert")
        call printToScreen(A_skew, "skew part")
        error stop "you can not use lapackInverseSymmetric for a non symmetric matrix"
      endif
#endif

    ! Store A in Ainv to prevent it from being overwritten by LAPACK
    Ainv = A
    n = size(A,1)

    ! DGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
    ! call DGETRF(n, n, Ainv, n, ipiv, info)
    call DPOTRF("U", n, Ainv, n, info)

    if (info /= 0) then
      call printToScreen(A, "lapackInverseSymmetric: Singualr Matrix")
      write(*,*) "lapackInverseSymmetric: lapack info code: ", info
    end if

    ! DGETRI computes the inverse of a matrix using the LU factorization
    ! computed by DGETRF.
    ! call DGETRI(n, Ainv, n, ipiv, work, n, info)
    call DPOTRI("U", n, Ainv, n, info)

    do j=1,n-1
      do i=j+1,n
        Ainv(i,j) = Ainv(j,i)
    enddo;enddo

    if (info /= 0) then
      call printToScreen(A, "LU factorization failed")
    end if
    if (present(info_out)) info_out = info
  end subroutine

  subroutine lapackInverseNonSymmetric(A, Ainv, info_out)
    use tensor_math_mod, only : rotateTensor4, getSymmetricPart
    use print_utils_mod, only : printToScreen
    real(k_real), dimension(:,:), intent(in) :: A
    real(k_real), dimension(size(A,1),size(A,2)), intent(inout) :: Ainv
    integer, optional :: info_out

    real(k_real), dimension(size(A,1)) :: work  ! work array for LAPACK
    integer, dimension(size(A,1)) :: ipiv   ! pivot indices
    integer :: n, info

    ! External procedures defined in LAPACK
    external DGETRF
    external DGETRI

#ifdef __DEBUG__
  ! check if matrix is symmetric
  real(k_real), dimension(size(A,1),size(A,2)) :: A_skew
  integer :: i,j
  A_skew = A-transpose(A)
  do j=1,size(A,2); do i=j,size(A,1)
    if (A(i,j).ne.0._k_real) A_skew(i,j) = A_skew(i,j)/A(i,j)
  enddo; enddo
  if (all(abs(A_skew).lt.1e-6)) then
    call printToScreen(A, "Matrix to invert is symmetric. You might want to consider using lapackInverseSymmetric")
    call printToScreen(A-transpose(A), "skew part")
  endif
#endif

    ! Store A in Ainv to prevent it from being overwritten by LAPACK
    Ainv = A
    n = size(A,1)

    ! DGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
    call DGETRF(n, n, Ainv, n, ipiv, info)

    if (info /= 0) then
      call printToScreen(A, "lapackInverseNonSymmetric Singualr Matrix")
    end if

    ! DGETRI computes the inverse of a matrix using the LU factorization
    ! computed by DGETRF.
    call DGETRI(n, Ainv, n, ipiv, work, n, info)

    if (info /= 0) then
      call printToScreen(A, "lapackInverseNonSymmetric LU factorization failed")
    end if
    if (present(info_out)) info_out = info
  end subroutine

  subroutine matinv2(A,Ainv)
    !! Performs a direct calculation of the inverse of a 2×2 matrix.
    real(k_real), intent(in) :: A(2,2)      !! Matrix
    real(k_real)             :: Ainv(2,2)   !! Inverse matrix
    real(k_real)             :: detinv

    ! Calculate the 1 / determinant of the matrix
    detinv = 1./(A(1,1)*A(2,2) - A(1,2)*A(2,1))
    call errorSingularMatrix(detinv, A)
    ! Calculate the inverse of the matrix
    Ainv(1,1) = +detinv * A(2,2)
    Ainv(2,1) = -detinv * A(2,1)
    Ainv(1,2) = -detinv * A(1,2)
    Ainv(2,2) = +detinv * A(1,1)
  end subroutine

  subroutine matinv3(A,Ainv)
  !! Performs a direct calculation of the inverse of a 3×3 matrix.
  real(k_real), intent(in) :: A(3,3)   !! Matrix
  real(k_real), intent(inout):: Ainv(3,3)   !! Inverse matrix
  real(k_real) :: detinv

  ! Calculate the inverse determinant of the matrix
  detinv = 1./(A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)&
            - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)&
            + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))
  call errorSingularMatrix(detinv, A)

  ! Calculate the inverse of the matrix
  Ainv(1,1) = +detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
  Ainv(2,1) = -detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
  Ainv(3,1) = +detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
  Ainv(1,2) = -detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
  Ainv(2,2) = +detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
  Ainv(3,2) = -detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
  Ainv(1,3) = +detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
  Ainv(2,3) = -detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
  Ainv(3,3) = +detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))
end subroutine

subroutine matinv4(A, Ainv)
  !! Performs a direct calculation of the inverse of a 4×4 matrix.
  real(k_real), intent(in) :: A(4,4)   !! Matrix
  real(k_real), intent(inout) :: Ainv(4,4)   !! Inverse matrix
  real(k_real)             :: detinv

  ! Calculate the inverse determinant of the matrix
  detinv = &
    1./(A(1,1)*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-&
                      A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))&
     - A(1,2)*(A(2,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,1)-&
     A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))&
     + A(1,3)*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-&
     A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))&
     - A(1,4)*(A(2,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(2,2)*(A(3,3)*A(4,1)-&
     A(3,1)*A(4,3))+A(2,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1))))

  call errorSingularMatrix(detinv, A)

  ! Calculate the inverse of the matrix
  Ainv(1,1) = detinv*(A(2,2)*(A(3,3)*A(4,4)-&
  A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+&
  A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))
  Ainv(2,1) = detinv*(A(2,1)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+&
  A(2,3)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(2,4)*(A(3,3)*A(4,1)-A(3,1)*A(4,3)))
  Ainv(3,1) = detinv*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+&
  A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
  Ainv(4,1) = detinv*(A(2,1)*(A(3,3)*A(4,2)-A(3,2)*A(4,3))+&
  A(2,2)*(A(3,1)*A(4,3)-A(3,3)*A(4,1))+A(2,3)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
  Ainv(1,2) = detinv*(A(1,2)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+&
  A(1,3)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(1,4)*(A(3,3)*A(4,2)-A(3,2)*A(4,3)))
  Ainv(2,2) = detinv*(A(1,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+&
  A(1,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(1,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))
  Ainv(3,2) = detinv*(A(1,1)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+&
  A(1,2)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(1,4)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
  Ainv(4,2) = detinv*(A(1,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+&
  A(1,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(1,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
  Ainv(1,3) = detinv*(A(1,2)*(A(2,3)*A(4,4)-A(2,4)*A(4,3))+&
  A(1,3)*(A(2,4)*A(4,2)-A(2,2)*A(4,4))+A(1,4)*(A(2,2)*A(4,3)-A(2,3)*A(4,2)))
  Ainv(2,3) = detinv*(A(1,1)*(A(2,4)*A(4,3)-A(2,3)*A(4,4))+&
  A(1,3)*(A(2,1)*A(4,4)-A(2,4)*A(4,1))+A(1,4)*(A(2,3)*A(4,1)-A(2,1)*A(4,3)))
  Ainv(3,3) = detinv*(A(1,1)*(A(2,2)*A(4,4)-A(2,4)*A(4,2))+&
  A(1,2)*(A(2,4)*A(4,1)-A(2,1)*A(4,4))+A(1,4)*(A(2,1)*A(4,2)-A(2,2)*A(4,1)))
  Ainv(4,3) = detinv*(A(1,1)*(A(2,3)*A(4,2)-A(2,2)*A(4,3))+&
  A(1,2)*(A(2,1)*A(4,3)-A(2,3)*A(4,1))+A(1,3)*(A(2,2)*A(4,1)-A(2,1)*A(4,2)))
  Ainv(1,4) = detinv*(A(1,2)*(A(2,4)*A(3,3)-A(2,3)*A(3,4))+&
  A(1,3)*(A(2,2)*A(3,4)-A(2,4)*A(3,2))+A(1,4)*(A(2,3)*A(3,2)-A(2,2)*A(3,3)))
  Ainv(2,4) = detinv*(A(1,1)*(A(2,3)*A(3,4)-A(2,4)*A(3,3))+&
  A(1,3)*(A(2,4)*A(3,1)-A(2,1)*A(3,4))+A(1,4)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)))
  Ainv(3,4) = detinv*(A(1,1)*(A(2,4)*A(3,2)-A(2,2)*A(3,4))+&
  A(1,2)*(A(2,1)*A(3,4)-A(2,4)*A(3,1))+A(1,4)*(A(2,2)*A(3,1)-A(2,1)*A(3,2)))
  Ainv(4,4) = detinv*(A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))+&
  A(1,2)*(A(2,3)*A(3,1)-A(2,1)*A(3,3))+A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1)))
end subroutine

subroutine errorSingularMatrix(det_inv, a)
  use print_utils_mod, only:  printToScreen
  real(k_real), dimension(:,:), intent(in) :: A
  real(k_real), intent(in) :: det_inv
  if (.not.ieee_is_finite(det_inv)) then
    write(*,*) "matrix is singular ", det_inv
    call printToScreen(A, "singular matrix")
    stop
  end if
end subroutine

end module matrix_inversion
