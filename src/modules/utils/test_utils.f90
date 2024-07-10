module test_utils_mod
  use, intrinsic :: IEEE_ARITHMETIC
  use kinds
  implicit none

  private :: test_is_finite_scalar, test_is_finite_array, test_is_finite_matrix
  private :: fuzzy_equal_scalar, fuzzy_equal_vectors, fuzzy_equal_matrix
  private :: equal_scalar, equal_vectors, equal_matrix

  public :: test_is_finite, fuzzy_equal, equal

  interface test_is_finite
    procedure test_is_finite_scalar
    procedure test_is_finite_array
    procedure test_is_finite_matrix
  end interface

  interface fuzzy_equal
    procedure fuzzy_equal_scalar
    procedure fuzzy_equal_vectors
    procedure fuzzy_equal_matrix
  end interface

  interface equal
    procedure equal_scalar
    procedure equal_vectors
    procedure equal_matrix
  end interface

  real(k_real), parameter :: fuzzy_abs_tol=1e-12, fuzzy_rel_tol=1e-6
contains

  subroutine test_is_finite_scalar(x, name)
    implicit none
    real(k_real), intent(in) :: x
    character(len=*), intent(in) :: name
#ifdef __DEBUG__
    if (IEEE_IS_NAN(x)) then
    write(*,*) "variable ", name, " contains nans"
    write(*,*) x
    stop
    endif
    if (.not.IEEE_IS_FINITE(x)) then
    write(*,*) "variable ", name, " is not finite"
    write(*,*) x
    stop
    endif
#endif
  end subroutine

  subroutine test_is_finite_array(x, name)
    implicit none
    real(k_real), intent(in) :: x(:)
    character(len=*), intent(in) :: name
#ifdef __DEBUG__
    if (any(IEEE_IS_NAN(x))) then
    write(*,*) "variable ", name, " contains nans"
    write(*,*) x
    stop
    endif
    if (.not.any(IEEE_IS_FINITE(x))) then
    write(*,*) "variable ", name, " is not finite"
    write(*,*) x
    stop
    endif
#endif
  end subroutine

  subroutine test_is_finite_matrix(x, name)
    implicit none
    real(k_real), intent(in) :: x(:,:)
    character(len=*), intent(in) :: name
#ifdef __DEBUG__
    if (any(IEEE_IS_NAN(x))) then
    write(*,*) "variable ", name, " contains nans"
    write(*,*) x
    stop
    endif
    if (.not.any(IEEE_IS_FINITE(x))) then
    write(*,*) "variable ", name, " is not finite"
    write(*,*) x
    stop
    endif
#endif
  end subroutine

  subroutine fuzzy_equal_scalar(x1, x2, name_x1, name_x2)
    implicit none
    real(k_real), intent(in) :: x1, x2
    character(len=*), intent(in) :: name_x1, name_x2
#ifdef __DEBUG__
    real(k_real) ::  abs_err, rel_err, min_err

    call test_is_finite(x1, name_x1)
    call test_is_finite(x2, name_x2)

    abs_err = abs(x1-x2)

    if (x1.ne.0._k_real) then
      rel_err = abs_err/abs(x1)
    else
      rel_err = 1e6
    endif

    min_err = min(abs_err, rel_err)

    if (min_err>fuzzy_abs_tol) then
      write(*,*) "var ", name_x1, " != ", "var ", name_x2
      write(*,*) name_x1, ": ", x1
      write(*,*) name_x2, ": ", x2
      write(*,*) "absolute error: ", abs_err
      write(*,*) "relative error: ", rel_err
      write(*,*) "min error: ", min_err
      stop
    endif
#endif
  end subroutine

  subroutine fuzzy_equal_vectors(x1, x2, name_x1, name_x2)
    implicit none
    real(k_real), intent(in) :: x1(:), x2(:)
    character(len=*), intent(in) :: name_x1, name_x2
#ifdef __DEBUG__
    real(k_real), dimension(:), ALLOCATABLE :: abs_err, rel_err, min_err
    integer :: i

    call test_is_finite(x1, name_x1)
    call test_is_finite(x2, name_x2)

    allocate(abs_err(size(x1)), rel_err(size(x1)), min_err(size(x1)))

    abs_err = abs(x1-x2)

      do i=1,size(x1)
        if (x1(i).ne.0._k_real) then
          rel_err(i) = abs_err(i)/abs(x1(i))
        else
          rel_err(i) = 1e6
        endif
      enddo

    min_err = min(abs_err, rel_err)

    if (any(min_err>fuzzy_abs_tol)) then
      write(*,*) "var ", name_x1, " != ", "var ", name_x2
      write(*,*) name_x1, ": ", x1
      write(*,*) name_x2, ": ", x2
      write(*,*) "absolute error: ", abs_err
      write(*,*) "relative error: ", rel_err
      write(*,*) "min error: ", min_err
      stop
    endif
#endif
  end subroutine

  subroutine fuzzy_equal_matrix(x1, x2, name_x1, name_x2)
    implicit none
    real(k_real), intent(in) :: x1(:,:), x2(:,:)
    character(len=*), intent(in) :: name_x1, name_x2
#ifdef __DEBUG__
    real(k_real), dimension(:,:), ALLOCATABLE :: abs_err, rel_err, min_err
    integer, dimension(2) :: nx_ny
    integer :: i,j, nx, ny

    call test_is_finite(x1, name_x1)
    call test_is_finite(x2, name_x2)
    nx_ny = shape(x1)
    nx = nx_ny(1)
    ny = nx_ny(2)

    allocate(abs_err(nx,ny), rel_err(nx,ny), min_err(nx,ny))

    abs_err = abs(x1-x2)

    do j=1,ny
      do i=1,nx
        if (x1(i,j).ne.0._k_real) then
          rel_err(i,j) = abs_err(i,j)/abs(x1(i,j))
        else
          rel_err(i,j) = 1e6
        endif
      enddo
    enddo

    min_err = min(abs_err, rel_err)

    if (any(min_err>fuzzy_abs_tol)) then
      write(*,*) "var ", name_x1, " != ", "var ", name_x2
      write(*,*) name_x1, ": ", x1
      write(*,*) name_x2, ": ", x2
      write(*,*) "absolute error: ", abs_err
      write(*,*) "relative error: ", rel_err
      write(*,*) "min error: ", min_err
      stop
    endif
#endif
  end subroutine


  subroutine equal_scalar(x1, x2, name_x1, name_x2)
    implicit none
    real(k_real), intent(in) :: x1, x2
    character(len=*), intent(in) :: name_x1, name_x2
#ifdef __DEBUG__
    real(k_real) ::  abs_err

    call test_is_finite(x1, name_x1)
    call test_is_finite(x2, name_x2)

    abs_err = abs(x1-x2)

    if (abs_err.ne.0._k_real) then
      write(*,*) "var ", name_x1, " != ", "var ", name_x2
      write(*,*) name_x1, ": ", x1
      write(*,*) name_x2, ": ", x2
      write(*,*) "absolute error: ", abs_err
      stop
    endif
#endif
  end subroutine

  subroutine equal_vectors(x1, x2, name_x1, name_x2)
    implicit none
    real(k_real), intent(in) :: x1(:), x2(:)
    character(len=*), intent(in) :: name_x1, name_x2
#ifdef __DEBUG__
    real(k_real), dimension(:), ALLOCATABLE :: abs_err

    call test_is_finite(x1, name_x1)
    call test_is_finite(x2, name_x2)

    allocate(abs_err(size(x1)))

    abs_err = abs(x1-x2)

    if (any(abs_err.ne.0._k_real)) then
      write(*,*) "var ", name_x1, " != ", "var ", name_x2
      write(*,*) name_x1, ": ", x1
      write(*,*) name_x2, ": ", x2
      write(*,*) "absolute error: ", abs_err
      stop
    endif
#endif
  end subroutine

  subroutine equal_matrix(x1, x2, name_x1, name_x2)
    implicit none
    real(k_real), intent(in) :: x1(:,:), x2(:,:)
    character(len=*), intent(in) :: name_x1, name_x2
#ifdef __DEBUG__
    real(k_real), dimension(:,:), ALLOCATABLE :: abs_err
    integer, dimension(2) :: nx_ny
    integer :: nx, ny

    call test_is_finite(x1, name_x1)
    call test_is_finite(x2, name_x2)
    nx_ny = shape(x1)
    nx = nx_ny(1)
    ny = nx_ny(2)

    allocate(abs_err(nx,ny))

    abs_err = abs(x1-x2)

    if (any(abs_err.ne.0._k_real)) then
      write(*,*) "var ", name_x1, " != ", "var ", name_x2
      write(*,*) name_x1, ": ", x1
      write(*,*) name_x2, ": ", x2
      write(*,*) "absolute error: ", abs_err
      stop
    endif
#endif
  end subroutine

end module
