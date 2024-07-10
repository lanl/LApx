module linear_interpolation_mod
  use kinds
  use polymorphic_dtype_array_mod, only : dtype_array_ptr
  implicit none

  type :: piecewise_linear_interpolation
    integer :: n_vals = 0
    integer :: n_bins = 0
    real(k_real), pointer, dimension(:) :: x_data => null()
    real(k_real), pointer, dimension(:) :: y_data => null()
    real(k_real), pointer, dimension(:,:) :: interp_coeff => null()
    real(k_real) :: old_x_value
    real(k_real) :: old_y_value
    logical :: initialized = .FALSE.
  contains

    procedure :: init => initLinearInterpolator
    procedure :: checkInRange => LinearInterpolatorCheckInRange
    procedure :: interpolate => LinearInterpolatorInterpolate
    procedure :: initCoefficients => LinearInterpolatorInitCoefficients
    procedure :: computeY => LinearInterpolatorComputeY
    procedure :: FindXBin => LinearInterpolatorFindXBin
  end type piecewise_linear_interpolation

  type, extends(dtype_array_ptr) :: piecewise_linear_interpolation_array_type
  contains
    procedure :: addElement => addLinearInterpolator
    procedure :: getElementPtr => getLinearInterpolatorPtr
  end type
contains


  subroutine initLinearInterpolator(this, n_vals, x_vals, y_vals)
    implicit none
    class(piecewise_linear_interpolation), intent(inout) :: this
    integer, intent(in) ::  n_vals
    real(k_real), dimension(:), intent(in) :: x_vals, y_vals
    integer :: i

    if (n_vals< 2) error stop "initLinearInterpolator: n_vals < 2 . Abort!"
    if (size(x_vals).ne.n_vals)  error stop "initLinearInterpolator: size(x_vals).ne.n_vals . Abort!"
    if (size(y_vals).ne.n_vals)  error stop "initLinearInterpolator: size(y_vals).ne.n_vals . Abort!"

    do i=2, n_vals
      if (x_vals(i).lt.x_vals(i-1)) error stop "initLinearInterpolator: x_vals are not sorted . Abort!"
      if (x_vals(i).eq.x_vals(i-1)) error stop "initLinearInterpolator: there are two or more x_vals with the same value. Abort!"
    enddo

    allocate(this%x_data(n_vals), this%y_data(n_vals)  )

    this%n_vals = n_vals
    this%x_data = x_vals
    this%y_data = y_vals
    this%old_x_value = x_vals(1)
    this%old_y_value = y_vals(1)

    this%n_bins = n_vals-1
    call this%initCoefficients()
    this%initialized = .TRUE.
  end subroutine

  subroutine LinearInterpolatorCheckInRange(this, x)
    class(piecewise_linear_interpolation), intent(in) :: this
    real(k_real), intent(in) :: x
    real(k_Real), parameter :: xtol = 1e-6
    if (x.lt.(this%x_data(1)-xtol)) error stop "LinearInterpolatorCheckInRange the provided x is below the interpolation limit "
    if (x.gt.(this%x_data(this%n_vals)+xtol)) error stop "LinearInterpolatorCheckInRange the provided x is greater than interpolation limit "
  end subroutine

  subroutine LinearInterpolatorInterpolate(this, x, y)
    implicit none
    class(piecewise_linear_interpolation), intent(inout) :: this
    real(k_real), intent(in) :: x
    real(k_real), intent(out) :: y
    integer :: bin

    if (x.eq.this%old_x_value ) then
      y = this%old_y_value
    else
      call this%checkInRange(x)
      call this%FindXBin(x, bin)
      call this%computeY(x, bin, y)
      this%old_y_value = y
      this%old_x_value = x
    endif

  end subroutine

  subroutine LinearInterpolatorComputeY(this, x, bin, y)
    implicit none
    class(piecewise_linear_interpolation), intent(inout) :: this
    real(k_real), intent(in) :: x
    integer, intent(in) :: bin
    real(k_real), intent(out) :: y

    y = this%interp_coeff(bin,1) * x + this%interp_coeff(bin,2)

  end subroutine


  subroutine LinearInterpolatorInitCoefficients(this)
    implicit none
    class(piecewise_linear_interpolation), intent(inout) :: this
    real(k_real) :: dx, dy
    integer :: i

    allocate(this%interp_coeff(this%n_bins, 2))

    do i=1, this%n_bins
      dx = this%x_data(i+1)-this%x_data(i)
      dy = this%y_data(i+1)-this%y_data(i)
      this%interp_coeff(i,1) = dy/dx
      this%interp_coeff(i,2) = -this%x_data(i)*this%interp_coeff(i,1)+this%y_data(i)
    enddo

  end subroutine


  subroutine LinearInterpolatorFindXBin(this, x, bin)
    implicit none
    class(piecewise_linear_interpolation), intent(in) :: this
    real(k_real), intent(in) :: x
    integer, intent(out) :: bin
    integer :: i

    bin = -1

    do i=1,this%n_bins
      if (x.lt.this%x_data(i+1)) then
        bin = i
        exit
      endif
    enddo

    if (bin<1) error stop "LinearInterpolatorFindXBin cannot find bin"

  end subroutine

  !------------------------------------------------------------------------------!
  !                  LINEAR INTERPOLATOR ARRAY SUBROUTINE                          !
  !------------------------------------------------------------------------------!
  subroutine addLinearInterpolator(this, new_element)
    class(piecewise_linear_interpolation_array_type), intent(inout) :: this
    class(piecewise_linear_interpolation), pointer, intent(inout) :: new_element

    call this%extend()
    this%all_pt(this%n_elements)%pt => new_element
    nullify(new_element)
  end subroutine

  subroutine getLinearInterpolatorPtr(this, idx, ptr)
    class(piecewise_linear_interpolation_array_type), intent(in) :: this
    integer, intent(in) :: idx
    class(piecewise_linear_interpolation), pointer, intent(out) :: ptr

    call this%checkElementExist(idx)

    associate(elem => this%all_pt(idx)%pt)
    select type(elem)
      class is (piecewise_linear_interpolation)
        ptr => elem
      class default
        error stop "getLinearInterpolatorPtr wrong class"
    end select
    end associate
  end subroutine


end module
