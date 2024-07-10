module polymorphic_dtype_array_mod

! This dtype_array_ptr type allows to create and manage arrays containing ANY
! derived type pointer, and uses unlimited polymorphism (i.e. the they array can
! be composed by any derived type)
! we use it because we want to reuse as much code as possible at least until
! we have a reasonable framework. In the future we might want to write specific
! types for each required type.

public dtype_array_ptr
private dtype_ptr

! an unlimited polymorphic data type
type :: dtype_ptr
  class(*), pointer :: pt => null()
end type

! an array of unlimited polymorphic data
type :: dtype_array_ptr
  type(dtype_ptr), dimension(:), pointer :: all_pt => null()
  integer :: n_elements = 0
  logical :: initialized=.false.
contains
  ! initializes the Dtype array
  procedure :: init => InitDTypeArray
  ! extend (i.e. resize teh DType array)
  procedure :: extend => ExtendDTypeArray
  ! same as above but leaves the free elements are at the beginning
  procedure :: extendFirstFree => ExtendDTypeArrayFirstFree
  ! nullify and deallocate a Dtype array
  procedure :: reset => resetDTypeArray
  ! check taht an element with the given index exists (plus other checks)
  procedure :: checkElementExist
end type

contains
  subroutine checkElementExist(this, idx)
    class(dtype_array_ptr), intent(in) :: this
    integer :: idx

    if (.not.(this%initialized)) error stop "array not initialized"
    if (idx<1) error stop "the provided index is < 1. Abort!"
    if (idx>this%n_elements) error stop "the provided index is > than the number of elements. Abort!"
    if (.not.(associated(this%all_pt(idx)%pt))) error stop "pointer not associated. Abort"

  end subroutine

  subroutine InitDTypeArray(this, n_add_in)
    class(dtype_array_ptr), intent(inout) :: this
    integer, optional :: n_add_in
    integer :: n_add

    if (this%initialized) error stop "DTypeArray already initialized. Abort!"

    if (present(n_add_in)) then
      n_add = n_add_in
    else
      n_add = 1
    end if

    allocate(this%all_pt(n_add))
    this%n_elements = n_add
    this%initialized =.true.
  end subroutine

  subroutine resetDTypeArray(this)
    class(dtype_array_ptr), intent(inout) :: this
    integer :: i
    if (.not.(this%initialized)) error stop "DTypeArray cannot reset a non-initialized DTypeArray"
    do i=1,this%n_elements
      nullify(this%all_pt(i)%pt)
    enddo
    deallocate(this%all_pt)
    this%n_elements = 0
    this%initialized = .false.
  end subroutine

  subroutine ExtendDTypeArray(this, n_add_in)
    class(dtype_array_ptr), intent(inout) :: this
    integer, intent(in), optional :: n_add_in
    integer :: i, n_add, n_elem_old
    type(dtype_array_ptr) :: dtype_array_temp

    if (present(n_add_in)) then
      n_add = n_add_in
    else
      n_add = 1
    end if

    if (.not.(this%initialized)) then
      call this%init(n_add)
    else
      call dtype_array_temp%init(this%n_elements)
      do i=1,this%n_elements
        dtype_array_temp%all_pt(i)%pt => this%all_pt(i)%pt
      enddo
      n_elem_old = this%n_elements

      call this%reset()
      call this%init(n_add+n_elem_old)

      do i=1, n_elem_old
        this%all_pt(i)%pt => dtype_array_temp%all_pt(i)%pt
      enddo
      call dtype_array_temp%reset()
    endif

  end subroutine

  subroutine ExtendDTypeArrayFirstFree(this, n_add_in)
    class(dtype_array_ptr), intent(inout) :: this
    integer, intent(in), optional :: n_add_in
    integer :: i, n_add, n_elem_old
    type(dtype_array_ptr) :: dtype_array_temp

    if (present(n_add_in)) then
      n_add = n_add_in
    else
      n_add = 1
    end if

    if (.not.(this%initialized)) then
      call this%init(n_add)
    else
      call dtype_array_temp%init(this%n_elements)
      do i=1,this%n_elements
        dtype_array_temp%all_pt(i)%pt => this%all_pt(i)%pt
      enddo
      n_elem_old = this%n_elements

      call this%reset()
      call this%init(n_add+n_elem_old)

      do i=1, n_elem_old
        this%all_pt(n_add+i)%pt => dtype_array_temp%all_pt(i)%pt
      enddo
      call dtype_array_temp%reset()
    endif

  end subroutine

end module
