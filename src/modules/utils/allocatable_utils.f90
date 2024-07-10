module allocatable_utils_module
contains
  subroutine concatenateIntegerAllocatableArrays(a, b, concatenated)
    integer, dimension(:), allocatable, intent(in) :: a, b
    integer, dimension(:), allocatable, intent(inout) :: concatenated
    integer :: size_a, size_b, size_a_plus_b
    integer :: i
    size_a=0
    size_b=0
    if (allocated(a)) size_a = size(a)
    if (allocated(b)) size_b = size(b)
    size_a_plus_b = size_a + size_b
    if (allocated(concatenated)) deallocate(concatenated)
    allocate(concatenated(size_a_plus_b))

    if (size_a_plus_b > 0) then
      if (size_a > 0) then
        do i=1,size_a
          concatenated(i) = a(i)
        end do
      end if

      if (size_b > 0) then
        do i=1,size_b
          concatenated(i+size_a) = b(i)
        end do
      end if
    end if
  end subroutine

  subroutine copyIntegerAllocatableArrays(a, a_copy)
    integer, dimension(:), allocatable, intent(in) :: a
    integer, dimension(:), allocatable, intent(inout) :: a_copy
    integer :: i, size_a
    if (allocated(a)) then
      size_a = size(a)
      allocate(a_copy(size_a))
      if (size_a > 0) then
        do i=1, size_a
          a_copy(i) = a(i)
        end do
      end if
    else
      allocate(a_copy(0))
    end if
  end subroutine
end module
