module string_module
  use kinds
  implicit none


  ! A smart string type always allocating a string of the proper length.
  ! It discards trailing and leading spaces.
  type :: string_type
      character(len=:), allocatable, private :: str
      logical, private :: is_empty = .true.
      integer :: string_length = -1
  contains
    ! procedure to set the string
    procedure :: setString => stringTypeSetString
    ! procedure to get the string
    procedure :: getString => stringTypeGetString
    !get the stored string length
    procedure :: getStringLength 
    !set the stored string length
    procedure, private :: setStringLength 
    ! procedure to know if the string is empty
    procedure :: isEmpty => stringTypeIsEmpty
    ! procedure to know if the string is empty
    procedure :: resetString => stringTypeResetString
    ! procedure to know if the string starts width certain characters
    procedure :: startsWith => stringTypeStartsWith
    ! procedure to know if the string ends width certain characters
    procedure :: endsWith => stringTypeEndsWith
    ! procedure to know if the string ends width certain characters
    procedure :: equal => stringTypeIsEqual
    ! procedure splitting a string based on a separetor character
    procedure :: splitString => stringTypeSplitString
    !convert a string to an integer, logical, or real number
    generic, public :: stringToVar => stringToVarLogical, stringToVarReal, stringToVarInteger
    procedure, private :: stringToVarLogical
    procedure, private :: stringToVarReal
    procedure, private :: stringToVarInteger
  end type string_type

  ! A container of smart strings
  type :: string_array
    ! we are using an array, hence we preallocate a given size.
    ! If space is not enough we extend the array. The standard allocate batch size is 10
    type(string_type), dimension(:), allocatable :: strings
    integer, private :: num_string = 0
  contains
    ! Procedure to add a new string to the list
    procedure :: addString => stringArrayAddString
    ! procedure to ask the container how many strings it holds
    procedure :: getNumStrings => stringArrayGetNumString
    ! procedure to retrieve a string given its position in the array
    procedure :: getStringByIndex => stringArrayGetStringByIndex
    ! procedure to retrieve a string given its position in the array
    procedure :: setStringByIndex => stringArraySetStringByIndex
    ! reset
    procedure :: reset => stringArrayReset
  end type

contains
  subroutine stringToVarLogical(this, var)
    class(string_type), intent(in) :: this
    logical, intent(inout) :: var
    read(this%str,*) var
  end subroutine

  subroutine stringToVarReal(this, var)
    class(string_type), intent(in) :: this
    real(k_real), intent(inout) :: var
    integer :: ierr
    read(this%str,*, iostat=ierr) var
    if (ierr.ne.0) then
      write(*,*) "I should find a real number but instead i find ``", this%str, "''"
      stop
    endif

  end subroutine

  subroutine stringToVarInteger(this, var)
    class(string_type), intent(in) :: this
    integer, intent(inout) :: var
    read(this%str,*) var
  end subroutine

  function stringTypeIsEmpty(this) result(is_empty)
    class(string_type), intent(in) :: this
    logical :: is_empty
    is_empty = this%is_empty
  end function

  subroutine stringTypeResetString(this)
    class(string_type), intent(inout) :: this
    this%str=""
    this%is_empty=.true.
    this%string_length = -1
  end subroutine

  function stringTypeStartsWith(this, str) result(starts_with)
    class(string_type), intent(in) :: this
    character(len=*) :: str
    logical :: starts_with
    integer :: n_char_to_compare

    starts_with = .false.
    n_char_to_compare = len(trim(adjustl(str)))
    if (n_char_to_compare<=0) error stop "nothing to compare the string with, check your input!"

    if (.not.(this%is_empty)) then
      if (len(this%str)>= n_char_to_compare) then
        if (this%str(1:n_char_to_compare)==str(1:n_char_to_compare)) starts_with = .true.
      endif
    endif
  end function

  function stringTypeIsEqual(this, str) result(strings_are_equal)
    class(string_type), intent(in) :: this
    character(len=*) :: str
    logical :: strings_are_equal
    integer :: len_input_string

    strings_are_equal = .false.
    len_input_string = len(trim(adjustl(str)))
    if (len_input_string<0) error stop "nothing to compare the string with, check your input!"

    if (len_input_string==0) then
      if (this%is_empty)  strings_are_equal = .true.
      return
    endif

    if (.not.(this%is_empty)) then
      if (len(this%str)== len_input_string) then
        if (this%str==trim(adjustl(str))) strings_are_equal = .true.
      endif
    endif
  end function

  function stringTypeEndsWith(this, str) result(ends_with)
    class(string_type), intent(in) :: this
    character(len=*) :: str
    logical :: ends_with
    integer :: n_char_to_compare

    ends_with = .false.
    n_char_to_compare = len(trim(adjustl(str)))
    if (n_char_to_compare<=0) error stop "nothing to compare the string with, check your input!"

    if (.not.(this%is_empty)) then
      if (len(this%str)>= n_char_to_compare) then
        if (this%str(len(this%str)-n_char_to_compare+1:len(this%str))==str(1:n_char_to_compare)) ends_with = .true.
      endif
    endif
  end function

  function stringTypeSplitString(this, separator) result(split_strings_array)
    class(string_type), intent(in) :: this
    character(len=1), intent(in) :: separator
    type(string_array) :: split_strings_array
    integer :: i, str_start, str_end

    str_start = 1
    str_end = 0

    if (this%isEmpty()) then
      call split_strings_array%addString("")
      return
    endif

    do i=2,len(this%str)
      if ((this%str(i-1:i-1) == separator).and.(this%str(i:i) /= separator)) then
        str_end = i-1
        call split_strings_array%addString(this%str(str_start:str_end))
        str_start = i
      end if
    end do

    if (str_end<len(this%str)) then
      str_end = len(this%str)
      call split_strings_array%addString(this%str(str_start:str_end))
    end if
  end function

  subroutine stringTypeSetString(this, str)
    class(string_type), intent(inout) :: this
    character(len=*), intent(in) :: str
    this%str = trim(adjustl(str))
    call this%setStringLength() 
    this%is_empty=.false.
  end subroutine

  function stringTypeGetString(this) result(str)
    class(string_type), intent(in) :: this
    character(len=:), allocatable :: str
    if (allocated(str)) deallocate(str)
    allocate( character(len=len(this%str)) :: str )
    str = this%str
  end function

  function getStringLength(this) result(res)
    class(string_type), intent(in) :: this
    integer :: res
    res = this%string_length
  end function
  
  subroutine setStringLength(this)
    class(string_type), intent(inout) :: this
    this%string_length =  len(this%str)
  end subroutine

  ! add string to the list
  subroutine stringArrayAddString(this, str)
    class(string_array), intent(inout) :: this
    character(len=*), intent(in) :: str
    integer, parameter :: allocation_batch = 1
    type(string_type), dimension(:), allocatable :: strings_temp
    integer :: i

    ! check if allocated
    if (allocated(this%strings)) then
      ! if we don't have enough space we need to play the tree card game:
      ! 1 copy the strings in an temporary container:
      ! 2 deallocate the current container and allocate it with more space
      ! 3 copy back the string to the extended container
      if (size(this%strings) == (this%num_string)) then
        allocate(strings_temp(this%num_string))
        do i=1,this%num_string
          strings_temp(i)=this%strings(i)
        end do
        deallocate(this%strings)
        allocate(this%strings(this%num_string+allocation_batch))
        do i=1,this%num_string
          this%strings(i)=strings_temp(i)
        end do
        deallocate(strings_temp)
      end if
    ! if not allocated, allolcate
    else
      allocate(this%strings(allocation_batch))
    end if

    ! increase the counter and store the string
    this%num_string = this%num_string+1
    call this%strings(this%num_string)%setString(str)
  end subroutine

  ! get how many strings we have stored in the string_array
  integer function stringArrayGetNumString(this) result(n_string)
    class(string_array), intent(in) :: this
    n_string = this%num_string
  end function

  ! get a string given the index
  function stringArrayGetStringByIndex(this, idx, expected_string, error) result(str)
    use number_to_string_mod
    class(string_array), intent(in) :: this
    character(:), allocatable :: str
    integer, intent(in) :: idx
    integer, intent(inout), optional :: error
    character(len=*), intent(in), optional :: expected_string
    integer :: i
    if (present(error)) error = 0
    
    if (this%GetNumStrings()==0) then
       if (present(expected_string)) then
        write(*,*) "I was expecting the following string: "
        write(*,*) trim(adjustl(expected_string))
        write(*,*) "but the provided line in the input file is empty"
       endif
       if (present(error)) then
          error = 1
          str = "string is empty"
          return
       else
         error stop "stringArrayGetStringByIndex string is empty"
       endif
    endif
    if (idx > this%GetNumStrings()) then
      if (present(error)) then
        error = 1
        str = "stringArrayGetStringByIndex: index " // int2string(idx) // " out of range!"
        return
      else 
        write(*,*) "stringArrayGetStringByIndex: index ", idx, " out of range. Abort!"
        do i=1, this%GetNumStrings()
          write(*,*) "stringArrayGetStringByIndex, index ", i , " content ", this%strings(idx)%getString()
        end do
        stop
      endif
    end if
    str = this%strings(idx)%getString()
  end function

  ! set a string given the index
  subroutine stringArraySetStringByIndex(this, idx, str)
    class(string_array), intent(inout) :: this
    character(len=*) :: str
    integer, intent(in) :: idx
    integer :: i
    if (this%GetNumStrings()==0) error stop "stringArraySetStringByIndex string is empty"
    if (idx > this%GetNumStrings()) then
      write(*,*) "stringArraySetStringByIndex: index ", idx, " out of range. Abort!"
      do i=1, this%GetNumStrings()
        write(*,*) "stringArraySetStringByIndex, index ", i , " content ", this%strings(idx)%getString()
      end do
      stop
    end if
    call this%strings(idx)%setString(str)
  end subroutine

  subroutine stringArrayReset(this)
    class(string_array), intent(inout) :: this
    if (allocated(this%strings)) deallocate(this%strings)
    this%num_string = 0
  end subroutine

end module string_module
