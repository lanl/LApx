module read_from_file_utils
  use kinds
  use string_module, only : string_type, string_array

  implicit none

  type file_reader
    integer :: file_id, line_number, iostat
    type(string_type) :: file_name
    character(len=1000) :: line
    type(string_type) :: smart_string
    type(string_array) :: smart_string_array
  contains 
    procedure :: openReadTextFile
    procedure :: closeTextFile
    procedure :: checkError
    procedure :: skipEmptyLine
    procedure :: readLine
    procedure :: skipNEmptyLines
    procedure :: printErrorHeader
    procedure :: readLineAndCheckStringsAreEqual
    procedure :: getFileID
    generic, public :: readParameter => readParameterLogical, readParameterReal, readParameterInteger, readParameterString
    procedure, private :: readParameterLogical
    procedure, private :: readParameterReal
    procedure, private :: readParameterInteger
    procedure, private :: readParameterString
    generic, public :: readVectorParameter => readVectorParameterLogical, readVectorParameterReal, readVectorParameterInteger, readVectorParameterString
    procedure, private :: readVectorParameterLogical
    procedure, private :: readVectorParameterReal
    procedure, private :: readVectorParameterInteger
    procedure, private :: readVectorParameterString
    procedure :: readLinearInterpParameter
    procedure :: readMultiLinearInterpParameter
    generic, public :: readStress6FromFile => readStress6FromFileReal, readStress6FromFileInteger
    procedure, private :: readStress6FromFileReal
    procedure, private :: readStress6FromFileInteger
    generic, public :: readDisplacement33FromFile => readDisplacement33FromFileReal, readDisplacement33FromFileInteger
    procedure, private :: readDisplacement33FromFileReal
    procedure, private :: readDisplacement33FromFileInteger
    generic, public :: readScalar => readScalarReal, readScalarInteger
    procedure, private :: readScalarReal
    procedure, private :: readScalarInteger
    generic, public :: readVector => readVectorReal, readVectorInteger
    procedure, private :: readVectorReal
    procedure, private :: readVectorInteger
    generic, public :: readMatrix => readMatrixReal
    procedure, private :: readMatrixReal
  end type

contains

function getFileID(this) result(file_id)
  class(file_reader), intent(inout) :: this
  integer :: file_id
  this%iostat = 0
  file_id = this%file_id

end function

subroutine openReadTextFile(this, filename, file_id)
  class(file_reader), intent(inout) :: this
  character(len=*), intent(in) :: filename
  integer, intent(out), optional :: file_id
  this%iostat = 0
  call this%file_name%setString(filename)
  this%line_number = 0

  open( newunit=this%file_id , file = filename , status='old' , action='read', iostat=this%iostat)
  if (present(file_id)) file_id = this%file_id
  call this%checkError()
end subroutine

subroutine closeTextFile(this)
  class(file_reader), intent(inout) :: this
  close(this%file_id)

  this%iostat = 0
  this%file_id = -1
  call this%file_name%resetString()
  this%line_number = 0
end subroutine

subroutine checkError(this)
  class(file_reader), intent(in) :: this

  if (this%iostat.ne.0) then
    call this%printErrorHeader()
    write(*,*) "IOstate error detected on file ", this%file_name%getString(), "at line number", this%line_number
    write(*,*) "IOstate code ", this%iostat
  endif

end subroutine

subroutine printErrorHeader(this)
  class(file_reader), intent(in) :: this

  write(*,*) "error in file ", this%file_name%getString(), " at line number ", this%line_number

end subroutine

subroutine readLine(this)
  implicit none
  class(file_reader), intent(inout) :: this

  this%line_number = this%line_number+1
  READ(this%file_id,fmt='(A)', iostat=this%iostat) this%line
  call this%checkError()

end subroutine

subroutine skipEmptyLine(this)
  implicit none
  class(file_reader), intent(inout) :: this

  call this%readLine()

  if (this%line.ne."") then
    call this%printErrorHeader()
    write(*,*) "I'm expecting an empty line, instead I'm skipping ``", this%line, "''"
    error stop "ABORT!"
  endif

end subroutine

subroutine skipNEmptyLines(this, n_lines_to_skip)
  implicit none
  class(file_reader), intent(inout) :: this
  integer, intent(in) :: n_lines_to_skip
  integer :: i
  
  if (n_lines_to_skip<1) error stop "skipNLines: n_lines_to_skip<1. Abort!"
  do i=1,n_lines_to_skip
    call this%skipEmptyLine()
  enddo
 end subroutine

subroutine readLineAndCheckStringsAreEqual(this, expected_string, smart_string_array, extra_error_string)
  use string_module, only : string_type
  implicit none
  class(file_reader), intent(inout) :: this
  type(string_array), intent(inout), optional :: smart_string_array
  character(len=*), intent(in) :: expected_string
  character(len=*), intent(in), optional :: extra_error_string
  integer :: smart_string_error
  logical :: string_are_equal

  smart_string_error = 0

  call this%readLine()
  call this%smart_string%resetString()
  call this%smart_string%setString(this%line)

  if (present(smart_string_array)) then
    call smart_string_array%reset()
    smart_string_array = this%smart_string%splitString(" ")
    call this%smart_string%setString(smart_string_array%getStringByIndex(1, expected_string=expected_string, error = smart_string_error))
    if (smart_string_error>0) then
      call this%printErrorHeader()
      error stop "abort"
    endif
  endif

  string_are_equal = this%smart_string%equal(expected_string)
  if(.not.string_are_equal) then
    call this%printErrorHeader()
    write(*,*) "I was expecting the following parameter: ", expected_string
    write(*,*) "instead the first word in the string is: ", this%smart_string%getString()
    if (present(extra_error_string)) write(*,*) extra_error_string
    error stop "Abort!"
  endif
end subroutine


subroutine readParameterLogical(this, parameter_name, param_value)
  implicit none
  class(file_reader), intent(inout) :: this
  character(len=*), intent(in) :: parameter_name
  logical, intent(out) :: param_value


  call this%readLineAndCheckStringsAreEqual(parameter_name, this%smart_string_array)
  if (this%smart_string_array%getNumStrings() < 2) then
    call this%printErrorHeader()
    error stop "there isn't any logical parameter to read. Abort!"
  endif
  call this%smart_string_array%strings(2)%stringToVar(param_value)

end subroutine

subroutine readParameterReal(this, parameter_name, param_value)
  implicit none
  class(file_reader), intent(inout) :: this
  character(len=*), intent(in) :: parameter_name
  real(k_real) :: param_value


  call this%readLineAndCheckStringsAreEqual(parameter_name,this%smart_string_array)
  if (this%smart_string_array%getNumStrings() < 2) then
    call this%printErrorHeader()
    error stop "there isn't any real parameter to read. Abort!"
  endif
  call this%smart_string_array%strings(2)%stringToVar(param_value)

end subroutine

subroutine readParameterInteger(this, parameter_name, param_value)
  implicit none
  class(file_reader), intent(inout) :: this
  character(len=*), intent(in) :: parameter_name
  integer :: param_value


  call this%readLineAndCheckStringsAreEqual(parameter_name, this%smart_string_array)
  if (this%smart_string_array%getNumStrings() < 2) then
    call this%printErrorHeader()
    error stop "there isn't any integer parameter to read. Abort!"
  endif
  call this%smart_string_array%strings(2)%stringToVar(param_value)

end subroutine

subroutine readParameterString(this, parameter_name, param_value)
  use string_module
  implicit none
  class(file_reader), intent(inout) :: this
  character(len=*), intent(in) :: parameter_name
  type(string_type), intent(inout) :: param_value
  integer :: smart_string_error

  smart_string_error = 0
  call this%readLineAndCheckStringsAreEqual(parameter_name, this%smart_string_array)
  if (this%smart_string_array%getNumStrings() < 2) then
    call this%printErrorHeader()
    error stop "there isn't any string parameter to read. Abort!"
  endif

  call param_value%resetString()
  call param_value%setString(this%smart_string_array%getStringByIndex(2,  error = smart_string_error))

    if (smart_string_error>0) then
      call this%printErrorHeader()
      error stop "abort"
    endif

end subroutine

subroutine readVectorParameterLogical(this, parameter_name, n_values, param_values)
  implicit none
  class(file_reader), intent(inout) :: this
  integer, intent(in) :: n_values
  character(len=*), intent(in) :: parameter_name
  logical, pointer, dimension(:), intent(inout) :: param_values
  integer :: i

  if (n_values<1) error stop "readVectorParameterLogical n_values < 1. Abort!"
  if (associated(param_values)) error stop "readVectorParameterLogical param_values already associated. Abort!"
  allocate(param_values(n_values))
  call this%readLineAndCheckStringsAreEqual(parameter_name, this%smart_string_array)
  if (this%smart_string_array%getNumStrings() < 2) then
    call this%printErrorHeader()
    error stop "there isn't any vector-logical parameter to read. Abort!"
  endif
  if (this%smart_string_array%getNumStrings() < 1+n_values) then
    call this%printErrorHeader()
    error stop "there aren't enough vector-logical parameter to read. Abort!"
  endif
  do i = 2,n_values+1
    call this%smart_string_array%strings(i)%stringToVar(param_values(i-1))
  enddo
end subroutine

subroutine readVectorReal(this, n_values, vector)
  implicit none
  class(file_reader), intent(inout) :: this
  integer, intent(in) :: n_values
  real(k_real), dimension(:), intent(inout) :: vector
  integer :: i
  
  call this%readLine()
  call this%smart_string%resetString()
  call this%smart_string%setString(this%line)
  call this%smart_string_array%reset()
  this%smart_string_array = this%smart_string%splitString(" ")

  if (this%smart_string_array%getNumStrings() < n_values) then
    call this%printErrorHeader()
    write(*,*) "write ", this%smart_string%getString()
    error stop "there aren't enough vector-real numbers to read. Abort!"
  endif
  do i = 1,n_values
    call this%smart_string_array%strings(i)%stringToVar(vector(i))
  enddo
end subroutine

subroutine readMatrixReal(this, n_row, n_col, matrix)
  implicit none
  class(file_reader), intent(inout) :: this
  integer, intent(in) :: n_row, n_col
  real(k_real), dimension(:,:), intent(inout) :: matrix
  integer :: row
  
  do row = 1, n_row
    call this%readVectorReal(n_col, matrix(row,:))
  enddo
end subroutine

subroutine readScalarReal(this, scalar)
  implicit none
  class(file_reader), intent(inout) :: this
  real(k_real), intent(inout) :: scalar
  
  call this%readLine()
  call this%smart_string%resetString()
  call this%smart_string%setString(this%line)
  call this%smart_string_array%reset()
  this%smart_string_array = this%smart_string%splitString(" ")

  if (this%smart_string_array%getNumStrings() < 1) then
    call this%printErrorHeader()
    error stop "string is empty. Abort!"
  endif
  call this%smart_string_array%strings(1)%stringToVar(scalar)
end subroutine

subroutine readScalarInteger(this, scalar)
  implicit none
  class(file_reader), intent(inout) :: this
  integer, intent(inout) :: scalar
  
  call this%readLine()
  call this%smart_string%resetString()
  call this%smart_string%setString(this%line)
  call this%smart_string_array%reset()
  this%smart_string_array = this%smart_string%splitString(" ")

  if (this%smart_string_array%getNumStrings() < 1) then
    call this%printErrorHeader()
    error stop "string is empty. Abort!"
  endif
  call this%smart_string_array%strings(1)%stringToVar(scalar)
end subroutine

subroutine readVectorInteger(this, n_values, vector)
  implicit none
  class(file_reader), intent(inout) :: this
  integer, intent(in) :: n_values
  integer, dimension(:), intent(inout) :: vector
  integer :: i
  
  call this%readLine()
  call this%smart_string%resetString()
  call this%smart_string%setString(this%line)
  call this%smart_string_array%reset()
  this%smart_string_array = this%smart_string%splitString(" ")

  if (this%smart_string_array%getNumStrings() < n_values) then
    call this%printErrorHeader()
    error stop "there aren't enough vector-integer numbers to read. Abort!"
  endif
  do i = 1,n_values
    call this%smart_string_array%strings(i)%stringToVar(vector(i))
  enddo
end subroutine

subroutine readVectorParameterReal(this, parameter_name, n_values, param_values)
  implicit none
  class(file_reader), intent(inout) :: this
  integer, intent(in) :: n_values
  character(len=*), intent(in) :: parameter_name
  real(k_real), pointer, dimension(:), intent(inout) :: param_values
  type(string_array) :: smart_string_array
  integer :: i

  if (n_values<1) error stop "readVectorParameterReal n_values < 1. Abort!"
  if (associated(param_values)) error stop "readVectorParameterReal param_values already associated. Abort!"
  allocate(param_values(n_values))
  call this%readLineAndCheckStringsAreEqual(parameter_name, smart_string_array)
  if (smart_string_array%getNumStrings() < 2) then
    call this%printErrorHeader()
    error stop "there isn't any vector-real parameter to read. Abort!"
  endif
  if (smart_string_array%getNumStrings() < 1+n_values) then
    call this%printErrorHeader()
    error stop "there aren't enough vector-real parameter to read. Abort!"
  endif
  do i = 2,n_values+1
    call smart_string_array%strings(i)%stringToVar(param_values(i-1))
  enddo
end subroutine

subroutine readVectorParameterInteger(this, parameter_name, n_values, param_values)
  implicit none
  class(file_reader), intent(inout) :: this
  integer, intent(in) :: n_values
  character(len=*), intent(in) :: parameter_name
  integer, pointer, dimension(:), intent(inout) :: param_values
  type(string_array) :: smart_string_array
  integer :: i

  if (n_values<1) error stop "readVectorParameterInteger n_values < 1. Abort!"
  if (associated(param_values)) error stop "readVectorParameterInteger param_values already associated. Abort!"
  allocate(param_values(n_values))
  call this%readLineAndCheckStringsAreEqual(parameter_name, smart_string_array)
  if (smart_string_array%getNumStrings() < 2) then
    call this%printErrorHeader()
    error stop "there isn't any vector-integer parameter to read. Abort!"
  endif
  if (smart_string_array%getNumStrings() < 1+n_values) then
    call this%printErrorHeader()
    error stop "there aren't enough vector-integer parameter to read. Abort!"
  endif
  do i = 2,n_values+1
    call smart_string_array%strings(i)%stringToVar(param_values(i-1))
  enddo
end subroutine

subroutine readVectorParameterString(this, parameter_name, n_values, param_values)
  use string_module
  implicit none
  class(file_reader), intent(inout) :: this
  integer, intent(in) :: n_values
  character(len=*), intent(in) :: parameter_name
  type(string_array), intent(inout) :: param_values
  type(string_array) :: smart_string_array
  integer :: smart_string_error
  integer :: i

  if (n_values<1) error stop "readVectorParameterString n_values < 1. Abort!"
  call param_values%reset()

  call this%readLineAndCheckStringsAreEqual(parameter_name, smart_string_array)
  if (smart_string_array%getNumStrings() < 2) then
    call this%printErrorHeader()
    error stop "there isn't any vector-string parameter to read. Abort!"
  endif
  if (smart_string_array%getNumStrings() < 1+n_values) then
    call this%printErrorHeader()
    error stop "there aren't enough vector-integer parameter to read. Abort!"
  endif

  do i = 2,n_values+1
    call param_values%addString(smart_string_array%getStringByIndex(i, error= smart_string_error))
    if (smart_string_error>0) then
      call this%printErrorHeader()
      error stop "abort"
    endif
  enddo
end subroutine

subroutine readLinearInterpParameter(this, parameter_name, value_unit, temperature_interpolator)
  use linear_interpolation_mod, only : piecewise_linear_interpolation
  implicit none
  class(piecewise_linear_interpolation), pointer, intent(inout) :: temperature_interpolator
  class(file_reader), intent(inout) :: this
  character(len=*), intent(in) :: parameter_name, value_unit
  integer :: n_values
  real(k_real), dimension(:), pointer :: x_values, y_values

  nullify(x_values, y_values)
  call this%readParameter(parameter_name//"-n-values", n_values)
  call this%readVectorParameter(parameter_name//"-temperature[K]", n_values, x_values)
  call this%readVectorParameter(parameter_name//"-values["//trim(adjustl((value_unit)))//"]", n_values, y_values)

  allocate(temperature_interpolator)
  call temperature_interpolator%init(n_values, x_values, y_values)
  deallocate(x_values, y_values)
  nullify(x_values, y_values)

end subroutine

subroutine readMultiLinearInterpParameter(this, n_interpolators, parameter_name, value_unit, mode_string, multi_linear_interpolator)
  use linear_interpolation_mod, only : piecewise_linear_interpolation_array_type, piecewise_linear_interpolation
  use number_to_string_mod
  use string_module
  implicit none
  class(piecewise_linear_interpolation_array_type), intent(inout) :: multi_linear_interpolator
  class(piecewise_linear_interpolation), pointer :: new_interpolator => null()
  class(file_reader), intent(inout) :: this
  integer, intent(in) :: n_interpolators
  character(len=*), intent(in) :: parameter_name, value_unit, mode_string
  type(string_array) :: dummy_string_array
  integer i

  do i=1, n_interpolators
    call this%readLineAndCheckStringsAreEqual(trim(adjustl(mode_string))//"-"//trim(adjustl(int2string(i))), dummy_string_array)
    call this%readLinearInterpParameter(parameter_name, value_unit, new_interpolator)
    call multi_linear_interpolator%addElement(new_interpolator)
    nullify(new_interpolator)
  enddo
end subroutine

subroutine readStress6FromFileReal(this, stress6Voigt)
  implicit none
  class(file_reader), intent(inout) :: this
  real(k_real), dimension(6), intent(out) :: stress6Voigt
  real(k_real) :: stress6(6)
  
  call this%readVector(3, stress6(1:3))
  call this%readVector(2, stress6(4:5))
  call this%readScalar(stress6(6))

  stress6Voigt(1) = stress6(1)
  stress6Voigt(2) = stress6(4)
  stress6Voigt(3) = stress6(6)
  stress6Voigt(4) = stress6(5)
  stress6Voigt(5) = stress6(3)
  stress6Voigt(6) = stress6(2)
end subroutine

subroutine readStress6FromFileInteger(this, stress6Voigt)
  implicit none
  class(file_reader), intent(inout) :: this
  integer, dimension(6), intent(out) :: stress6Voigt
  integer :: stress6(6)

  call this%readVector(3, stress6(1:3))
  call this%readVector(2, stress6(4:5))
  call this%readScalar(stress6(6))

  stress6Voigt(1) = stress6(1)
  stress6Voigt(2) = stress6(4)
  stress6Voigt(3) = stress6(6)
  stress6Voigt(4) = stress6(5)
  stress6Voigt(5) = stress6(3)
  stress6Voigt(6) = stress6(2)
end subroutine

subroutine readDisplacement33FromFileReal(this, displacement33)
  implicit none
  class(file_reader), intent(inout) :: this
  real(k_real), dimension(3,3), intent(out) :: displacement33
  integer :: i
  do i=1,3
      call this%readVector(3, displacement33(:,i))
      call this%checkError()
      if (this%iostat.ne.0) error stop "Abort"
  enddo
end subroutine

subroutine readDisplacement33FromFileInteger(this, displacement33)
  implicit none
  class(file_reader), intent(inout) :: this
  integer, dimension(3,3), intent(out) :: displacement33
  integer :: i
  do i=1,3
      call this%readVector(3, displacement33(:,i))
      call this%checkError()
      if (this%iostat.ne.0) error stop "Abort"
  enddo
end subroutine

end module
