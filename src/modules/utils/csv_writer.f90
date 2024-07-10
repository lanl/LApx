module csv_writer_mod
  use kinds
  use string_module, only : string_type, string_array

#include "macro_debug.fpp"

  implicit none

  type csv_writer
    integer :: file_id, iostat
    character(len=256) :: error_msg
    type(string_type) :: file_name
    type(string_type) :: smart_string
    type(string_array) :: smart_string_array
    logical :: is_open = .FALSE.
    logical :: new_line = .TRUE.

__DECL_CLASS_UNUSED_THIS__
  contains 
    procedure :: createNewCSV
    procedure :: openToAppendCSV
    procedure, private :: AppendScalarHeader
    procedure, private :: AppendScalarValue
    procedure :: AppendScalar
    procedure, private :: AppendVectorHeader
    procedure, private :: AppendVectorValues
    procedure :: AppendVector
    procedure, private :: AppendTensorHeader
    procedure, private :: AppendTensorValues
    procedure :: AppendTensor
    procedure :: closeCSVFile
    procedure :: checkError
    procedure :: printErrorHeader
    procedure :: getFileID
    procedure :: checkFileIsOpen
    procedure :: errorIfNotOpen
    procedure :: closeFileAndAbortOnError
    procedure :: newLine
  end type

contains


subroutine closeFileAndAbortOnError(this, extra_string)
  use mpi_variables_mod, only : i_am_mpi_master
  implicit none
  class(csv_writer), intent(inout) :: this
  character(len=*), intent(in) :: extra_string

  if (i_am_mpi_master) then
    call this%checkError()
    if (this%iostat.ne.0) then
      call stderr(extra_string)
      call this%printErrorHeader()
      call this%closeCSVFile(new_line=.false.)
      error stop "Abort"
    endif
endif
end subroutine

function checkFileIsOpen(this) result(is_open)
  use mpi_variables_mod, only : i_am_mpi_master
  implicit none
  class(csv_writer), intent(in) :: this
  logical :: is_open
  if (i_am_mpi_master) then
    is_open = this%is_open
  endif
end function

subroutine errorIfNotOpen(this,extra_string)
  use mpi_variables_mod, only : i_am_mpi_master
  implicit none
  class(csv_writer), intent(in) :: this
  character(len=*), intent(in) :: extra_string

  if (i_am_mpi_master) then
    if (.not.(this%checkFileIsOpen())) then
      call stderr(trim(adjustl(extra_string)))
      call stderr("errorIfNotOpen: The file named "// this%file_name%getString()// "is not open")
      error stop "errorIfNotOpen Abort"
    endif
  endif
end subroutine

function getFileID(this) result(file_id)
  use mpi_variables_mod, only : i_am_mpi_master
  implicit none
  class(csv_writer), intent(in) :: this
  integer :: file_id
  if (i_am_mpi_master) then
    file_id = this%file_id
  endif
end function

subroutine createNewCSV(this, filename)
  use mpi_variables_mod, only : i_am_mpi_master
  implicit none
  class(csv_writer), intent(inout) :: this
  character(len=*), intent(in) :: filename
  logical :: csv_exists

  if (i_am_mpi_master) then

    this%iostat = 0
    call this%file_name%setString(filename)
    inquire(file=filename, exist=csv_exists)
    if (csv_exists) then
      open( newunit=this%file_id , file = filename , status='replace' , action='write', iostat=this%iostat, iomsg=this%error_msg)
    else 
      open( newunit=this%file_id , file = filename , status='new' , action='write', iostat=this%iostat, iomsg=this%error_msg)
    endif
    call this%closeFileAndAbortOnError("createNewCSV Abort")
    this%is_open=.TRUE.
    call this%closeCSVFile(new_line=.false.)
  endif
end subroutine


subroutine openToAppendCSV(this)
  use mpi_variables_mod, only : i_am_mpi_master
  implicit none
  class(csv_writer), intent(inout) :: this
  if (i_am_mpi_master) then
    if (this%is_open) then
      call stderr("openToAppendCSV: Trying to open "// this%file_name%getString()// " which is already open")
      error stop "Abort!"
    endif
    open( newunit=this%file_id , file = this%file_name%getString() , status='old' , action='write', position="append", iostat=this%iostat, iomsg=this%error_msg)
    call this%closeFileAndAbortOnError("openToAppendCSV Abort")
    this%is_open=.TRUE.
  endif
end subroutine

subroutine closeCSVFile(this, new_line)
  use mpi_variables_mod, only : i_am_mpi_master
  implicit none
  class(csv_writer), intent(inout) :: this
  logical, optional, intent(in) :: new_line
  logical :: new_line_
  new_line_ = .TRUE.
  if (present(new_line)) new_line_=new_line
  if (i_am_mpi_master) then
    IF (new_line_) call this%NewLine()
    call this%errorIfNotOpen("closeCSVFile")
    close(this%file_id, iostat=this%iostat, iomsg=this%error_msg)
    this%is_open=.FALSE.
  endif
end subroutine

subroutine checkError(this)
  use number_to_string_mod, only : int2string
  use mpi_variables_mod, only : i_am_mpi_master
  implicit none
  class(csv_writer), intent(in) :: this
  if (i_am_mpi_master) then
    if (this%iostat.ne.0) then
      call this%printErrorHeader()
      call stderr("error in file " // this%file_name%getString())
      call stderr("IOstate error detected on file "// this%file_name%getString())
      call stderr("IOstate code " // int2string(this%iostat))
      call stderr("Error message " // this%error_msg)
    endif
  endif
end subroutine

subroutine printErrorHeader(this)
  use mpi_variables_mod, only : i_am_mpi_master
  implicit none
  class(csv_writer), intent(in) :: this
  if (i_am_mpi_master) call stderr("error in file " // this%file_name%getString())
end subroutine

subroutine AppendScalarHeader(this, scalar_name)
  use mpi_variables_mod, only : i_am_mpi_master
  implicit none
  class(csv_writer), intent(inout) :: this
  character(len=*), intent(in) :: scalar_name
  if (i_am_mpi_master) then
    call this%errorIfNotOpen("AppendScalarHeader")

    if (.not.(this%new_line))  write(this%file_id, fmt='(*(1A))', advance='no', iostat=this%iostat, iomsg = this%error_msg) ","
    write(this%file_id,fmt='(*(A))', advance='no',  iostat=this%iostat, iomsg=this%error_msg) scalar_name
  
    call this%closeFileAndAbortOnError("AppendScalarHeader")
    this%new_line=.FALSE.
  endif
end subroutine

subroutine AppendScalarValue(this, scalar_value)
  use mpi_variables_mod, only : i_am_mpi_master
  implicit none
  class(csv_writer), intent(inout) :: this
  real(k_real), intent(in) :: scalar_value
  
  if (i_am_mpi_master) then
    call this%errorIfNotOpen("AppendScalarValue")

    if (.not.(this%new_line))  write(this%file_id, fmt='(*(1A))', advance='no', iostat=this%iostat, iomsg=this%error_msg) ","
    write(this%file_id,fmt='(*(G0))', advance='no',  iostat=this%iostat, iomsg=this%error_msg) scalar_value
    call this%closeFileAndAbortOnError("AppendScalarValue")
    this%new_line=.FALSE.
  endif
end subroutine

subroutine AppendScalar(this, scalar_value, scalar_name, write_header)
  implicit none
  class(csv_writer), intent(inout) :: this
  real(k_real), intent(in) :: scalar_value
  character(len=*), intent(in) :: scalar_name
  logical, intent(in) :: write_header
  
  if (write_header) then
    call this%AppendScalarHeader(scalar_name)
  else
    call this%AppendScalarValue(scalar_value)
  endif

end subroutine

subroutine AppendVectorHeader(this, n_components, vector_name)
  use write_to_file_utils_mod, only : writeTensorIndecesToFile
  use mpi_variables_mod, only : i_am_mpi_master
  use number_to_string_mod, only : int2string
  implicit none
  class(csv_writer), intent(inout) :: this
  character(len=*), intent(in) :: vector_name
  integer, intent(in) :: n_components
  integer :: idx
  
  if (i_am_mpi_master) then
    call this%errorIfNotOpen("AppendVectorHeader")
    
    do idx=1,n_components
      call this%smart_string%setString(trim(adjustl(vector_name))//trim(adjustl(int2string(idx))))
      call this%AppendScalarHeader(this%smart_string%getString())
      this%new_line=.FALSE.
    enddo
  endif
end subroutine

subroutine AppendVectorValues(this, n_components, vector_values)
  use write_to_file_utils_mod, only : writeTensorIndecesToFile
  use mpi_variables_mod, only : i_am_mpi_master
  implicit none
  class(csv_writer), intent(inout) :: this
  integer, intent(in) :: n_components
  real(k_real), intent(in) :: vector_values(n_components)
    integer :: idx
  
  if (i_am_mpi_master) then
    call this%errorIfNotOpen("AppendVectorValues")
    do idx=1,n_components
      call this%AppendScalarValue(vector_values(idx))
      this%new_line=.FALSE.
    enddo
  endif
end subroutine

subroutine AppendVector(this, n_components, vector_values, vector_name, write_header)
  use write_to_file_utils_mod, only : write2DArrayToFileInline
  implicit none
  class(csv_writer), intent(inout) :: this
  integer, intent(in) :: n_components
  real(k_real), dimension(:), intent(in) :: vector_values
  character(len=*), intent(in) :: vector_name
  logical, intent(in) :: write_header

  if (write_header) then
    call this%AppendVectorHeader( n_components, vector_name)
  else
    call this%AppendVectorValues( n_components, vector_values)
  endif
end subroutine

subroutine NewLine(this)
  use mpi_variables_mod, only : i_am_mpi_master
  implicit none
  class(csv_writer), intent(inout) :: this
  if (i_am_mpi_master) then
    call this%errorIfNotOpen("AppendScalarHeader")
    write(this%file_id,fmt='(*(A))', advance='yes',  iostat=this%iostat, iomsg=this%error_msg) ""
    call this%closeFileAndAbortOnError("AppendScalarHeader")
    this%new_line=.TRUE.
  endif

end subroutine

subroutine AppendTensorHeader(this, tensor_name)
  use write_to_file_utils_mod, only : writeTensorIndecesToFile
  use mpi_variables_mod, only : i_am_mpi_master
  implicit none
  class(csv_writer), intent(inout) :: this
  character(len=*), intent(in) :: tensor_name

  if (i_am_mpi_master) then
    call this%errorIfNotOpen("AppendTensorHeader")
    call writeTensorIndecesToFile(this%file_id, tensor_name, new_line = this%new_line, my_stat_in = this%iostat, my_err_msg_in=this%error_msg)  
    call this%closeFileAndAbortOnError("AppendTensorHeader")
    this%new_line=.FALSE.
  endif
end subroutine

subroutine AppendTensorValues(this, tensor, stress_strain)
  use write_to_file_utils_mod, only : write2DArrayToFileInline
  use mpi_variables_mod, only : i_am_mpi_master
  implicit none
  class(csv_writer), intent(inout) :: this
  real(k_real), dimension(:,:), intent(in) :: tensor
  character(len=*), intent(in) :: stress_strain
  if (i_am_mpi_master) then
    call this%errorIfNotOpen("AppendTensorValues")
    call write2DArrayToFileInline(this%file_id, tensor, stress_strain, new_line = this%new_line, my_stat_in = this%iostat, my_err_msg_in=this%error_msg)
    call this%closeFileAndAbortOnError("AppendTensorValues")
    this%new_line=.FALSE.
  endif
end subroutine

subroutine AppendTensor(this, tensor, stress_strain, tensor_name, write_header)
  use write_to_file_utils_mod, only : write2DArrayToFileInline
  implicit none
  class(csv_writer), intent(inout) :: this
  real(k_real), dimension(:,:), intent(in) :: tensor
  character(len=*), intent(in) :: stress_strain
  character(len=*), intent(in) :: tensor_name
  logical, intent(in) :: write_header

  if (write_header) then
    call this%AppendTensorHeader(tensor_name)
  else
    call this%AppendTensorValues( tensor, stress_strain)
  endif
end subroutine

end module
