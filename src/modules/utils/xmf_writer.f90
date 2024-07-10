module xmf_writer_mod
  use kinds
  use string_module, only : string_type, string_array
  use mpi_variables_mod, only : i_am_mpi_master

#include "macro_debug.fpp"

  implicit none

  type xmf_writer_type
    integer :: file_id, iostat
    character(len=256) :: error_msg
    type(string_type) :: file_name
    type(string_type) :: smart_string
    type(string_array) :: smart_string_array
    logical :: is_open = .FALSE.
    integer, dimension(3) :: nx_ny_nz

__DECL_CLASS_UNUSED_THIS__
  contains 
    procedure :: createNewXMF
    procedure :: openToAppendXMF
    procedure :: closeXMFFile
    procedure :: checkError
    procedure :: printErrorHeader
    procedure :: getFileID
    procedure :: checkFileIsOpen
    procedure :: errorIfNotOpen
    procedure :: closeFileAndAbortOnError
    procedure :: writeLine
    procedure :: writeXMFHeader
    procedure :: writeXMFEndOfFile
    procedure :: writeVectorInteger
    procedure :: write3Dscalar
    procedure :: write3DVector
    procedure :: write3DMatrix
  end type

contains


subroutine closeFileAndAbortOnError(this, extra_string)
  
  implicit none
  class(xmf_writer_type), intent(inout) :: this
  character(len=*), intent(in) :: extra_string

  if (i_am_mpi_master) then
    call this%checkError()
    if (this%iostat.ne.0) then
      call stderr(extra_string)
      call this%printErrorHeader()
      error stop "Abort"
    endif
endif
end subroutine

function checkFileIsOpen(this) result(is_open)
  implicit none
  class(xmf_writer_type), intent(in) :: this
  logical :: is_open
  if (i_am_mpi_master) then
    is_open = this%is_open
  endif
end function

subroutine errorIfNotOpen(this,extra_string)
  implicit none
  class(xmf_writer_type), intent(in) :: this
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
  implicit none
  class(xmf_writer_type), intent(in) :: this
  integer :: file_id
  if (i_am_mpi_master) then
    file_id = this%file_id
  endif
end function

subroutine createNewXMF(this, filename, nx_ny_nz)
  implicit none
  class(xmf_writer_type), intent(inout) :: this
  character(len=*), intent(in) :: filename
  integer, intent(in),dimension(3) :: nx_ny_nz
  logical :: csv_exists
  this%nx_ny_nz = nx_ny_nz
  if (i_am_mpi_master) then

    this%iostat = 0
    call this%file_name%setString(filename)
    inquire(file=filename, exist=csv_exists)
    if (csv_exists) then
      open( newunit=this%file_id , file = filename , status='replace' , action='write', iostat=this%iostat, iomsg=this%error_msg)
    else 
      open( newunit=this%file_id , file = filename , status='new' , action='write', iostat=this%iostat, iomsg=this%error_msg)
    endif
    call this%closeFileAndAbortOnError("createNewXMF Abort")
    this%is_open=.TRUE.
  endif
end subroutine


subroutine openToAppendXMF(this)
  implicit none
  class(xmf_writer_type), intent(inout) :: this
  if (i_am_mpi_master) then
    if (this%is_open) then
      call stderr("openToAppendXMF: Trying to open "// this%file_name%getString()// " which is already open")
      error stop "Abort!"
    endif
    open( newunit=this%file_id , file = this%file_name%getString() , status='old' , action='write', position="append", iostat=this%iostat, iomsg=this%error_msg)
    call this%closeFileAndAbortOnError("openToAppendXMF Abort")
    this%is_open=.TRUE.
  endif
end subroutine

subroutine closeXMFFile(this)
  implicit none
  class(xmf_writer_type), intent(inout) :: this

  if (i_am_mpi_master) then
    call this%errorIfNotOpen("closeXMFFile")
    close(this%file_id, iostat=this%iostat, iomsg=this%error_msg)
    this%is_open=.FALSE.
  endif
end subroutine

subroutine checkError(this)
  use number_to_string_mod, only : int2string
  implicit none
  class(xmf_writer_type), intent(in) :: this
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
  implicit none
  class(xmf_writer_type), intent(in) :: this

  if (i_am_mpi_master) call stderr("error in file " // this%file_name%getString())
end subroutine

subroutine writeLine(this,line, new_line)
  implicit none
  class(xmf_writer_type), intent(in) :: this
  character(len=*), intent(in) :: line
  logical, optional :: new_line
  logical :: new_line_

  if (i_am_mpi_master) then
    call this%errorIfNotOpen("")
    new_line_ = .TRUE.
    if (present(new_line)) new_line_ = new_line

    if (new_line_) then
      write(this%file_id,fmt="(A)", advance='yes') line
    else 
      write(this%file_id,fmt="(A)", advance='no') line
    endif
  endif

end subroutine

subroutine writeVectorInteger(this, vector_int, new_line)
  use number_to_string_mod, only : int2string
  implicit none
  class(xmf_writer_type), intent(in) :: this
  integer, intent(in), dimension(:) :: vector_int
  integer :: n_vals, i
  logical, optional :: new_line
  logical :: new_line_

  n_vals = product(shape(vector_int))
  if (i_am_mpi_master) then
    call this%errorIfNotOpen("")
    new_line_ = .TRUE.
    if (present(new_line)) new_line_ = new_line

    do i=1,n_vals
      write(this%file_id,fmt="(A)", advance='no') trim(adjustl(int2string(vector_int(i))))
      if (i < n_vals) write(this%file_id,fmt="(A)", advance='no') " "
    enddo

    if (new_line_) then
      write(this%file_id,fmt="(A)", advance='yes') ""
    endif
  endif

end subroutine

subroutine writeXMFHeader(this, time)
  use number_to_string_mod
  IMPLICIT NONE
  class(xmf_writer_type), intent(in) :: this
  real(k_real), intent(in) :: time
  integer :: i
  if (i_am_mpi_master) then
  call this%errorIfNotOpen("")
  
  call this%writeLine('<?xml version="1.0"?>')
  call this%writeLine('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd"[]>')
  call this%writeLine('<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="2.2">')
  call this%writeLine('<Domain>')
  call this%writeLine('<Grid Name="MyGrid" GridType="Uniform">')
  call this%writeLine('<Topology TopologyType="3DCoRectMesh" Dimensions="', new_line=.FALSE.)
  do i=3,1,-1
  call this%writeLine(trim(adjustl((int2string(this%nx_ny_nz(i)+1))))//" ", new_line=.FALSE.)
  enddo
  call this%writeLine('"></Topology>')
  call this%writeLine( '<Geometry Type="ORIGIN_DXDYDZ">')
  call this%writeLine( '<!-- Origin -->')
  call this%writeLine( '<DataItem Format="XML" Dimensions="3">0 0 0</DataItem>')
  call this%writeLine( '<!-- DxDyDz (Spacing/Resolution)-->')
  call this%writeLine( '<DataItem Format="XML" Dimensions="3">1 1 1</DataItem>')
  call this%writeLine( '</Geometry>')
  call this%writeLine( '<Time Value="'//trim(adjustl(real2string(time)))//'"/>')
  endif
end subroutine

subroutine writeXMFEndOfFile(this)
  use number_to_string_mod
  IMPLICIT NONE
  class(xmf_writer_type), intent(in) :: this
  if (i_am_mpi_master) then
  call this%errorIfNotOpen("")
  
  call this%writeLine('</Grid>')
  call this%writeLine('</Domain>')
  call this%writeLine('</Xdmf>')
  endif
end subroutine

subroutine write3Dscalar(this, hdf5_fname, field_full_path, att_name, n_stateful_level, suffix)
  use number_to_string_mod
  IMPLICIT NONE
  class(xmf_writer_type), intent(in) :: this
  character(len=*), intent(in) :: hdf5_fname, att_name, field_full_path
  character(len=*), intent(in), optional :: suffix
  integer, intent(in) :: n_stateful_level

  if (i_am_mpi_master) then
  call this%writeLine( '<Attribute Name="'//trim(adjustl(att_name)), new_line=.FALSE. )
  if (present(suffix)) call this%writeLine(trim(adjustl(suffix)), new_line=.FALSE. )
  call this%writeLine( '" AttributeType="Scalar" Center="Cell">')
  call this%writeLine( '<DataItem ItemType="HyperSlab"')
  call this%writeLine( ' Dimensions="', new_line=.FALSE. )
  call this%writeVectorInteger((/this%nx_ny_nz(3), this%nx_ny_nz(2), this%nx_ny_nz(1)/) , new_line=.FALSE. )
  call this%writeLine('"')
  call this%writeLine( 'Type="HyperSlab">')

  call this%writeLine( '<DataItem Dimensions="3 4"')
  call this%writeLine( 'Format="XML">')
  call this%writeVectorInteger( (/0, 0, 0, 0/) )
  call this%writeVectorInteger( (/n_stateful_level, 1, 1, 1/))
  call this%writeVectorInteger(  (/1, this%nx_ny_nz(3), this%nx_ny_nz(2), this%nx_ny_nz(1)/))
  call this%writeLine(  '</DataItem>' )

  call this%writeLine(  '<DataItem')
  call this%writeLine(  'Name="'//trim(adjustl(att_name)), new_line=.FALSE. )
  if (present(suffix)) call this%writeLine(trim(adjustl(suffix)), new_line=.FALSE. )
  call this%writeLine(  '"')
  call this%writeLine(  'Dimensions="', new_line=.FALSE. )
  call this%writeVectorInteger((/n_stateful_level,  this%nx_ny_nz(3), this%nx_ny_nz(2), this%nx_ny_nz(1)/), new_line=.FALSE.)
  call this%writeLine('"')
  call this%writeLine(  'Format="HDF">')
  call this%writeLine(  trim(adjustl(hdf5_fname))//":"//trim(adjustl(field_full_path))//trim(adjustl(att_name)))
  call this%writeLine(  '</DataItem>')
  call this%writeLine(  '</DataItem>')
  call this%writeLine(  '</Attribute>')
  endif
end subroutine

subroutine write3DVector(this, hdf5_fname, field_full_path, att_name, n_components, n_stateful_level, suffix)
  use number_to_string_mod
  IMPLICIT NONE
  class(xmf_writer_type), intent(in) :: this
  character(len=*), intent(in) :: hdf5_fname, att_name, field_full_path
  integer, intent(in) :: n_components
  character(len=*), intent(in), optional :: suffix
  integer, intent(in) :: n_stateful_level
  integer :: i
  if (i_am_mpi_master) then
    do i = 1,n_components
  call this%writeLine( '<Attribute Name="'//trim(adjustl(att_name))//'_i'//trim(adjustl(int2string(i))), new_line=.FALSE. )
  if (present(suffix)) call this%writeLine(trim(adjustl(suffix)), new_line=.FALSE. )
  call this%writeLine( '" AttributeType="Scalar" Center="Cell">')
  call this%writeLine( '<DataItem ItemType="HyperSlab"')
  call this%writeLine( ' Dimensions="', new_line=.FALSE. )
  call this%writeVectorInteger((/this%nx_ny_nz(3), this%nx_ny_nz(2), this%nx_ny_nz(1)/) , new_line=.FALSE. )
  call this%writeLine('"')
  call this%writeLine( 'Type="HyperSlab">')

  call this%writeLine( '<DataItem Dimensions="3 5"')
  call this%writeLine( 'Format="XML">')
  call this%writeVectorInteger( (/0, 0, 0, 0, i-1/) )
  call this%writeVectorInteger( (/n_stateful_level, 1, 1, 1, n_components/))
  call this%writeVectorInteger(  (/1, this%nx_ny_nz(3), this%nx_ny_nz(2), this%nx_ny_nz(1), 1/))
  call this%writeLine(  '</DataItem>' )

  call this%writeLine(  '<DataItem')
  call this%writeLine(  'Name="'//trim(adjustl(att_name)), new_line=.FALSE. )
  if (present(suffix)) call this%writeLine(trim(adjustl(suffix)), new_line=.FALSE. )
  call this%writeLine(  '"')
  call this%writeLine(  'Dimensions="', new_line=.FALSE. )
  call this%writeVectorInteger((/n_stateful_level,  this%nx_ny_nz(3), this%nx_ny_nz(2), this%nx_ny_nz(1), n_components/), new_line=.FALSE.)
  call this%writeLine('"')
  call this%writeLine(  'Format="HDF">')
  call this%writeLine(  trim(adjustl(hdf5_fname))//":"//trim(adjustl(field_full_path))//trim(adjustl(att_name)))
  call this%writeLine(  '</DataItem>')
  call this%writeLine(  '</DataItem>')
  call this%writeLine(  '</Attribute>')
    enddo
  endif
end subroutine

subroutine write3DMatrix(this, hdf5_fname, field_full_path, att_name, n_comp_i, n_comp_j, n_stateful_level, suffix)
  use number_to_string_mod
  IMPLICIT NONE
  class(xmf_writer_type), intent(in) :: this
  character(len=*), intent(in) :: hdf5_fname, att_name, field_full_path
  integer, intent(in) :: n_comp_i, n_comp_j
  character(len=*), intent(in), optional :: suffix
  integer, intent(in) :: n_stateful_level
  integer :: i, j
  if (i_am_mpi_master) then
  do i = 1,n_comp_i
    do j = 1,n_comp_j
  call this%writeLine( '<Attribute Name="'//trim(adjustl(att_name))//'_i'//trim(adjustl(int2string(i)))//'_j'//trim(adjustl(int2string(j))), new_line=.FALSE. )
  if (present(suffix)) call this%writeLine(trim(adjustl(suffix)), new_line=.FALSE. )
  call this%writeLine( '" AttributeType="Scalar" Center="Cell">')
  call this%writeLine( '<DataItem ItemType="HyperSlab"')
  call this%writeLine( ' Dimensions="', new_line=.FALSE. )
  call this%writeVectorInteger((/this%nx_ny_nz(3), this%nx_ny_nz(2), this%nx_ny_nz(1)/) , new_line=.FALSE. )
  call this%writeLine('"')
  call this%writeLine( 'Type="HyperSlab">')

  call this%writeLine( '<DataItem Dimensions="3 6"')
  call this%writeLine( 'Format="XML">')
  call this%writeVectorInteger( (/0, 0, 0, 0, j-1, i-1/) )
  call this%writeVectorInteger( (/n_stateful_level, 1, 1, 1, n_comp_j, n_comp_i/))
  call this%writeVectorInteger(  (/1, this%nx_ny_nz(3), this%nx_ny_nz(2), this%nx_ny_nz(1), 1, 1/))
  call this%writeLine(  '</DataItem>' )

  call this%writeLine(  '<DataItem')
  call this%writeLine(  'Name="'//trim(adjustl(att_name)), new_line=.FALSE. )
  if (present(suffix)) call this%writeLine(trim(adjustl(suffix)), new_line=.FALSE. )
  call this%writeLine(  '"')
  call this%writeLine(  'Dimensions="', new_line=.FALSE. )
  call this%writeVectorInteger((/n_stateful_level,  this%nx_ny_nz(3), this%nx_ny_nz(2), this%nx_ny_nz(1), n_comp_j, n_comp_i/), new_line=.FALSE.)
  call this%writeLine('"')
  call this%writeLine(  'Format="HDF">')
  call this%writeLine(  trim(adjustl(hdf5_fname))//":"//trim(adjustl(field_full_path))//trim(adjustl(att_name)))
  call this%writeLine(  '</DataItem>')
  call this%writeLine(  '</DataItem>')
  call this%writeLine(  '</Attribute>')
  enddo
  enddo
  endif
end subroutine

end module
