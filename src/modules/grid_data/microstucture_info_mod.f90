module microstructure_info_mod
    use kinds
    use string_module, only : string_type, string_array
    implicit none

type microstructure_info_type
    integer :: num_phases, n_points
    integer, pointer, dimension(:) :: num_voxel
    real(k_real), dimension(:), pointer :: voxel_physical_dimension => null()
    real(k_real) :: voxel_volume
    type(string_type) :: microstucture_fname
    type(string_array) :: single_crystal_file_name
    logical, dimension(:), pointer :: is_gas_phase => null()

    contains
    
    procedure :: readMicorstructureInfoFromOptionFile
    procedure :: getNumPhases
    procedure :: getPhaseFileNameByIndex
    procedure :: getIsGasPhaseByIndex
    procedure :: getIsGasPhaseArrayPtr
    procedure :: getGlobalDimensions
    procedure :: getMicrostructureFileName
    procedure :: getVoxelPhysicalDimension
end type

contains
subroutine readMicorstructureInfoFromOptionFile(this, optf_reader)
  use read_from_file_utils, only : file_reader
  use number_to_string_mod, only : int2string
  implicit none
  class(microstructure_info_type), intent(inout) :: this
  class(file_reader), intent(inout) :: optf_reader
  type(string_type) :: read_string_buffer
  type(string_array) :: dummy_string_array
  logical :: logical_buffer
  integer :: iph
  call optf_reader%readVectorParameter("-voxel-physical-dimensions[Dx,Dy,Dz[m]]", 3, this%voxel_physical_dimension)
  call optf_reader%readParameter("-microstructure-file-name", this%microstucture_fname)
  this%voxel_volume = product(this%voxel_physical_dimension)


  call optf_reader%readVectorParameter("-microstructure-voxels[nx,ny,nz]", 3, this%num_voxel)
  this%n_points = product(this%num_voxel)
  call optf_reader%readParameter("-number-of-phases", this%num_phases)
  if (this%num_phases<=0) error stop "number-of-phases must be a postivie integer. abort"

  ! allocate per phase required variables
  allocate(this%is_gas_phase(this%num_phases))
  call read_string_buffer%resetString()

  do iph=1,this%num_phases
    call optf_reader%readLineAndCheckStringsAreEqual("-phase-"//trim(adjustl(int2string(iph)))//"-info", dummy_string_array)
    call optf_reader%readParameter("gas-phase[TRUE/FALSE]", logical_buffer)
    this%is_gas_phase(iph) = logical_buffer
    if(.not.(this%is_gas_phase(iph))) then
      call optf_reader%readParameter("single-crystal-phase-file", read_string_buffer)
      call this%single_crystal_file_name%addString(read_string_buffer%getString())
      call read_string_buffer%resetString()
    endif
  enddo
end subroutine

function getVoxelPhysicalDimension(this) result(voxel_size)
  implicit none
  class(microstructure_info_type), intent(in) :: this
  real(k_real), dimension(3) :: voxel_size
  voxel_size = this%voxel_physical_dimension
end function

function getNumPhases(this) result(num_phases)
  implicit none
  class(microstructure_info_type), intent(in) :: this
  integer :: num_phases
  num_phases = this%num_phases
end function

function getIsGasPhaseByIndex(this,ph_idx) result(is_gas)
  implicit none
  class(microstructure_info_type), intent(in) :: this
  integer, intent(in) :: ph_idx
  logical :: is_gas
  is_gas = this%is_gas_phase(ph_idx)
end function

subroutine getIsGasPhaseArrayPtr(this, is_gas_array_ptr)
  implicit none
  class(microstructure_info_type), intent(in) :: this
  logical, dimension(:), pointer,  intent(inout) :: is_gas_array_ptr
  if (associated(is_gas_array_ptr))  error stop "getIsGasPhaseArrayPtr: is_gas_array_ptr already associated"
  is_gas_array_ptr => this%is_gas_phase

end subroutine

function getPhaseFileNameByIndex(this, ph_idx) result(fname)
  implicit none
  class(microstructure_info_type), intent(in) :: this
  integer, intent(in) :: ph_idx
  character(:), allocatable :: fname
  integer :: error

  fname = this%single_crystal_file_name%getStringByIndex(ph_idx, error=error)
  if (error>0) error stop "MicrostructureInfo:getPhaseFileNameByIndex string is empty"
end function

function getMicrostructureFileName(this) result(fname)
  implicit none
  class(microstructure_info_type), intent(in) :: this
  character(:), allocatable :: fname

  fname = this%microstucture_fname%getString()
end function

subroutine getGlobalDimensions(this, n_voxels_xyz)
  implicit none
  class(microstructure_info_type), intent(in) :: this
  integer, intent(out) :: n_voxels_xyz(3)

  n_voxels_xyz = this%num_voxel
end subroutine

end module