module phase_info_mod
    use kinds
    use string_module, only : string_type, string_array
    use all_grid_data_mod , only : all_grid_data
    use phase_material_mod, only : phase_material_type
    use polymorphic_dtype_array_mod, only : dtype_array_ptr
    implicit none

type phase_info_type
    logical :: is_gas_phase = .false.

    logical ::   use_isotropic_plasticity=.false., &
                 use_crystal_plasticity=.false., &
                 use_glide=.false., &
                 use_climb=.false., &
                 use_diffusion= .false., &
                 use_porosity=.false., &
                 use_damage=.false.

    integer :: num_slip_system = 0
    integer :: phase_id = 0
    logical :: is_initialized = .false.
    real(k_real), dimension(:), pointer :: burger_vector_length => null(), &
                                           single_crystal_direction => null(), &
                                           single_crystal_normal => null()
    real(k_real), dimension(:,:), pointer :: single_crystal_schmid => null(), &
                                             single_crystal_climb => null()

    type(phase_material_type), pointer :: phase_material_ptr => null()
    type(all_grid_data), pointer :: phase_grid_data => null()

    contains
    procedure :: init
    procedure :: getNumSS
    procedure :: setNumSS
    procedure :: checkInit
    
end type

type, extends(dtype_array_ptr) :: phase_info_array_type
contains
  procedure :: addElement => addPhaseInfoElem
  procedure :: getElementPtr => getPhaseInfoPtr
  procedure :: getNumSSPerPhase => getNumSSPhaseInfo
end type


contains
    subroutine init(this, phase_id, is_gas_phase, phase_grid_data_ptr, phase_material_ptr)
      implicit none
      class(phase_info_type), intent(inout) :: this
      integer, intent(in) :: phase_id
      logical, intent(in) :: is_gas_phase
      class(all_grid_data), intent(in), pointer :: phase_grid_data_ptr
      class(phase_material_type), intent(in), pointer :: phase_material_ptr


      this%phase_id = phase_id
      this%is_gas_phase = is_gas_phase
      this%phase_material_ptr => phase_material_ptr
      this%phase_grid_data => phase_grid_data_ptr
      this%is_initialized = .true.

    end subroutine

    subroutine setNumSS(this, n_ss)
      implicit none
      class(phase_info_type), intent(inout) :: this
      integer, intent(in) :: n_ss
      this%num_slip_system= n_ss
      this%use_crystal_plasticity=.true.
    end subroutine 

    function getNumSS(this) result(n_ss)
      implicit none
      class(phase_info_type), intent(in) :: this
      integer :: n_ss
      call this%checkInit()

      n_ss = this%num_slip_system
    end function 

    subroutine checkInit(this)
      implicit none
      class(phase_info_type), intent(in) :: this
      if (.not.this%is_initialized) then
        error stop "phase info not initialized" 
      endif
    end subroutine

subroutine readPhaseFile(this, matf_reader)
  use read_from_file_utils, only : file_reader
  use number_to_string_mod, only : int2string
  implicit none
  class(phase_info_type), intent(inout) :: this
  class(file_reader), intent(inout) :: matf_reader
  ! type(string_type) :: read_string_buffer
  ! type(string_array) :: dummy_string_array
  ! real(k_real), pointer, dimension(:) :: real_pointer_read_buffer =>null()
  ! logical :: logical_buffer

  call matf_reader%readParameter("use-isotropic-plasticity", this%use_isotropic_plasticity)
  call matf_reader%readParameter("use-crystal-plasticity", this%use_crystal_plasticity)
  call matf_reader%readParameter("use-glide", this%use_glide)
  call matf_reader%readParameter("use-climb", this%use_climb)
  call matf_reader%readParameter("use-diffusion", this%use_diffusion)
  call matf_reader%readParameter("use-porosity", this%use_porosity)
  call matf_reader%readParameter("use-damage", this%use_damage)

end subroutine

  !------------------------------------------------------------------------------!
  !                  PHASE INFO ARRAY SUBROUTINE                          !
  !------------------------------------------------------------------------------!
subroutine addPhaseInfoElem(this, new_element)
  class(phase_info_array_type), intent(inout) :: this
  class(phase_info_type), pointer, intent(inout) :: new_element

  call this%extend()
  this%all_pt(this%n_elements)%pt => new_element
  nullify(new_element)
end subroutine

subroutine getPhaseInfoPtr(this, idx, ptr)
  class(phase_info_array_type), intent(in) :: this
  integer, intent(in) :: idx
  class(phase_info_type), pointer, intent(out) :: ptr

  call this%checkElementExist(idx)

  associate(elem => this%all_pt(idx)%pt)
  select type(elem)
    class is (phase_info_type)
      ptr => elem
    class default
      error stop "getPhaseInfoPtr wrong class"
  end select
  end associate
end subroutine

subroutine getNumSSPhaseInfo(this, num_ss)
  class(phase_info_array_type), intent(in) :: this
  integer, dimension(:), pointer, intent(inout) :: num_ss
  integer :: idx
  class(phase_info_type), pointer :: ptr => null()

  if (associated(num_ss)) error stop "getNumSSPerPhase: num_ss already associated"
  allocate(num_ss(this%n_elements))

  do idx = 1, this%n_elements
    call this%getElementPtr(idx, ptr)
    num_ss(idx) = max(ptr%getNumSS(),1)
    nullify(ptr)
  enddo
end subroutine

end module