module isotropic_diffusion_mod
use kinds
use diffusion_base_mod, only : diffusion_base
use all_grid_data_mod, only : all_grid_data
use common_material_parameter_mod, only : common_material_parameter

implicit none

type, extends(diffusion_base) :: nabarro_herring_plus_coble_diffusion
  real(k_real), pointer, dimension(:) :: geometric_factor => null() ! a pointer to the gb diffusivity coefficients
  integer, pointer, dimension(:,:,:) :: gb_type_id_grid => null()
  integer, pointer :: gb_type_id_ptr=> null()
  logical :: use_RuleOfMixture_Diffusion = .TRUE.
contains
  procedure :: readMaterialParametersFromFile => readMaterialParametersFromFileNabarroCobleDiffusion
  procedure :: initParameters => initParametersNabarroCobleDiffusion
  procedure :: addFieldVariablesToGrid => addFieldVariablesToGridNabarroCobleDiffusion
  procedure :: initGridPointers => initGridPointersNabarroCobleDiffusion
  procedure :: setPointData => setPointDataNabarroCobleDiffusion
  procedure :: initStateVariablesAtMaterialPoint => initStateVariablesAtMaterialPointNabarroCobleDiffusion
  procedure :: updateStateVariablesAtMaterialPointInnerLoop => updateStateVarsAtMaterialPointInnerLoopNabarroCobleDiffusion
  procedure :: updateStateVariablesAtMaterialPointStaggered => updateStateVarsAtMaterialPointStaggeredNabarroCobleDiffusion
  procedure :: computeEpsilonDotAndDepsilonDotDStress => computeEpsilonDotAndDepsilonDotDStressNabarroCobleDiffusion

  procedure :: computeDiffusivity
  procedure :: computeVacancyThermalEquilibriumConcentration
  procedure :: acceptRejectSolution => acceptRejectSolutionNabarroCobleDiffusion
end type

interface
  module subroutine readMaterialParametersFromFileNabarroCobleDiffusion(this, matf_reader)
    use read_from_file_utils, only : file_reader
    implicit none
    class(nabarro_herring_plus_coble_diffusion), intent(inout) :: this
    class(file_reader), intent(inout) :: matf_reader
  end subroutine

  module subroutine initParametersNabarroCobleDiffusion(this, phase_id, common_material_parameter_ptr, use_damage)
    implicit none
    class(nabarro_herring_plus_coble_diffusion), intent(inout) :: this
    integer, intent(in) :: phase_id
    class(common_material_parameter), pointer, intent(in) :: common_material_parameter_ptr
    logical, intent(in) :: use_damage
  end subroutine

  module subroutine  addFieldVariablesToGridNabarroCobleDiffusion(this)
    use grid_data_var_type_mod
    implicit none
    class(nabarro_herring_plus_coble_diffusion), intent(inout) :: this
  end subroutine

  module subroutine initGridPointersNabarroCobleDiffusion(this)
    implicit none
    class(nabarro_herring_plus_coble_diffusion), intent(inout) :: this
  end subroutine

  module subroutine setPointDataNabarroCobleDiffusion(this, ix, iy, iz)
    implicit none
    class(nabarro_herring_plus_coble_diffusion), intent(inout) :: this
    integer, intent(in) ::  ix, iy, iz
  end subroutine

  module subroutine computeEpsilonDotAndDepsilonDotDStressNabarroCobleDiffusion(this, stress6, epsilon_dot, depsilon_dot_dstress)
    implicit none
    class(nabarro_herring_plus_coble_diffusion), intent(inout) :: this
    real(k_real), target, intent(in) ::stress6(6)
    real(k_real), target, intent(out) :: epsilon_dot(6), depsilon_dot_dstress(6,6)
  end subroutine

  module subroutine initStateVariablesAtMaterialPointNabarroCobleDiffusion(this)
    use test_utils_mod
    implicit none
    class(nabarro_herring_plus_coble_diffusion), intent(inout) :: this
  end subroutine

  module subroutine updateStateVarsAtMaterialPointInnerLoopNabarroCobleDiffusion(this)
    implicit none
    class(nabarro_herring_plus_coble_diffusion), intent(inout) :: this
  end subroutine

  module subroutine updateStateVarsAtMaterialPointStaggeredNabarroCobleDiffusion(this)
    implicit none
    class(nabarro_herring_plus_coble_diffusion), intent(inout) :: this
  end subroutine

  module subroutine computeDiffusivity(this, gb_id, D)
    implicit none
    class(nabarro_herring_plus_coble_diffusion), intent(inout) :: this
    integer, intent(in) :: gb_id
    real(k_real), intent(out) :: D
  end subroutine

  module subroutine computeVacancyThermalEquilibriumConcentration(this, gb_id, C)
    implicit none
    class(nabarro_herring_plus_coble_diffusion), intent(inout) :: this
    integer, intent(in) :: gb_id
    real(k_real), intent(out) :: C
  end subroutine

  module subroutine acceptRejectSolutionNabarroCobleDiffusion(this, dt_max, accept_solution_flag)
    implicit none
    class(nabarro_herring_plus_coble_diffusion), intent(inout) :: this
    real(k_real), intent(out) :: dt_max
    logical, intent(out) :: accept_solution_flag
  end subroutine

end interface

contains 
subroutine readMaterialParametersFromFileISOdiffusion(matf_reader, phase_id, use_damage, &
                                                      all_mighty_grid_in, sim_all_macro_data, the_bc_object, &
                                                        common_material_parameter_ptr, inelastic_strain_base_ptr)
  use read_from_file_utils, only : file_reader
  use all_mighty_grid_mod, only : all_mighty_grid_type
  use simumaltion_macro_data_mod, only : sim_all_macro_data_obj
  use bc_objects_mod, only : boundary_condition_array_type
  use common_material_parameter_mod, only : common_material_parameter
  use inelastic_strain_mod, only : inelastic_strain_base
  use string_module, only : string_type
  use log_file_mod, only : write_detailed_log_to_screen
  ! use macro
  implicit none
  type(file_reader), intent(inout) :: matf_reader
  class(inelastic_strain_base), pointer, intent(inout) :: inelastic_strain_base_ptr
  class(common_material_parameter), pointer, intent(inout) :: common_material_parameter_ptr
  type(string_type) :: material_model_name
  integer, intent(in) :: phase_id
  logical, intent(in) :: use_damage
  class(all_mighty_grid_type), intent(in), target :: all_mighty_grid_in
  class(sim_all_macro_data_obj), intent(in), target :: sim_all_macro_data
  class(boundary_condition_array_type), intent(in), target :: the_bc_object
  class(nabarro_herring_plus_coble_diffusion), pointer :: nabarro_herring_plus_coble_diffusion_temp  => null()

  call matf_reader%readParameter("--Diffusion-model", material_model_name)
  select case(material_model_name%getString())
  case ("coble-nabarro-diffusion")
    allocate(nabarro_herring_plus_coble_diffusion_temp)
    inelastic_strain_base_ptr => nabarro_herring_plus_coble_diffusion_temp
    call nabarro_herring_plus_coble_diffusion_temp%linkBaseObjects(phase_id, all_mighty_grid_in, sim_all_macro_data, the_bc_object)
    call nabarro_herring_plus_coble_diffusion_temp%initParameters(phase_id, common_material_parameter_ptr, use_damage)

  case default
    if (write_detailed_log_to_screen) write(*,*) "unrecognized material model for climb", material_model_name%getString()
    if (write_detailed_log_to_screen) write(*,*) "available material models are: "
    if (write_detailed_log_to_screen) write(*,*) "nabarro_herring_plus_coble_diffusion_temp "
    error stop "abort"
  end select

  call inelastic_strain_base_ptr%readMaterialParametersFromFile(matf_reader)

  select case(material_model_name%getString())
  case ("coble-nabarro-diffusion")
    nullify(nabarro_herring_plus_coble_diffusion_temp)
  case default
    if (write_detailed_log_to_screen) write(*,*) "somthing is worng when trying nullifying the pointer for ", material_model_name%getString()
    error stop "abort"
  end select
  
end subroutine
end module isotropic_diffusion_mod
