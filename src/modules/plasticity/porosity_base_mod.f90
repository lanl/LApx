module porosity_base_mod
use kinds
use material_base_mod, only : material_base
use common_material_parameter_mod, only : common_material_parameter

implicit none

type, extends(material_base) :: porosity_base

  real(k_real), pointer, dimension(:,:,:,:,:) :: stress_grid=> null() ! a pointer to the stress vector
  real(k_real), pointer, dimension(:,:,:) :: porosity_grid=> null() ! a pointer to the porosity, only used for damage
  real(k_real), pointer, dimension(:,:,:) :: porosity_grid_old=> null() ! a pointer to the porosity, only used for damage
  real(k_real), pointer, dimension(:,:,:) :: porosity_rate_grid=> null() ! a pointer to the porosity, only used for damage
  real(k_real), pointer, dimension(:,:,:) :: porosity_rate_grid_old=> null() ! a pointer to the porosity, only used for damage
  real(k_real), pointer, dimension(:,:,:) :: porosity_nucleation_grid=> null()
  real(k_real), pointer, dimension(:,:,:) :: porosity_nucleation_grid_old=> null()
  real(k_real), pointer, dimension(:,:,:) :: porosity_growth_grid=> null()
  real(k_real), pointer, dimension(:,:,:) :: porosity_growth_grid_old=> null()
  real(k_real), pointer, dimension(:,:,:) :: effective_porosity_grid=> null() ! a pointer to the effectvie porosity
  real(k_real), pointer, dimension(:,:,:) :: effective_porosity_grid_old=> null() ! a pointer to the effectvie porosity
  real(k_real), pointer, dimension(:,:,:) :: effective_porosity_rate_grid=> null() ! a pointer to the porosity, only used for damage
  real(k_real), pointer, dimension(:,:,:) :: effective_porosity_rate_grid_old=> null() ! a pointer to the porosity, only used for damage

  real(k_real), pointer, dimension(:,:,:,:,:) :: inelastic_strain_rate_grid=> null()
  real(k_real), pointer, dimension(:,:,:,:,:) :: inelastic_strain_grid => null()
  real(k_real), pointer, dimension(:,:,:,:) :: gb_normals_grid => null()
  integer, pointer, dimension(:,:,:) :: gb_id_grid => null()
  real(k_real), pointer, dimension(:,:,:,:) :: void_radius_grid => null()
  real(k_real), pointer, dimension(:,:,:,:) :: void_density_grid => null()
  real(k_real), pointer, dimension(:,:,:,:) :: void_radius_old_grid => null()
  real(k_real), pointer, dimension(:,:,:,:) :: void_density_old_grid => null()

  real(k_real), pointer, dimension(:,:) :: stress_ptr => null() ! a pointer to the stress vector
  real(k_real), pointer :: porosity_ptr => null() ! a pointer to the porosity, only used for damage
  real(k_real), pointer :: porosity_ptr_old => null() ! a pointer to the porosity, only used for damage
  real(k_real), pointer :: porosity_rate_ptr => null() ! a pointer to the porosity, only used for damage
  real(k_real), pointer :: porosity_rate_ptr_old => null() ! a pointer to the porosity, only used for damage
  real(k_real), pointer :: porosity_nucleation_ptr => null()
  real(k_real), pointer :: porosity_nucleation_ptr_old=> null()
  real(k_real), pointer :: porosity_growth_ptr => null()
  real(k_real), pointer :: porosity_growth_ptr_old=> null()
  real(k_real), pointer :: effective_porosity_ptr => null() ! a pointer to the porosity, only used for damage
  real(k_real), pointer :: effective_porosity_ptr_old => null()
  real(k_real), pointer :: effective_porosity_rate_ptr => null() ! a pointer to the porosity, only used for damage
  real(k_real), pointer :: effective_porosity_rate_ptr_old => null() ! a pointer to the porosity, only used for damage

  real(k_real), dimension(:,:), pointer :: inelastic_strain_rate_ptr => null()
  real(k_real), dimension(:,:), pointer :: inelastic_strain_ptr => null()
  real(k_real), pointer, dimension(:) :: gb_normals_ptr => null()
  integer, pointer :: gb_id_ptr => null()
  real(k_real), pointer, dimension(:) :: void_radius_ptr => null()
  real(k_real), pointer, dimension(:) :: void_density_ptr => null()
  real(k_real), pointer, dimension(:) :: void_radius_old_ptr => null()
  real(k_real), pointer, dimension(:) :: void_density_old_ptr => null()
  real(k_real), dimension(13) :: void_radius_bin_limits = (/1e-12_k_real, 1e-11_k_real, 1e-10_k_real, 1e-9_k_real, 1e-8_k_real,&
                                                1e-7_k_real, 1e-6_k_real, 1e-5_k_real, 1e-4_k_real, 1e-3_k_real, &
                                                1e-2_k_real, 1e-1_k_real, 1._k_real/)
  integer :: n_bins = 14
  real(k_real) :: total_porosity
  ! model parameters
  real(k_real) :: min_void_size, & ! the minimum allowed void size after nucleation
                  max_porosity_nucleation, & ! when using chu-needlamn 
                  critical_nucleation_strain, &
                  nucleation_std_dev, &
                  h, & ! the void shape function
                  sin_theta, & !
                  coalescence_critical_porosity, &
                  saturation_porosity, &
                  fustar, & ! coalescence paramter
                  surface_tension, &
                  initial_void_radius, &
                  initial_void_volume, &
                  initial_porosity

  real(k_real), pointer :: omega => null(), &
                           delta_gb => null()
  !
  real(k_real), pointer, dimension(:) :: E_migration => null(), &
                                         E_formation => null(), &
                                         vacancy_diffusivity_coefficient => null(), &
                                         porosity_abs_rel_tol => null()
  ! other variables
  real(k_real) :: stress6(6) ! the stress vector in B basis

  ! use coalescence
  logical :: use_coalescence = .TRUE., &
             use_diffusive_growth =.FALSE., &
             use_bulk_porosity_evolution = .FALSE., &
             use_porosity_nucleation = .FALSE., &
             use_preseeded_porosity = .FALSE., &
             use_plasticity_constrained_growth = .FALSE.
contains

  procedure :: initParameters !-> base function initializing grid_data pointers and allocating the required space
  procedure :: initGridPointers => initGridPointersPorosityBase !-> base function initializing grid_data pointers and allocating the required space

  ! setPointData must always be carryed over trough classes. This is used to set pointers for doing local calculations
  procedure :: setPointData => setPointersPorosityBase
  procedure :: addFieldVariablesToGrid => addFieldVariablesToGridPorosityBase
  procedure :: updatePorosity => updatePorosityBase
  procedure :: computeDiffusivityPorosityBase
  procedure :: updateStateVariablesAtMaterialPointStaggered => updateStateVarsAtMaterialPointStaggeredPorosityBase
  procedure :: updateStateVariablesAtMaterialPointOuterLoop => updateStateVarsAtMaterialPointOuterLoopPorosityBase
  procedure :: readMaterialParametersFromFile => readMaterialParametersFromFilePorosityBase
  procedure :: initStateVariablesAtMaterialPoint => initStateVariablesAtMaterialPointPorosityBase
  procedure :: VoidGrowthGursonPlusDiffusion ! the void growth law
  procedure :: VoidNucleationChuNeedlemann ! the void nucelation rate law
  procedure :: VoidNucleationChuNeedlemannNSites ! the void nucelation rate law
  procedure :: VoidNucleationSize ! the nucleation void size
  procedure :: FindVoidRadiusBin ! the nucleation void size
  procedure :: updateStateVariablesStaggered => updateStateVarsStaggeredPorosityBase
  procedure :: updateStateVariablesOuterLoop => updateStateVarsOuterLoopPorosityBase
  procedure :: acceptRejectSolution => acceptRejectSolutionPorosityBase
  procedure :: computeVoidVolume
  procedure :: computeVoidVolumeRate
  procedure :: writeAverageQuantitiesCSV => writeAverageQuantitiesPorosityBase

end type porosity_base

interface
  module subroutine readMaterialParametersFromFilePorosityBase(this, matf_reader)
    use read_from_file_utils, only : file_reader
    implicit none
    class(porosity_base), intent(inout) :: this
    class(file_reader), intent(inout) :: matf_reader
  end subroutine

  module subroutine addFieldVariablesToGridPorosityBase(this)
    implicit none
    class(porosity_base), intent(inout) :: this
  end subroutine

  module subroutine initParameters(this, phase_id, common_material_parameter_ptr)
    implicit none
    class(porosity_base), intent(inout) :: this
    integer, intent(in) :: phase_id
    class(common_material_parameter), pointer, intent(in) :: common_material_parameter_ptr
  end subroutine

  module subroutine initGridPointersPorosityBase(this)
    implicit none
    class(porosity_base), intent(inout) :: this
  end subroutine

  module subroutine setPointersPorosityBase(this, ix, iy, iz)
    implicit none
    class(porosity_base), intent(inout) :: this
    integer, intent(in) ::  ix, iy, iz
  end subroutine

  module subroutine updatePorosityBase(this)
    implicit none
    class(porosity_base), intent(inout) :: this
  end subroutine

  module subroutine computeDiffusivityPorosityBase(this, D)
    use math_constants, only : kBSI
    implicit none
    class(porosity_base), intent(inout) :: this
    real(k_real), intent(out) :: D
  end subroutine

  module subroutine updateStateVarsAtMaterialPointStaggeredPorosityBase(this)
    implicit none
    class(porosity_base), intent(inout) :: this
  end subroutine

  module subroutine updateStateVarsAtMaterialPointOuterLoopPorosityBase(this)
    implicit none
    class(porosity_base), intent(inout) :: this
  end subroutine

  module subroutine VoidNucleationChuNeedlemann(this, Tn, e_critical, e_std_dev, eq_plastic_strain, eq_plastic_strain_rate, df)
    implicit none
    class(porosity_base), intent(inout) :: this
    real(k_real), intent(in) :: Tn, e_critical, e_std_dev, eq_plastic_strain, eq_plastic_strain_rate
    real(k_real), intent(out) :: df
  end subroutine

  module subroutine VoidNucleationChuNeedlemannNSites(this, Tn, e_critical, e_std_dev, eq_plastic_strain, eq_plastic_strain_rate, dn, df)
    implicit none
    class(porosity_base), intent(inout) :: this
    real(k_real), intent(in) :: Tn, e_critical, e_std_dev, eq_plastic_strain, eq_plastic_strain_rate
    real(k_real), intent(out) :: dn, df
  end subroutine

  module subroutine VoidNucleationSize(this, Tn,  a0)
    implicit none
    class(porosity_base), intent(inout) :: this
    real(k_real), intent(in) :: Tn
    real(k_real), intent(out) :: a0
  end subroutine

  module subroutine computeVoidVolume(this, r, V)
    implicit none
    class(porosity_base), intent(in) :: this
    real(k_real), intent(in) :: r
    real(k_real), intent(out) :: V
  end subroutine
  
  module subroutine computeVoidVolumeRate(this, r, rdot, Vdot)
    implicit none
    class(porosity_base), intent(in) :: this
    real(k_real), intent(in) :: r, rdot
    real(k_real), intent(out) :: Vdot
  end subroutine

  module subroutine VoidGrowthGursonPlusDiffusion(this, a_old, n_density, D, s_VM, Tn, edotp_VM, edotp_H, df_gurson, df_growth)
    use test_utils_mod, only : test_is_finite
    use print_utils_mod, only : printToScreen
    use math_constants, only : PI
    use, intrinsic :: IEEE_ARITHMETIC
    implicit none
    class(porosity_base), intent(inout) :: this
    real(k_real), intent(in) :: a_old, & !-> old void size
                                n_density, & !-> void number density
                                D, &     ! GB diffusivity
                                s_VM,  & ! von mises stress
                                Tn,  &   ! normal traction
                                edotp_VM, &  ! equivalent strain rate
                                edotp_H, & ! hidrostatic plastic strain rate
                                df_gurson
    real(k_real), intent(out) :: df_growth

  end subroutine

  module subroutine FindVoidRadiusBin(this, a, bin_idx)
    implicit none
    class(porosity_base), intent(in) :: this
    real(k_real), intent(in) :: a
    integer, intent(out) :: bin_idx
  end subroutine

  module subroutine initStateVariablesAtMaterialPointPorosityBase(this)
    implicit none
    class(porosity_base), intent(inout) :: this
  end subroutine

  module subroutine updateStateVarsStaggeredPorosityBase(this)
    implicit none
    class(porosity_base), intent(inout) :: this
  end subroutine

  module subroutine updateStateVarsOuterLoopPorosityBase(this)
    implicit none
    class(porosity_base), intent(inout) :: this
  end subroutine

  module subroutine acceptRejectSolutionPorosityBase(this, dt_max, accept_solution_flag)
    implicit none
    class(porosity_base), intent(inout) :: this
    real(k_real), intent(out) :: dt_max
    logical, intent(out) :: accept_solution_flag
  end subroutine

  module subroutine writeAverageQuantitiesPorosityBase(this, csv_writer_obj, write_headers)
    use csv_writer_mod, only : csv_writer
    implicit none
    class(porosity_base), intent(inout) :: this
    class(csv_writer), intent(inout) :: csv_writer_obj
    logical, intent(in) :: write_headers
  end subroutine

end interface

contains

subroutine readMaterialParametersFromFilePorosity(matf_reader, phase_id, &
                                                  all_mighty_grid_in, sim_all_macro_data, the_bc_object, &
                                                  common_material_parameter_ptr, porosity_base_ptr)
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
class(porosity_base), pointer, intent(inout) :: porosity_base_ptr
class(common_material_parameter), pointer, intent(inout) :: common_material_parameter_ptr
type(string_type) :: material_model_name
integer, intent(in) :: phase_id
class(all_mighty_grid_type), intent(in), target :: all_mighty_grid_in
class(sim_all_macro_data_obj), intent(in), target :: sim_all_macro_data
class(boundary_condition_array_type), intent(in), target :: the_bc_object
class(porosity_base), pointer :: porosity_base_temp => null()

      call matf_reader%readParameter("--Porosity-model", material_model_name)
      select case(material_model_name%getString())
      case ("BTKTLC")
        allocate(porosity_base_temp)
        porosity_base_ptr => porosity_base_temp
        call porosity_base_temp%linkBaseObjects(phase_id, all_mighty_grid_in, sim_all_macro_data, the_bc_object)
        call porosity_base_temp%initParameters(phase_id, common_material_parameter_ptr)
  
      case default
        if(write_detailed_log_to_screen) write(*,*) "unrecognized material model for porosity ", material_model_name%getString()
        if(write_detailed_log_to_screen) write(*,*) "available porosity models are: "
        if(write_detailed_log_to_screen) write(*,*) "BTKTLC "
        error stop "abort"
      end select
  
      call porosity_base_ptr%readMaterialParametersFromFile(matf_reader)
      select case(material_model_name%getString())
      case ("BTKTLC")
  
      case default
        if (write_detailed_log_to_screen) write(*,*) "somthing is worng when trying nullifying the pointer for ", material_model_name%getString()
        error stop "abort"
      end select

end subroutine


end module porosity_base_mod
