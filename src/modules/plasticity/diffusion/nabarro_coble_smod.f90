submodule(isotropic_diffusion_mod) nabarro_coble_smod
implicit none

#include "macro_debug.fpp"

contains

module subroutine readMaterialParametersFromFileNabarroCobleDiffusion(this, matf_reader)
  use read_from_file_utils, only : file_reader
  use string_module , only : string_array
  use units_conversion_mod, only : eV2Joule
  implicit none
  class(nabarro_herring_plus_coble_diffusion), intent(inout) :: this
  class(file_reader), intent(inout) :: matf_reader

  call matf_reader%readVectorParameter("diffusion-geometric-factor[bulk,GB-[unitless]]", 2, this%geometric_factor)
  call matf_reader%readParameter("use-RuleOfMixture-Diffusion[TRUE/FALSE]", this%use_RuleOfMixture_Diffusion)
end subroutine

module subroutine initParametersNabarroCobleDiffusion(this, phase_id, common_material_parameter_ptr, use_damage)
  implicit none
  class(nabarro_herring_plus_coble_diffusion), intent(inout) :: this
  integer, intent(in) :: phase_id
  class(common_material_parameter), pointer, intent(in) :: common_material_parameter_ptr
  logical, intent(in) :: use_damage

  call this%initParametersDiffusionBase(phase_id, common_material_parameter_ptr, use_damage)

end subroutine

module subroutine  addFieldVariablesToGridNabarroCobleDiffusion(this)
  use grid_data_var_type_mod
  implicit none
  class(nabarro_herring_plus_coble_diffusion), intent(inout) :: this

  call this%diffusion_base%addFieldVariablesToGrid()
  associate (all_grid_data_vars => this%grid_data)

  end associate
end subroutine

module subroutine initGridPointersNabarroCobleDiffusion(this)
  implicit none
  class(nabarro_herring_plus_coble_diffusion), intent(inout) :: this

  call this%diffusion_base%initGridPointers()

  call this%common_grid_data%getScalarIntegerDataPointerByName("gb_id", this%gb_type_id_grid)

end subroutine

module subroutine initStateVariablesAtMaterialPointNabarroCobleDiffusion(this)
  use test_utils_mod
  implicit none
  class(nabarro_herring_plus_coble_diffusion), intent(inout) :: this
  __DECL_CLASS_UNUSED_THIS__
  ! this model does not update any state variable
  __SUPPRESS_CLASS_UNUSED_THIS__
end subroutine

module subroutine updateStateVarsAtMaterialPointInnerLoopNabarroCobleDiffusion(this)
  implicit none
  class(nabarro_herring_plus_coble_diffusion), intent(inout) :: this
  __DECL_CLASS_UNUSED_THIS__
  ! this model does not update any state variable
  __SUPPRESS_CLASS_UNUSED_THIS__
end subroutine

module subroutine updateStateVarsAtMaterialPointStaggeredNabarroCobleDiffusion(this)
  implicit none
  class(nabarro_herring_plus_coble_diffusion), intent(inout) :: this
  __DECL_CLASS_UNUSED_THIS__
  ! this model does not update any state variable
  __SUPPRESS_CLASS_UNUSED_THIS__
end subroutine

module subroutine setPointDataNabarroCobleDiffusion(this, ix, iy, iz)
  implicit none
  class(nabarro_herring_plus_coble_diffusion), intent(inout) :: this
  integer, intent(in) ::  ix, iy, iz

  call this%diffusion_base%setPointData(ix, iy, iz)
  this%gb_type_id_ptr => this%gb_type_id_grid(ix, iy, iz)

end subroutine

module subroutine computeEpsilonDotAndDepsilonDotDStressNabarroCobleDiffusion(this, stress6, epsilon_dot, depsilon_dot_dstress)
  use math_constants, only : kBSI
  use units_conversion_mod, only : MPA2Pa
  use print_utils_mod
  use global, only : voxel_size
  implicit none
  class(nabarro_herring_plus_coble_diffusion), intent(inout) :: this
  real(k_real), target, intent(in) ::stress6(6)
  real(k_real), target, intent(out) :: epsilon_dot(6), depsilon_dot_dstress(6,6)
  integer :: i , j, gb_id
  real(k_real) :: D_Bulk, C_Bulk, D_GB, C_GB, diff_prefactor_Bulk, diff_prefactor_GB
  real(k_real) :: GB_fraction_for_ROM, Bulk_fraction_for_ROM

  gb_id = min(this%gb_type_id_ptr, 2)
  call this%computeDiffusivity(1, D_Bulk)
  call this%computeVacancyThermalEquilibriumConcentration(1, C_Bulk)
  call this%computeDiffusivity(2, D_GB)
  call this%computeVacancyThermalEquilibriumConcentration(2, C_GB)

  associate( omega => this%common_material_parameter_ptr%atomic_volume, &
             grain_diameter => this%common_material_parameter_ptr%grain_diameter, &
             gb_thickness => this%common_material_parameter_ptr%gb_thickness, &
             T => this%temperature)

  diff_prefactor_Bulk = this%geometric_factor(1) *omega/grain_diameter**2
  diff_prefactor_GB = this%geometric_factor(2) *omega/grain_diameter**3*gb_thickness

  select case(gb_id)
     case (1) !bulk
        GB_fraction_for_ROM = 0._k_real
        Bulk_fraction_for_ROM = 1._k_real - GB_fraction_for_ROM
     case (2) !grain boundary
        GB_fraction_for_ROM = 1._k_real
        if(this%use_RuleOfMixture_Diffusion) GB_fraction_for_ROM = gb_thickness/voxel_size(1)
        ! TODO: The above equation is valid only for cubic voxels. It needs to be revised for general parallelepiped.
        Bulk_fraction_for_ROM = 1._k_real - GB_fraction_for_ROM
     case default
       write(*,*) "Unrecognized grain boundary id: ", gb_id
       error stop "abort"
     end select

  do j=1,6
    epsilon_dot(j) = (stress6(j)*MPa2Pa/(kBSI*T))*(Bulk_fraction_for_ROM*diff_prefactor_Bulk*C_Bulk*D_Bulk + GB_fraction_for_ROM*diff_prefactor_GB*C_GB*D_GB) 
    do i=1,6
      if (i==j) then
        depsilon_dot_dstress(i,j) = (MPa2Pa/(kBSI*T))*(Bulk_fraction_for_ROM*diff_prefactor_Bulk*C_Bulk*D_Bulk + GB_fraction_for_ROM*diff_prefactor_GB*C_GB*D_GB)
      else
        depsilon_dot_dstress(i,j) = 0._k_real
      endif
    enddo
  enddo
  end associate

end subroutine

module subroutine computeDiffusivity(this, gb_id, D)
  use math_constants, only : kBSI
  use units_conversion_mod, only : MPA2Pa
  implicit none
  class(nabarro_herring_plus_coble_diffusion), intent(inout) :: this
  integer, intent(in) :: gb_id
  real(k_real), intent(out) :: D

  associate(Diff => this%common_material_parameter_ptr%vacancy_diffusivity_coefficient, &
            Em => this%common_material_parameter_ptr%vacancy_migration_energy(gb_id), &
            T => this%temperature)
  D = Diff(gb_id) * exp(-Em/(kBSI*T))

  end associate
end subroutine

module subroutine computeVacancyThermalEquilibriumConcentration(this, gb_id, C)
  use math_constants, only : kBSI
  use units_conversion_mod, only : MPA2Pa
  implicit none
  class(nabarro_herring_plus_coble_diffusion), intent(inout) :: this
  integer, intent(in) :: gb_id
  real(k_real), intent(out) :: C

  associate(Ef => this%common_material_parameter_ptr%vacancy_formation_energy(gb_id), &
            T => this%temperature )
  C = exp(-Ef/(kBSI*T))
  end associate
end subroutine

module subroutine acceptRejectSolutionNabarroCobleDiffusion(this, dt_max, accept_solution_flag)
  use log_file_mod, only : writeToScreen
  implicit none
  class(nabarro_herring_plus_coble_diffusion), intent(inout) :: this
  real(k_real), intent(out) :: dt_max
  logical, intent(out) :: accept_solution_flag

  dt_max = this%dt_max
  accept_solution_flag = .true.

  call writeToScreen("acceptRejectSolutionNabarroCobleDiffusion, nothing to check")

end subroutine

end submodule
