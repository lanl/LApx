submodule(climb_mod) wkkct_climb_smod
use test_utils_mod
implicit none

#include "macro_debug.fpp"

contains

module subroutine initParametersWKKCTClimbBase(this, phase_id, common_material_parameter_ptr,&
                                              use_damage, n_gauss, n_std_dev_gauss_integration, &
                                              elasticity_obj, crystal_paremeters_ptr)
  use stiffness_base_mod, only : stiffness_base         
  use cp_base_mod, only : crystal_paremeters_type                              
  implicit none
  class(wkkct_climb_base), intent(inout) :: this
  integer, intent(in) :: phase_id
  class(common_material_parameter), pointer, intent(in) :: common_material_parameter_ptr
  logical, intent(in) :: use_damage
  integer, intent(in):: n_gauss
  real(k_real), intent(in) :: n_std_dev_gauss_integration
  class(stiffness_base), pointer, intent(in) :: elasticity_obj
  class(crystal_paremeters_type), pointer, intent(in) :: crystal_paremeters_ptr

  call this%initParametersClimbBase(phase_id, common_material_parameter_ptr,&
                                    use_damage, n_gauss, n_std_dev_gauss_integration, &
                                    elasticity_obj, crystal_paremeters_ptr)

  allocate(this%burgN(this%n_slip_modes))
  this%burgN = crystal_paremeters_ptr%burgVectorL

  this%elasticity_obj => elasticity_obj
  this%shear_mod_phase_ptr => this%elasticity_obj%G_avg
end subroutine

module subroutine  addFieldVariablesToGridWKKCTClimbBase(this)
  use grid_data_var_type_mod
  implicit none
  class(wkkct_climb_base), intent(inout) :: this

  call this%cp_climb_base%addFieldVariablesToGrid()

  associate (all_grid_data_vars => this%grid_data)

  ! while this model uses the grid variable "rho_m", it is not its responsability
  ! to add it to the grid. in fact this model should not work if a glide material
  ! is not present
  call all_grid_data_vars%addVar("climb_velocity", ss_generic_vector, &
                                 additional_var_dimensions=(/this%n_gauss/), stateful_level=2)
  call all_grid_data_vars%addVar("climb_velocity_CellWall", ss_generic_vector, &
                                 additional_var_dimensions=(/this%n_gauss/), stateful_level=2)
  end associate

end subroutine

module subroutine initGridPointersWKKCTClimbBase(this)
  implicit none
  class(wkkct_climb_base), intent(inout) :: this

  call this%cp_climb_base%initGridPointers()
  call this%grid_data%getSSScalarDataPointerByName("rho_m", this%rho_m_grid)
  call this%grid_data%getSSScalarDataPointerByName("rho_cw", this%rho_cw_grid)
  call this%grid_data%getSSGenericVectorDataPointerByName("climb_velocity", this%climb_velocity_grid)
  call this%grid_data%getSSGenericVectorDataPointerByName("climb_velocity_CellWall", this%climb_velocity_CellWall_grid)
  call this%grid_data%getSSScalarDataPointerByName("sub_grain_size", this%sub_grain_size_grid)

end subroutine

module subroutine setPointDataWKKCTClimbBase(this, ix,iy,iz)
  implicit none
  class(wkkct_climb_base), intent(inout) :: this
  integer, intent(in) :: ix,iy,iz
  associate(nss => this%n_ss, &
            ngauss => this%n_gauss)
  ! betadot is already set in the cp_climb_base class
  call this%cp_climb_base%setPointData(ix,iy,iz)
  this%rho_m_ptr => this%rho_m_grid(1:nss,ix,iy,iz)
  this%rho_cw_ptr => this%rho_cw_grid(1:nss,ix,iy,iz)
  this%climb_velocity_ptr => this%climb_velocity_grid(1:ngauss,1:nss,ix,iy,iz)
  this%climb_velocity_CellWall_ptr => this%climb_velocity_CellWall_grid(1:ngauss,1:nss,ix,iy,iz)
  this%sub_grain_size_ptr => this%sub_grain_size_grid(1:nss,ix,iy,iz)
  end associate

  this%shear_mod = this%shear_mod_phase_ptr

end subroutine

module subroutine initStateVariablesAtMaterialPointWKKCTClimbBase(this)
  implicit none
  class(wkkct_climb_base), intent(inout) :: this
  __DECL_CLASS_UNUSED_THIS__
  ! this model does not requires initializing any state variable
  __SUPPRESS_CLASS_UNUSED_THIS__
end subroutine

module subroutine updateStateVariablesAtMaterialPointInnerLoopWKKCTClimbBase(this)
  implicit none
  class(wkkct_climb_base), intent(inout) :: this
  __DECL_CLASS_UNUSED_THIS__
  ! this model does not update any state variable
  __SUPPRESS_CLASS_UNUSED_THIS__
end subroutine

module subroutine updateStateVariablesAtMaterialPointStaggeredWKKCTClimbBase(this)
  implicit none
  class(wkkct_climb_base), intent(inout) :: this
  __DECL_CLASS_UNUSED_THIS__
  ! this model does not update any state variable
  __SUPPRESS_CLASS_UNUSED_THIS__
end subroutine

module subroutine computeClimbVelocityWKKCTClimbBase(this, rss, vclimb, dvclimb_dtauclimb)
  implicit none
  class(wkkct_climb_base), intent(in) :: this
  real(k_real), intent(in) :: rss
  real(k_real), intent(out) :: vclimb, dvclimb_dtauclimb
  real(k_real) ::  ci_bulk, dci_bulk_dtauclimb, &
                   cv_bulk, dcv_bulk_dtauclimb, &
                   cv_core, dcv_core_dtauclimb, &
                   zi, dzi_dtauclimb, &
                   zv, dzv_dtauclimb, &
                   Di, Dv

  associate (b => this%burgN(this%ss_idx))

  call this%computeVacancyDiffusivity(Dv)
  call this%computeInterstialDiffusivity(Di)
  call this%computeZinterstitial(rss, zi, dzi_dtauclimb)
  call this%computeZVacancy(rss, zv, dzv_dtauclimb)
  call this%computeConcentrationInterstiatialBulk(rss, ci_bulk, dci_bulk_dtauclimb)
  call this%computeConcentrationVacancyBulk(rss, cv_bulk, dcv_bulk_dtauclimb)
  call this%computeConcentrationVacancyCore(rss, cv_core, dcv_core_dtauclimb)

  vclimb = 1._k_real/b*(zi*Di*ci_bulk - &
                    zv*Dv*(cv_bulk - cv_core) )

  dvclimb_dtauclimb = 1._k_real/b*(dzi_dtauclimb*Di*ci_bulk + zi*Di*dci_bulk_dtauclimb - &
                    (dzv_dtauclimb*Dv*(cv_bulk - cv_core) + &
                     zv*Dv*(dcv_bulk_dtauclimb - dcv_core_dtauclimb)) )

  !   if (this%ix*this%iy*this%iz.eq.1.and.i_am_mpi_master.and.this%ss_idx.eq.1.and.this%igauss.eq.1) then
  !   write(*,*)  "ss, igauss, Dv, Di, rss, zi, zv, ci_bulk, cv_bulk, cv_core, vclimb"
  !   write(*,*) this%ss_idx, this%igauss, Dv, Di, rss, zi, zv, ci_bulk, cv_bulk, cv_core, vclimb
  ! endif

  end associate

end subroutine

module subroutine computeClimbVelocityCellWallWKKCTClimbBase(this, rss, vclimb_CW, dvclimb_dtauclimb_CW)
  implicit none
  class(wkkct_climb_base), intent(in) :: this
  real(k_real), intent(in) :: rss
  real(k_real), intent(out) :: vclimb_CW, dvclimb_dtauclimb_CW
  real(k_real) ::  ci_bulk, dci_bulk_dtauclimb, &
                   cv_CW, dcv_CW_dtauclimb, &
                   cv_core, dcv_bulk_dtauclimb, &
                   zi, dzi_dtauclimb, &
                   zv, dzv_dtauclimb, &
                   Di, Dv

  associate (b => this%burgN(this%ss_idx), &
             rho_cw => this%rho_cw_ptr(this%ss_idx))

  call this%computeVacancyDiffusivity(Dv)
  call this%computeInterstialDiffusivity(Di)
  call this%computeZinterstitial(rss, zi, dzi_dtauclimb)
  call this%computeZVacancy(rss, zv, dzv_dtauclimb)
  call this%computeConcentrationInterstiatialBulk(rss, ci_bulk, dci_bulk_dtauclimb)
  call this%computeConcentrationVacancyCellWall(rho_cw, cv_CW, dcv_CW_dtauclimb)
  call this%computeConcentrationVacancyCore(rss, cv_core, dcv_bulk_dtauclimb)

  vclimb_CW = 1._k_real/b*(zi*Di*ci_bulk - &
                    zv*Dv*(cv_CW - cv_core) )

  dvclimb_dtauclimb_CW = 1._k_real/b*(dzi_dtauclimb*Di*ci_bulk + zi*Di*dci_bulk_dtauclimb - &
                          (dzv_dtauclimb*Dv*(cv_CW - cv_core) + &
                          zv*Dv*(dcv_bulk_dtauclimb - dcv_CW_dtauclimb)) )

  end associate

end subroutine

module subroutine computeZinterstitialWKKCTClimbBase(this, tau_climb, zi, dzi_dtauclimb)
  implicit none
  class(wkkct_climb_base), intent(in) :: this
  real(k_real), intent(in) :: tau_climb
  real(k_real), intent(out) :: zi, dzi_dtauclimb
  __DECL_CLASS_UNUSED_THIS__
  __DECL_UNUSED_REAL__

  __SUPPRESS_CLASS_UNUSED_THIS__
  __SUPPRESS_UNUSED_REAL__(tau_climb)
  __SUPPRESS_UNUSED_REAL_OUT__(zi)
  __SUPPRESS_UNUSED_REAL_OUT__(dzi_dtauclimb)

  error stop "if you end up here it meas you didn't override computeConcentrationInterstiatialBulkWKKCTBase"
end subroutine

module subroutine computeZVacancyWKKCTClimbBase(this, tau_climb, zv, dzv_dtauclimb)
  implicit none
  class(wkkct_climb_base), intent(in) :: this
  real(k_real), intent(in) :: tau_climb
  real(k_real), intent(out) :: zv, dzv_dtauclimb
  __DECL_CLASS_UNUSED_THIS__
  __DECL_UNUSED_REAL__

  __SUPPRESS_CLASS_UNUSED_THIS__
  __SUPPRESS_UNUSED_REAL__(tau_climb)
  __SUPPRESS_UNUSED_REAL_OUT__(zv)
  __SUPPRESS_UNUSED_REAL_OUT__(dzv_dtauclimb)

  error stop "if you end up here it meas you didn't override computeZVacancyWKKCTBase"
end subroutine

module subroutine computeConcentrationInterstiatialBulkWKKCTBase(this, tau_climb, ci_bulk, dci_bulk_dtauclimb)
  implicit none
  class(wkkct_climb_base), intent(in) :: this
  real(k_real), intent(in) :: tau_climb
  real(k_real), intent(out) :: ci_bulk, dci_bulk_dtauclimb
  __DECL_CLASS_UNUSED_THIS__
  __DECL_UNUSED_REAL__

  __SUPPRESS_CLASS_UNUSED_THIS__
  __SUPPRESS_UNUSED_REAL__(tau_climb)
  __SUPPRESS_UNUSED_REAL_OUT__(ci_bulk)
  __SUPPRESS_UNUSED_REAL_OUT__(dci_bulk_dtauclimb)

  error stop "if you end up here it meas you didn't override computeConcentrationInterstiatialBulkWKKCTBase"
end subroutine

module subroutine computeConcentrationVacancyBulkWKKCTBase(this, tau_climb, cv_bulk, dcv_bulk_dtauclimb)
  implicit none
  class(wkkct_climb_base), intent(in) :: this
  real(k_real), intent(in) :: tau_climb
  real(k_real), intent(out) :: cv_bulk, dcv_bulk_dtauclimb
  __DECL_CLASS_UNUSED_THIS__
  __DECL_UNUSED_REAL__

  __SUPPRESS_CLASS_UNUSED_THIS__
  __SUPPRESS_UNUSED_REAL__(tau_climb)
  __SUPPRESS_UNUSED_REAL_OUT__(cv_bulk)
  __SUPPRESS_UNUSED_REAL_OUT__(dcv_bulk_dtauclimb)

  error stop "if you end up here it meas you didn't override computeConcentrationVacancyBulkWKKCTBase"
end subroutine

module subroutine computeConcentrationVacancyCellWallWKKCTBase(this, rho_cw, cv_CW, dcv_CW_dtauclimb)
  implicit none
  class(wkkct_climb_base), intent(in) :: this
  real(k_real), intent(in) :: rho_cw
  real(k_real), intent(out) :: cv_CW, dcv_CW_dtauclimb
  __DECL_CLASS_UNUSED_THIS__
  __DECL_UNUSED_REAL__

  __SUPPRESS_CLASS_UNUSED_THIS__
  __SUPPRESS_UNUSED_REAL__(rho_cw)
  __SUPPRESS_UNUSED_REAL_OUT__(cv_CW)
  __SUPPRESS_UNUSED_REAL_OUT__(dcv_CW_dtauclimb)

  error stop "if you end up here it meas you didn't override computeConcentrationVacancyCellWallWKKCTBase"
end subroutine

module subroutine computeConcentrationVacancyCoreWKKCTBase(this, tau_climb, cv_core, dcv_core_dtauclimb)
  implicit none
  class(wkkct_climb_base), intent(in) :: this
  real(k_real), intent(in) :: tau_climb
  real(k_real), intent(out) :: cv_core, dcv_core_dtauclimb
  __DECL_CLASS_UNUSED_THIS__
  __DECL_UNUSED_REAL__

  __SUPPRESS_CLASS_UNUSED_THIS__
  __SUPPRESS_UNUSED_REAL__(tau_climb)
  __SUPPRESS_UNUSED_REAL_OUT__(cv_core)
  __SUPPRESS_UNUSED_REAL_OUT__(dcv_core_dtauclimb)

  error stop "if you end up here it meas you didn't override computeConcentrationVacancyCoreWKKCTBase"
end subroutine


module subroutine computeGammaDotAndDGammaDotDRssSSWKKCTClimbBase(this, rss, gdot, dgdot_drss)
  ! use print_utils_mod, only : printToScreen
  implicit none
  class(wkkct_climb_base), intent(inout) :: this
  real(k_real), intent(in) :: rss
  real(k_real), intent(out) :: gdot
  real(k_real), intent(out) :: dgdot_drss
  real(k_real) :: vclimb, dvclimb_dtauclimb,vclimb_CW,dvclimb_dtauclimb_CW
  associate( rho_ms => this%rho_m_ptr(this%ss_idx), &
             burgN => this%burgN(this%ss_idx), &
             climb_reduce_ratio => this%climb_reduce_ratio, &
             ratio_edge_screw=> this%ratio_edge_screw, &
             tau_climb => rss, &
             betadot => gdot, &
             dbetadot_dtauclimb => dgdot_drss)
  !compute the dislocation velocity
  call this%computeClimbVelocity(tau_climb, vclimb, dvclimb_dtauclimb)
  this%climb_velocity_ptr(this%igauss,this%ss_idx) = vclimb
  
  betadot = climb_reduce_ratio*ratio_edge_screw*rho_ms*burgN * vclimb
  dbetadot_dtauclimb = climb_reduce_ratio*ratio_edge_screw*rho_ms*burgN * dvclimb_dtauclimb

  !compute the climb velocity in the cell wall
  call this%computeClimbVelocityCellWall(tau_climb, vclimb_CW, dvclimb_dtauclimb_CW)
  this%climb_velocity_CellWall_ptr(this%igauss,this%ss_idx) = vclimb_CW
  end associate

end subroutine

module subroutine computeInterstialDiffusivity(this, Di)
  use math_constants, only : kBSI
  implicit none
  class(wkkct_climb_base), intent(in) :: this
  real(k_real), intent(out) :: Di

  associate( Di0 => this%Di0, &
             T => this%temperature, &
             Ei=> this%E_migration_interstitial )

  Di = Di0*exp(-Ei/(kBSI*T))
  end associate
end subroutine

  module subroutine resolved_stress_equilibrium_concentration(this, tau_climb, ctau, dctau_dtauclimb)
    use math_constants, only : kBSI
    use units_conversion_mod, only : MPa2Pa
    implicit none
    class(wkkct_climb_base), intent(in) :: this
    real(k_real), intent(in) :: tau_climb
    real(k_real), intent(out) :: ctau
    real(k_real), intent(out) :: dctau_dtauclimb

    associate( omega => this%burgN(this%ss_idx)**3, &
               T => this%temperature)
    ! stress is in megapascal so we need to move it back to Pascal
    ctau = exp(tau_climb*MPa2Pa*omega/(kBSI*T))
    dctau_dtauclimb =  ctau*omega*MPa2Pa/(kBSI*T)
    end associate
  end subroutine

  module subroutine cellwall_local_concentration(this, rho_cw, c_CW_rho, dcCW_dtauclimb)
    use math_constants, only : kBSI, PI
    use units_conversion_mod, only : MPa2Pa
    implicit none
    class(wkkct_climb_base), intent(in) :: this
    real(k_real), intent(in) :: rho_cw
    real(k_real), intent(out) :: c_CW_rho
    real(k_real), intent(out) :: dcCW_dtauclimb
    real(k_real) :: cell_volume, cellwall_volume

    associate( burgN => this%burgN(this%ss_idx), &
               omega => this%burgN(this%ss_idx)**3, &
               T => this%temperature, &
               shear_mod => this%shear_mod, &
               sub_grain_size => this%sub_grain_size_ptr(this%ss_idx), &
               cellwall_thickness => this%common_material_parameter_ptr%cellwall_thickness)

    cell_volume = (4._k_Real/3._k_Real)*PI*(sub_grain_size/2._k_Real)**3
    cellwall_volume = (4._k_Real/3._k_Real)*PI*( (sub_grain_size/2._k_Real)**3-((sub_grain_size-cellwall_thickness)/2._k_Real)**3 )
    c_CW_rho = exp(omega*burgN*shear_mod*sqrt(rho_cw*cell_volume/cellwall_volume)/(kBSI*T))
    dcCW_dtauclimb =  0._k_Real
    end associate
  end subroutine

  module subroutine hydrostatic_stress_equilibrium_concentration(this, csigma_h, dcsigma_h_dtauclimb)
    use math_constants, only : kBSI
    use units_conversion_mod, only : MPa2Pa
    use tensor_math_mod, only: computePressure, computeDPressureDstress
    implicit none
    class(wkkct_climb_base), intent(in) :: this
    real(k_real), intent(out) :: csigma_h, dcsigma_h_dtauclimb
    real(k_real) :: sh, dsh_dstress(6)

    associate( omega => this%burgN(this%ss_idx)**3, &
               T => this%temperature)

    sh = computePressure(this%stress6)
    dsh_dstress = computeDPressureDstress(this%stress6)

    csigma_h = exp(sh*MPa2Pa*omega/(kBSI*T))
    !TODO WE MISS dcsigmah_dsigma_h, we need to discuss how we want to implemment this
    dcsigma_h_dtauclimb = 0._k_real

    end associate

  end subroutine

!********************************************************************************!
! Subroutines implementation for  Wen et al. 2020 EQ. 20
!********************************************************************************!

module subroutine readMaterialParametersFromFileWKKCTClimb(this, matf_reader)
  use read_from_file_utils, only : file_reader
  use string_module , only : string_array
  use units_conversion_mod, only : eV2Joule
  implicit none
  class(wkkct_climb), intent(inout) :: this
  class(file_reader), intent(inout) :: matf_reader
  type(string_array) :: dummy_string_array
  real(k_real), dimension(this%n_slip_modes) :: BurgN_temp
  call matf_reader%readLineAndCheckStringsAreEqual("-shared-parameters", dummy_string_array)
  call matf_reader%readParameter( "climb-reduced-ratio", this%climb_reduce_ratio)
  this%vacancy_diffusivity_coefficient => this%common_material_parameter_ptr%vacancy_diffusivity_coefficient ! point to the bulk value
  this%E_formation_vacancy => this%common_material_parameter_ptr%vacancy_formation_energy ! point to the bulk value
  this%E_migration_vacancy => this%common_material_parameter_ptr%vacancy_migration_energy ! point to the bulk value
  call matf_reader%readParameter( "dislocation-capture-efficiency-vacancy[zv0]", this%zv0)
  call matf_reader%readParameter( "n-exponent-damage[unitless]", this%n_exp_damage)

  ! some sanity checks
  if (this%climb_reduce_ratio<0._k_real) error stop "WKKCTClimb climb_reduce_ratio <0."
  if (this%climb_reduce_ratio>1._k_real) error stop "WKKCTClimb climb_reduce_ratio >1."

  if (any(this%E_formation_vacancy<0._k_real)) error stop "WKKCTClimb G_f <=0."
  if (any(this%vacancy_diffusivity_coefficient<=0._k_real)) error stop "WKKCTClimb vacancy_diffusivity_coefficient <=0."
  if (any(this%E_migration_vacancy<=0._k_real)) error stop "WKKCTClimb E_migration_vacancy <=0."
  if (this%zv0<=0._k_real) error stop "WKKCTClimb zv0 <=0."

  ! since this model does not use interstitials we forcefull yset interstiail parameters to 0
  this%Di0 = 0._k_real
  this%zi0 = 0._k_real
  this%E_migration_interstitial = 0._k_real

  BurgN_temp = this%burgN
  deallocate(this%burgN); nullify(this%burgN)
  call this%convertModeToSlipSystemParameter(BurgN_temp, this%burgN)

end subroutine


module subroutine initParametersWKKCTClimb(this, phase_id, common_material_parameter_ptr,&
                                          use_damage, n_gauss, n_std_dev_gauss_integration, &
                                          elasticity_obj, crystal_paremeters_ptr)
    use stiffness_base_mod, only : stiffness_base         
    use cp_base_mod, only : crystal_paremeters_type                              
    implicit none
    class(wkkct_climb), intent(inout) :: this
    integer, intent(in) :: phase_id
    class(common_material_parameter), pointer, intent(in) :: common_material_parameter_ptr
    logical, intent(in) :: use_damage
    integer, intent(in):: n_gauss
    real(k_real), intent(in) :: n_std_dev_gauss_integration
    class(stiffness_base), pointer, intent(in) :: elasticity_obj
    class(crystal_paremeters_type), pointer, intent(in) :: crystal_paremeters_ptr

  call this%initParametersWKKCTClimbBase(phase_id, common_material_parameter_ptr, &
                                          use_damage, n_gauss, n_std_dev_gauss_integration, &
                                          elasticity_obj, crystal_paremeters_ptr)

end subroutine

module subroutine initGridPointersWKKCTClimb(this)
  implicit none
  class(wkkct_climb), intent(inout) :: this

  call this%wkkct_climb_base%initGridPointers()
end subroutine

module subroutine setPointDataWKKCTClimb(this, ix, iy, iz)
  implicit none
  class(wkkct_climb), intent(inout) :: this
  integer, intent(in) :: ix,iy,iz

  ! gammadot is already set in the cp_base_class
  call this%wkkct_climb_base%setPointData(ix,iy,iz)
end subroutine

module subroutine computeConcentrationInterstiatialBulkWKKCTClimb(this, tau_climb, ci_bulk, dci_bulk_dtauclimb)
  implicit none
  class(wkkct_climb), intent(in) :: this
  real(k_real), intent(in) :: tau_climb
  real(k_real), intent(out) :: ci_bulk, dci_bulk_dtauclimb
  ! this model does not consider interstitial
  __DECL_CLASS_UNUSED_THIS__
  __DECL_UNUSED_REAL__

  __SUPPRESS_CLASS_UNUSED_THIS__
  __SUPPRESS_UNUSED_REAL__(tau_climb)

  ci_bulk = 0._k_Real
  dci_bulk_dtauclimb = 0._k_real
end subroutine

module subroutine computeZinterstitialWKKCTClimb(this, tau_climb, zi, dzi_dtauclimb)
  implicit none
  class(wkkct_climb), intent(in) :: this
  real(k_real), intent(in) :: tau_climb
  real(k_real), intent(out) :: zi, dzi_dtauclimb
  __DECL_UNUSED_REAL__
  __DECL_CLASS_UNUSED_THIS__

  __SUPPRESS_CLASS_UNUSED_THIS__
  __SUPPRESS_UNUSED_REAL__(tau_climb)
  zi = 0._k_real
  dzi_dtauclimb = 0._k_real

end subroutine

module subroutine computeZVacancyWKKCTClimb(this, tau_climb, zv, dzv_dtauclimb)
  implicit none
  class(wkkct_climb), intent(in) :: this
  real(k_real), intent(in) :: tau_climb
  real(k_real), intent(out) :: zv, dzv_dtauclimb
  __DECL_UNUSED_REAL__

  __SUPPRESS_UNUSED_REAL__(tau_climb)
  zv = this%zv0
  dzv_dtauclimb = 0._k_real
end subroutine


module subroutine computeConcentrationVacancyBulkWKKCTClimb(this, tau_climb, cv_bulk, dcv_bulk_dtauclimb)
  implicit none
  class(wkkct_climb), intent(in) :: this
  real(k_real), intent(in) :: tau_climb
  real(k_real), intent(out) :: cv_bulk, dcv_bulk_dtauclimb
  real(k_real) :: c_th, dcth_dtauclimb
  __DECL_UNUSED_REAL__

  __SUPPRESS_UNUSED_REAL__(tau_climb)
  call this%thermal_equilibrium_concentration(c_th, dcth_dtauclimb)

  cv_bulk = c_th
  dcv_bulk_dtauclimb = dcth_dtauclimb

end subroutine

module subroutine computeConcentrationVacancyCellWallWKKCTClimb(this, rho_cw, cv_CW, dcv_CW_dtauclimb)
  implicit none
  class(wkkct_climb), intent(in) :: this
  real(k_real), intent(in) :: rho_cw
  real(k_real), intent(out) :: cv_CW, dcv_CW_dtauclimb
  real(k_real) :: c_th, dcth_dtauclimb, c_CW_rho, dcCW_dtauclimb

  call this%thermal_equilibrium_concentration(c_th, dcth_dtauclimb)
  call this%cellwall_local_concentration(rho_cw, c_CW_rho, dcCW_dtauclimb)
  cv_CW = c_th*c_CW_rho
  dcv_CW_dtauclimb = dcth_dtauclimb + c_th*dcCW_dtauclimb

end subroutine

module subroutine computeConcentrationVacancyCoreWKKCTClimb(this, tau_climb, cv_core, dcv_core_dtauclimb)
  implicit none
  class(wkkct_climb), intent(in) :: this
  real(k_real), intent(in) :: tau_climb
  real(k_real), intent(out) :: cv_core, dcv_core_dtauclimb
  real(k_real) :: c_th, dcth_dtauclimb, c_tau, dctau_dtauclimb

  call this%thermal_equilibrium_concentration(c_th, dcth_dtauclimb)
  call this%resolved_stress_equilibrium_concentration(tau_climb, c_tau, dctau_dtauclimb)

  cv_core = c_th*c_tau
  dcv_core_dtauclimb = dcth_dtauclimb*c_tau + c_th*dctau_dtauclimb

end subroutine

module subroutine acceptRejectSolutionWKKCTClimb(this, dt_max, accept_solution_flag)
  use log_file_mod, only : writeToScreen
  implicit none
  class(wkkct_climb), intent(inout) :: this
  real(k_real), intent(out) :: dt_max
  logical, intent(out) :: accept_solution_flag

  dt_max = this%dt_max
  accept_solution_flag = .true.

  call writeToScreen("acceptRejectSolutionWKKCTClimb, nothing to check")

end subroutine

!***********************************************************************************!
! WKKCT climb with irradiation. This model implements Wen et al. 2020, Eq 19, 23, 24
!***********************************************************************************!
module subroutine readMaterialParametersFromFileWKKCTClimbIrradiation(this, matf_reader)
  use read_from_file_utils, only : file_reader
  use string_module , only : string_array
  use units_conversion_mod, only : eV2Joule
  implicit none
  class(wkkct_climb_irradiation), intent(inout) :: this
  class(file_reader), intent(inout) :: matf_reader
  type(string_array) :: dummy_string_array
  real(k_real), dimension(this%n_slip_modes) :: BurgN_temp

  call matf_reader%readLineAndCheckStringsAreEqual("-shared-parameters", dummy_string_array)
  call matf_reader%readParameter( "climb-reduced-ratio", this%climb_reduce_ratio)
  call matf_reader%readParameter( "defect-pair-creation-rate[K0]", this%pK0)
  call matf_reader%readParameter( "n-exponent-damage[unitless]", this%n_exp_damage)

  ! sanity checks
  if (this%climb_reduce_ratio<0._k_real) error stop "WKKCTClimb climb_reduce_ratio <0."
  if (this%climb_reduce_ratio>1._k_real) error stop "WKKCTClimb climb_reduce_ratio >1."

  this%vacancy_diffusivity_coefficient => this%common_material_parameter_ptr%vacancy_diffusivity_coefficient ! point to the bulk value
  this%E_formation_vacancy => this%common_material_parameter_ptr%vacancy_formation_energy ! point to the bulk value
  this%E_migration_vacancy => this%common_material_parameter_ptr%vacancy_migration_energy ! point to the bulk value
  call matf_reader%readParameter( "dislocation-capture-efficiency-vacancy[zv0]", this%zv0)

  ! some sanity checks
  if (any(this%E_formation_vacancy<0._k_real)) error stop "WKKCTClimb G_f <=0."
  if (any(this%vacancy_diffusivity_coefficient<=0._k_real)) error stop "WKKCTClimb vacancy_diffusivity_coefficient <=0."
  if (any(this%E_migration_vacancy<0._k_real)) error stop "WKKCTClimb E_migration_vacancy <=0."
  if (this%zv0<=0._k_real) error stop "WKKCTClimb zv0 <=0."

  call matf_reader%skipEmptyLine()
  call matf_reader%readLineAndCheckStringsAreEqual("-interstitial-related-parameters", dummy_string_array)
  call matf_reader%readParameter( "interstitial-diffusivity-coefficient[Di0]", this%Di0)
  call matf_reader%readParameter( "interstitial-migration-energy[Ei]", this%E_migration_interstitial)
  this%E_migration_interstitial = this%E_migration_interstitial*eV2Joule
  call matf_reader%readParameter( "stress-free-capture-efficiency-interstitial[zi0]", this%zi0)
  call matf_reader%readParameter( "scaling-factor-stress-induced-capture-efficiency[Zsi]", this%Zsi)

  ! some sanity checks
  if (this%E_migration_interstitial<0._k_real) error stop "WKKCTClimb E_migration_interstitial <=0."
  if (this%Di0<=0._k_real) error stop "WKKCTClimb Di0 <=0."
  if (this%zi0<=0._k_real) error stop "WKKCTClimb zi0 <=0."
  if (this%Zsi<=0._k_real) error stop "WKKCTClimb Zsi <=0."

  BurgN_temp = this%burgN
  deallocate(this%burgN); nullify(this%burgN)
  call this%convertModeToSlipSystemParameter(BurgN_temp, this%burgN)

end subroutine

module subroutine initParametersWKKCTClimbIrradiation(this, phase_id, common_material_parameter_ptr,&
                                                      use_damage, n_gauss, n_std_dev_gauss_integration, &
                                                      elasticity_obj, crystal_paremeters_ptr)
    use stiffness_base_mod, only : stiffness_base         
    use cp_base_mod, only : crystal_paremeters_type                              
    implicit none
    class(wkkct_climb_irradiation), intent(inout) :: this
    integer, intent(in) :: phase_id
    class(common_material_parameter), pointer, intent(in) :: common_material_parameter_ptr
    logical, intent(in) :: use_damage
    integer, intent(in):: n_gauss
    real(k_real), intent(in) :: n_std_dev_gauss_integration
    class(stiffness_base), pointer, intent(in) :: elasticity_obj
    class(crystal_paremeters_type), pointer, intent(in) :: crystal_paremeters_ptr

  call this%initParametersWKKCTClimbBase(phase_id, common_material_parameter_ptr, &
                                          use_damage, n_gauss, n_std_dev_gauss_integration, &
                                          elasticity_obj, crystal_paremeters_ptr)

end subroutine

module subroutine computeZinterstitialWKKCTClimbIrradiation(this, tau_climb, zi, dzi_dtauclimb)
  implicit none
  class(wkkct_climb_irradiation), intent(in) :: this
  real(k_real), intent(in) :: tau_climb
  real(k_real), intent(out) :: zi, dzi_dtauclimb

  zi = this%zi0 + this%Zsi*tau_climb
  dzi_dtauclimb = this%Zsi

end subroutine

module subroutine computeZVacancyWKKCTClimbIrradiation(this, tau_climb, zv, dzv_dtauclimb)
  implicit none
  class(wkkct_climb_irradiation), intent(in) :: this
  real(k_real), intent(in) :: tau_climb
  real(k_real), intent(out) :: zv, dzv_dtauclimb
  __DECL_UNUSED_REAL__

  __SUPPRESS_UNUSED_REAL__(tau_climb)
  zv = this%zv0
  dzv_dtauclimb = 0._k_real
end subroutine

module subroutine  addFieldVariablesToGridWKKCTClimbIrradiation(this)
  use grid_data_var_type_mod
  implicit none
  class(wkkct_climb_irradiation), intent(inout) :: this

  call this%cp_climb_base%addFieldVariablesToGrid()

  associate (all_grid_data_vars => this%grid_data)

  ! while this model uses the grid variable "rho_m", it is not its responsability
  ! to add it to the grid. in fact this model should not work if a glide material
  ! is not present
  call all_grid_data_vars%addVar("con_ig", scalar)
  call all_grid_data_vars%addVar("con_vg", scalar)
  end associate

end subroutine

module subroutine initGridPointersWKKCTClimbIrradiation(this)
  implicit none
  class(wkkct_climb_irradiation), intent(inout) :: this

  call this%wkkct_climb_base%initGridPointers()
  call this%grid_data%getScalarDataPointerByName("con_ig", this%conc_interstitial_grid)
  call this%grid_data%getScalarDataPointerByName("con_vg", this%conc_vacancy_grid)

end subroutine


module subroutine setPointDataWKKCTClimbIrradiation(this, ix,iy,iz)
  implicit none
  class(wkkct_climb_irradiation), intent(inout) :: this
  integer, intent(in) :: ix,iy,iz

  call this%wkkct_climb_base%setPointData(ix,iy,iz)
  this%conc_interstitial_ptr => this%conc_interstitial_grid(ix,iy,iz)
  this%conc_vacancy_ptr => this%conc_vacancy_grid(ix,iy,iz)

end subroutine

module subroutine  computeSaturationConcentrationWKKCTClimbIrradiation(this)
  implicit none
  class(wkkct_climb_irradiation), intent(inout) :: this

  real(k_real) :: c_th, dcth_dtauclimb, &
                  c_tau, dctau_dtauclimb, &
                  x1, x2, x3, rho_ms, zi, zv, dummy, &
                  Di, Dv


  associate(c_i => this%conc_interstitial_ptr, &
            c_v => this%conc_vacancy_ptr, &
            rho_m =>this%rho_m_ptr, &
            clir => this%climb_reduce_ratio, &
            ratio_edge_screw => this%ratio_edge_screw, &
            tau_climb =>this%rss, &
            ss_idx => this%ss_idx)


    call this%computeResolvedStressAndReslovedStressDStress(this%stress6)
    call this%thermal_equilibrium_concentration(c_th, dcth_dtauclimb)
    call this%computeVacancyDiffusivity(Dv)
    call this%computeInterstialDiffusivity(Di)

    x1 = 0._k_real; x2 = 0._k_real; x3 = 0._k_real;
    do ss_idx=1,this%n_ss
      call this%computeZVacancy(tau_climb(ss_idx), zv, dummy)
      call this%computeZinterstitial(tau_climb(ss_idx), zi, dummy)
      call this%resolved_stress_equilibrium_concentration(tau_climb(ss_idx), c_tau, dctau_dtauclimb)
      rho_ms = clir*ratio_edge_screw*rho_m(ss_idx)
      x1 = x1 + zv*rho_ms
      x2 = x2 + zv*rho_ms*c_tau
      x3 = x3 + zi*rho_ms
    enddo
    c_i=this%pk0/(x3*Di)
    c_v=(this%pk0+x2*Dv*c_th)/(x1*Dv)
  end associate
end subroutine

  module subroutine initStateVarsAtMaterialPointWKKCTClimbIrradiation(this)
    implicit none
      class(wkkct_climb_irradiation), intent(inout) :: this
      call this%computeSaturationConcentrationWKKCTClimbIrradiation
  end subroutine

  module subroutine updateStateVarsAtMaterialPointInnerLoopWKKCTClimbIrradiation(this)
    implicit none
    class(wkkct_climb_irradiation), intent(inout) :: this
    __DECL_CLASS_UNUSED_THIS__
    ! we update everything explicitly after convergence, so nothing to do here
    __SUPPRESS_CLASS_UNUSED_THIS__
  end subroutine

  module subroutine updateStateVarsAtMaterialPointStaggeredWKKCTClimbIrradiation(this)
    implicit none
    class(wkkct_climb_irradiation), intent(inout) :: this
    call this%computeSaturationConcentrationWKKCTClimbIrradiation
  end subroutine


  module subroutine computeConcentrationInterstiatialBulkWKKCTClimbIrradiation(this, tau_climb, ci_bulk, dci_bulk_dtauclimb)
    implicit none
    class(wkkct_climb_irradiation), intent(in) :: this
    real(k_real), intent(in) :: tau_climb
    real(k_real), intent(out) :: ci_bulk, dci_bulk_dtauclimb
    __DECL_UNUSED_REAL__

    __SUPPRESS_UNUSED_REAL__(tau_climb)
    ci_bulk = this%conc_interstitial_ptr
    dci_bulk_dtauclimb = 0
  end subroutine

  module subroutine computeConcentrationVacancyBulkWKKCTClimbIrradiation(this, tau_climb, cv_bulk, dcv_bulk_dtauclimb)
    implicit none
    class(wkkct_climb_irradiation), intent(in) :: this
    real(k_real), intent(in) :: tau_climb
    real(k_real), intent(out) :: cv_bulk, dcv_bulk_dtauclimb
    __DECL_UNUSED_REAL__

    __SUPPRESS_UNUSED_REAL__(tau_climb)
    cv_bulk = this%conc_vacancy_ptr
    dcv_bulk_dtauclimb = 0
  end subroutine

  module subroutine computeConcentrationVacancyCoreWKKCTClimbIrradiation(this, tau_climb, cv_core, dcv_core_dtauclimb)
    implicit none
    class(wkkct_climb_irradiation), intent(in) :: this
    real(k_real), intent(in) :: tau_climb
    real(k_real), intent(out) :: cv_core, dcv_core_dtauclimb

    real(k_real) :: c_th, dcth_dtauclimb, &
                    c_tau, dctau_dtauclimb

    call this%thermal_equilibrium_concentration(c_th, dcth_dtauclimb)
    call this%resolved_stress_equilibrium_concentration(tau_climb, c_tau, dctau_dtauclimb)

    cv_core = c_th*c_tau
    dcv_core_dtauclimb = dcth_dtauclimb*c_tau + c_th*dctau_dtauclimb

  end subroutine





! module subroutine initParametersWKKCTClimbChemoMechanical(this, phase_id, n_ss_per_mode, use_damage, n_gauss, &
!                                         burgN, climb_reduce_ratio, ratio_edge_screw, zv0, diff_v, c_v0)
!   implicit none
!   class(wkkct_climb_chemo_mech), intent(inout) :: this
!   integer, intent(in) :: phase_id, n_ss_per_mode(:)
!   real(k_real), target, intent(in), dimension(:) :: burgN
!   real(k_real), target, intent(in) :: climb_reduce_ratio, ratio_edge_screw, zv0, diff_v, c_v0
!
!   logical, intent(in) :: use_damage
!   integer, intent(in):: n_gauss
!
!   if (n_gauss < 2 ) error stop "WKKCT climb rule n_gauss must be > 1"
!   call this%initParametersWKKCTClimbBase(phase_id, n_ss_per_mode, use_damage, n_gauss, &
!                                           burgN, climb_reduce_ratio, ratio_edge_screw )
!
!   this%z_vacancy0 = zv0
!   this%diff_vacancy = diff_v
!   this%conc_vacancy_0 = c_v0
!
! end subroutine
!
! module subroutine initGridPointersWKKCTClimbChemoMechanical(this)
!   implicit none
!   class(wkkct_climb_chemo_mech), intent(inout) :: this
!
!   call this%wkkct_climb_base%initGridPointers()
!   call this%grid_data%getScalarDataPointerByName("concentration", this%concentration_grid)
!
! end subroutine
!
! module subroutine setPointDataWKKCTClimbChemoMechanical(this, ix,iy,iz)
!   implicit none
!   class(wkkct_climb_chemo_mech), intent(inout) :: this
!   integer, intent(in) :: ix,iy,iz
!
!   ! gamamdot is already set in the cp_base_class
!   call this%wkkct_climb_base%setPointData(ix,iy,iz)
!   this%concentration_ptr => this%concentration_grid(ix,iy,iz)
!
! end subroutine
!
! module subroutine computeClimbVelocityWKKCTClimbChemoMechanical(this, rss)
!   use math_constants, only : kBOLTZMANN
!   use tensor_math_mod, only: computePressure
!   implicit none
!   class(wkkct_climb_chemo_mech), intent(inout) :: this
!   real(k_real), intent(in) :: rss
!   real(k_real) :: zi0
!   associate(v_climb => this%v_climb, &
!             burgN => this%burgN(this%ss_idx), &
!             omega=> this%burgN(this%ss_idx)**3, &
!             temperature => this%temperature, &
!             diff_v => this%diff_vacancy, &
!             c_v0 => this%conc_vacancy_0, &
!             z_v0 => this%z_vacancy0, &
!             sh => computePressure(this%stress6), &
!             conc => this%concentration_ptr &
!             )
!
!   v_climb = omega/burgN * (z_v0*diff_v* &
!           2._k_real/omega *(c_v0 * exp(-sh*omega/(kBOLTZMANN*temperature)) * &
!                                    exp(rss*omega/(kBOLTZMANN*temperature)) - conc ) )
!
!   end associate
! end subroutine
!
end submodule
