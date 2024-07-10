#include "macro_debug.fpp"

submodule(glide_mod) hutchinson_smod
implicit none
contains
!
module subroutine initParametersHutchinsonGlide(this, phase_id, common_material_parameter_ptr, &
                                             use_damage, n_gauss, n_std_dev_gauss_integration, &
                                            elasticity_obj, crystal_paremeters_ptr)
  use stiffness_base_mod, only : stiffness_base  
  use cp_base_mod, only : crystal_paremeters_type                                          
  implicit none
  class(hutchinson_glide), intent(inout) :: this
  integer, intent(in) :: phase_id 
  class(common_material_parameter), pointer, intent(in) :: common_material_parameter_ptr
  logical, intent(in) :: use_damage
  integer, intent(in):: n_gauss
  real(k_real), intent(in) :: n_std_dev_gauss_integration
  class(stiffness_base), pointer, intent(in) :: elasticity_obj
  class(crystal_paremeters_type), pointer, intent(in) :: crystal_paremeters_ptr
  call this%initParametersGlideBase(phase_id, common_material_parameter_ptr, &
                                                use_damage, n_gauss, n_std_dev_gauss_integration, &
                                                elasticity_obj, crystal_paremeters_ptr)
end subroutine

module subroutine readMaterialParametersFromFileHutchinsonGlide(this, matf_reader)
  use read_from_file_utils, only : file_reader
  use string_module , only : string_array
  implicit none
  class(hutchinson_glide), intent(inout) :: this
  class(file_reader), intent(inout) :: matf_reader
  real(k_real), dimension(:), pointer :: mode_dep_parameter => null()
  type(string_array) :: dummy_string_array

  call matf_reader%readVectorParameter("gamma-dot-0[unitless]", this%n_slip_modes, mode_dep_parameter)
  call this%convertModeToSlipSystemParameter(mode_dep_parameter, this%gdot0)
  deallocate(mode_dep_parameter); nullify(mode_dep_parameter)
  call matf_reader%readVectorParameter("n-exponent[unitless]", this%n_slip_modes, mode_dep_parameter)
  call this%convertModeToSlipSystemParameter(mode_dep_parameter, this%n_exp)
  deallocate(mode_dep_parameter); nullify(mode_dep_parameter)
  call matf_reader%readParameter("n-exponent-damage[unitless]", this%n_exp_damage)
  call matf_reader%readVectorParameter("initial-critical-resolved-shear-stress[MPa]", this%n_slip_modes, mode_dep_parameter)
  call this%convertModeToSlipSystemParameter(mode_dep_parameter, this%crss0)
  deallocate(mode_dep_parameter); nullify(mode_dep_parameter)
  call matf_reader%readParameter("rss-std-deviation[MPa]", this%rss_std_dev)
  call matf_reader%skipEmptyLine()

  call matf_reader%readLineAndCheckStringsAreEqual( "-hardening-parameters", dummy_string_array)
  call matf_reader%readParameter("use-hardening[TRUE/FALSE]", this%use_hardening)
  call matf_reader%readVectorParameter("tau1_hardening[MPa]", this%n_slip_modes, mode_dep_parameter)
  call this%convertModeToSlipSystemParameter(mode_dep_parameter, this%tau1_hardening)
  deallocate(mode_dep_parameter); nullify(mode_dep_parameter)
  call matf_reader%readVectorParameter("theta0_hardening[MPa]", this%n_slip_modes, mode_dep_parameter)
  call this%convertModeToSlipSystemParameter(mode_dep_parameter, this%theta0_hardening)
  deallocate(mode_dep_parameter); nullify(mode_dep_parameter)
  call matf_reader%readVectorParameter("theta1_hardening[MPa]", this%n_slip_modes, mode_dep_parameter)
  call this%convertModeToSlipSystemParameter(mode_dep_parameter, this%theta1_hardening)
  deallocate(mode_dep_parameter); nullify(mode_dep_parameter)
  call matf_reader%skipEmptyLine()

  call matf_reader%readLineAndCheckStringsAreEqual( "-time-march-tolerances", dummy_string_array)
  call matf_reader%readVectorParameter("crss-absolute-and-relative-tol[MPa,unitless]", 2, this%crss_abs_rel_tol)
end subroutine

module subroutine  addFieldVariablesToGridHutchinsonGlide(this)
  use grid_data_var_type_mod
  implicit none
  class(hutchinson_glide), intent(inout) :: this

  call this%cp_glide_base%addFieldVariablesToGrid()
  associate (all_grid_data_vars => this%grid_data)

  call all_grid_data_vars%addVar("crss", ss_scalar, stateful_level=2)
  call all_grid_data_vars%addVar("crss_rate", ss_scalar, stateful_level=2)
  call all_grid_data_vars%addVar("accumulated_gamma", scalar, stateful_level=2)
  end associate
end subroutine

module subroutine initGridPointersHutchinsonGlide(this)
  implicit none
  class(hutchinson_glide), intent(inout) :: this

  call this%cp_glide_base%initGridPointers()

  call this%grid_data%getSSScalarDataPointerByName("crss", this%crss_grid)
  call this%grid_data%getSSScalarDataPointerByName("crss_rate", this%crss_rate_grid)
  call this%grid_data%getSSScalarDataPointerByNameOld("crss", this%crss_old_grid)
  call this%grid_data%getSSScalarDataPointerByNameOld("crss_rate", this%crss_rate_old_grid)
  call this%grid_data%getScalarDataPointerByName("accumulated_gamma", this%gamma_accumulated_grid)
  call this%grid_data%getScalarDataPointerByNameOld("accumulated_gamma", this%gamma_accumulated_old_grid)
end subroutine

module subroutine setPointDataHutchinsonGlide(this, ix,iy,iz)
  implicit none
  class(hutchinson_glide), intent(inout) :: this
  integer, intent(in) :: ix,iy,iz

  !calling the parent class to do its job
  call this%cp_glide_base%setPointData(ix,iy,iz)
  this%crss_ptr => this%crss_grid(:,ix,iy,iz)
  this%crss_old_ptr => this%crss_old_grid(:,ix,iy,iz)
  this%crss_rate_ptr => this%crss_rate_grid(:,ix,iy,iz)
  this%crss_rate_old_ptr => this%crss_rate_old_grid(:,ix,iy,iz)
  this%accumulated_gamma_ptr => this%gamma_accumulated_grid(ix,iy,iz)
  this%accumulated_gamma_old_ptr => this%gamma_accumulated_old_grid(ix,iy,iz)
end subroutine
!
module subroutine computeGammaDotAndDGammaDotDRssSSHutchinsonGlide(this, rss, gdot, dgdot_drss)
  implicit none
  class(hutchinson_glide), intent(inout) :: this
  real(k_real), intent(in) :: rss
  real(k_real), intent(out) :: gdot
  real(k_real), intent(out) :: dgdot_drss

  associate (gdot0 => this%gdot0, &
             crss => this%crss_ptr, &
             n => this%n_exp, &
             ss_idx => this%ss_idx )


    gdot = gdot0(ss_idx)*(abs(rss/crss(ss_idx)))**n(ss_idx) * &
                            sign(1._k_real,rss)

    ! we are ssuming dcrss/drss = 0
    dgdot_drss = n(ss_idx) * gdot0(ss_idx)*&
                        (abs(rss/crss(ss_idx)))**(n(ss_idx)-1) / crss(ss_idx)

  end associate
end subroutine

module subroutine initStateVariablesAtMaterialPointHutchinsonGlide(this)
  use test_utils_mod
  implicit none
  class(hutchinson_glide), intent(inout) :: this

  this%crss_ptr = this%crss0(1)
  if (this%n_gauss > 1) then
  this%std_dev_ptr = this%rss_std_dev
  endif
end subroutine

module subroutine updateStateVariablesAtMaterialPointInnerLoopHutchinsonGlide(this)
  implicit none
  class(hutchinson_glide), intent(inout) :: this
  __DECL_CLASS_UNUSED_THIS__

  __SUPPRESS_CLASS_UNUSED_THIS__
  ! nothing to do here
end subroutine

module subroutine updateStateVariablesAtMaterialPointStaggeredHutchinsonGlide(this)
  implicit none
  class(hutchinson_glide), intent(inout) :: this
  real(k_real) :: delta_gamma, delta_gamma_H(this%n_ss), exp_gamma(this%n_ss), dtauvoce_ddeltagamma(this%n_ss)
  if (this%use_hardening) then
  !this works per material point wise
  delta_gamma = sum(abs(this%gammadot_ptr))*this%dt
  delta_gamma_H = matmul(this%h_alpha_beta, abs(this%gammadot_ptr))*this%dt
  this%accumulated_gamma_ptr = this%accumulated_gamma_old_ptr + delta_gamma
  
  if (any(this%tau1_hardening.ne.0._k_real)) then
    exp_gamma = exp(-this%accumulated_gamma_ptr*this%theta0_hardening/this%tau1_hardening)
    dtauvoce_ddeltagamma = this%theta1_hardening * (1._k_real - exp_gamma) +  &
    (this%tau1_hardening + this%theta1_hardening * this%accumulated_gamma_ptr) * exp_gamma * this%theta0_hardening/this%tau1_hardening
  else 
    ! this is the linear hardening case which happens when tau1 = 0.
    dtauvoce_ddeltagamma = this%theta0_hardening
  endif

  this%crss_ptr = this%crss_old_ptr + dtauvoce_ddeltagamma * delta_gamma_H
  this%crss_rate_ptr = dtauvoce_ddeltagamma * delta_gamma_H/this%dt
  endif
end subroutine

module subroutine acceptRejectSolutionHutchinsonGlide(this, dt_max, accept_solution_flag)
  use mpi_useful_routines_mod, only : MPIComputeQualitySSScalar, MPIMaxIncrementGridSSScalar
  use mpi_variables_mod, only : i_am_mpi_master
  use time_march_mod, only : computeMaxAllowedTimeStepQuality
  implicit none
  class(hutchinson_glide), intent(inout) :: this
  real(k_real), intent(out) :: dt_max
  logical, intent(out) :: accept_solution_flag
  real(k_real) :: max_crss_increment
  real(k_real) :: max_dt_var
  real(k_real) :: Q
  integer :: npoints, dims(4)
  logical :: accept_var

  if (i_am_mpi_master) write(*,*) "*********************************"
  if (i_am_mpi_master) write(*,*) "acceptRejectSolutionHutchinsonGlide"

  accept_solution_flag = .true.
  dt_max = this%dt_max

  npoints = this%grid_data%getGlobalGridNPoints()
  dims = shape(this%crss_grid)
  npoints = npoints * dims(1)
  call MPIComputeQualitySSScalar(this%crss_grid, this%crss_old_grid, this%crss_rate_grid, this%crss_rate_old_grid, & 
                npoints, this%dt,this%crss_abs_rel_tol(2), this%crss_abs_rel_tol(1), Q, max_crss_increment)
  call computeMaxAllowedTimeStepQuality(Q, this%dt, "CRSS ",  accept_var, max_crss_increment, max_dt_var)
  accept_solution_flag = accept_solution_flag.and.accept_var
  dt_max = min(dt_max, max_dt_var)

  if (accept_solution_flag.and.i_am_mpi_master) write(*,*) "solution might be ACCEPTED"
  if (.not.accept_solution_flag.and.i_am_mpi_master) write(*,*) "solution will be REJECTED"
  if (i_am_mpi_master) write(*,*) "new allowable dt is ", dt_max
  if (i_am_mpi_master) write(*,*) "*********************************"
  if (i_am_mpi_master) write(*,*) ""

end subroutine

end submodule
