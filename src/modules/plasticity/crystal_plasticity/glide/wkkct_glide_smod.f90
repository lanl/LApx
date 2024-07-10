submodule(glide_mod) wkkct_glide_smod
use test_utils_mod
implicit none

#include "macro_debug.fpp"


contains

module subroutine readMaterialParametersFromFileWKKCTGlide(this, matf_reader)
  use math_constants, only : PI
  use read_from_file_utils, only : file_reader
  use string_module , only : string_array
  use units_conversion_mod, only : eV2Joule, MPa2Pa
  use number_to_string_mod
  implicit none
  class(wkkct_glide), intent(inout) :: this
  class(file_reader), intent(inout) :: matf_reader
  type(string_array) :: dummy_string_array
  real(k_real), dimension(this%n_slip_modes) :: BurgN_temp
  real(k_real), dimension(:), pointer :: mode_dep_parameter => null()
  integer :: prec_id, solute_id
  real(k_real), dimension(:), pointer :: obstacle_type_variable_temp => null()
  real(k_real) :: c0_cd_tmp

  BurgN_temp = this%burgN
  deallocate(this%burgN); nullify(this%burgN)
  call this%convertModeToSlipSystemParameter(BurgN_temp, this%burgN)

  this%rho_mass => this%common_material_parameter_ptr%mass_density
  this%grain_diameter = this%common_material_parameter_ptr%grain_diameter
  !sanity checks
  if (this%rho_mass<=0._k_real) error stop "WKKCTGlide rho_mass <=0."

  call matf_reader%readLineAndCheckStringsAreEqual( "-microstructure-fingerprint-parameters", dummy_string_array)
  call matf_reader%readParameter("number-of-precipitate-types[n_prec_type]", this%n_prec_type)
  if (this%n_prec_type>0) then
  call matf_reader%readVectorParameter("precipitate-number-density[number_density_prec]", this%n_prec_type, this%number_density_prec_ptr)
  call matf_reader%readVectorParameter("initial-pricipitate-size[initial_diameter_prec]", this%n_prec_type, this%initial_diameter_prec_ptr)
  endif
  call matf_reader%readParameter("number-of-solute-types[n_solute_type]", this%n_solute_type)
  if(this%n_solute_type.gt.0) then
    call matf_reader%readParameter("solutes-as-backstress[TRUE/FALSE]", this%solutes_as_backstress)
  else 
    call matf_reader%readLineAndCheckStringsAreEqual("solutes-as-backstress[TRUE/FALSE]", dummy_string_array)
    this%solutes_as_backstress = .FALSE.
  endif
  if(this%n_dislocation_loop_type.gt.0) then
    write(*,*) "this%n_dislocation_loop_type ", this%n_dislocation_loop_type
    call matf_reader%readVectorParameter("dislocation-loops-number-density[1/m^2]", this%n_dislocation_loop_type, this%number_density_disloc_loops_ptr)
    call matf_reader%readVectorParameter("dislocation-loops-size[m]", this%n_dislocation_loop_type, this%initial_diameter_disloc_loops_ptr)
    if (ANY (this%number_density_disloc_loops_ptr .le. 0._k_real)) error stop "WKKCTGlide number_density_disloc_loops_ptr <=0."
    if (ANY (this%initial_diameter_disloc_loops_ptr .le. 0._k_real)) error stop "WKKCTGlide initial_diameter_disloc_loops_ptr <=0."
  else
    call matf_reader%readLineAndCheckStringsAreEqual("dislocation-loops-number-density[1/m^2]", dummy_string_array)
    call matf_reader%readLineAndCheckStringsAreEqual("dislocation-loops-size[m]", dummy_string_array)
  endif

  call matf_reader%readVectorParameter( "initial-cell-interior-density[rho_ci]", this%n_slip_modes, mode_dep_parameter)
  call this%convertModeToSlipSystemParameter(mode_dep_parameter, this%rho0_m)
  deallocate(mode_dep_parameter); nullify(mode_dep_parameter)
  call matf_reader%readParameter( "dislocation-cell-wall-present[TRUE/FALSE]", this%use_dislocation_cell_wall)
  if (this%use_dislocation_cell_wall) then
    call matf_reader%readVectorParameter( "initial-cell-wall-disl-density[rho_cw]", this%n_slip_modes, mode_dep_parameter)
    call this%convertModeToSlipSystemParameter(mode_dep_parameter, this%rho0_cw)
    deallocate(mode_dep_parameter); nullify(mode_dep_parameter)
    if (any(this%rho0_cw<0._k_real)) error stop "WKKCTGlide rho0_cw <0."
  else 
    call matf_reader%readLineAndCheckStringsAreEqual("initial-cell-wall-disl-density[rho_cw]", dummy_string_array)
  endif
  
  
  call matf_reader%readParameter( "standard-deviation-scaling-coeff[etav]", this%etav)
  call matf_reader%readParameter( "evolve-standard-deviation[TRUE/FALSE]", this%evolve_std_dev)
  call matf_reader%readParameter( "preaging-time[s]", this%pre_aging_time)
  call matf_reader%skipEmptyLine()

  !sanity checks
  if (this%n_prec_type>0) then
  if (ANY (this%number_density_prec_ptr<0._k_real)) error stop "WKKCTGlide number_density_prec <=0."
  if (ANY (this%initial_diameter_prec_ptr<=0._k_real) ) error stop "WKKCTGlide initial_diameter_prec <=0."
  endif

  if (any(this%rho0_m<0._k_real)) error stop "WKKCTGlide rho0_m <0."
  

  call matf_reader%readLineAndCheckStringsAreEqual( "-glide-velocity-parameters", dummy_string_array)
  call matf_reader%readVectorParameter("p_exp", this%n_slip_modes, mode_dep_parameter)
  call this%convertModeToSlipSystemParameter(mode_dep_parameter, this%p_exp)
  deallocate(mode_dep_parameter); nullify(mode_dep_parameter)
  call matf_reader%readVectorParameter("q_exp", this%n_slip_modes, mode_dep_parameter)
  call this%convertModeToSlipSystemParameter(mode_dep_parameter, this%q_exp)
  deallocate(mode_dep_parameter); nullify(mode_dep_parameter)
  call matf_reader%readParameter( "n-exponent-damage[unitless]", this%n_exp_damage)
  call matf_reader%readParameter("K-hall-petch[1/sqrt(m)]", this%k_hall_petch)
  call matf_reader%readVectorParameter("crss-superposition-exponent[n_crss]", this%n_slip_modes, mode_dep_parameter)
  call this%convertModeToSlipSystemParameter(mode_dep_parameter, this%n_crss)
  deallocate(mode_dep_parameter); nullify(mode_dep_parameter)
  call matf_reader%readVectorParameter("thermal-activation-energy-dislocation[deltaG0_disl]", this%n_slip_modes, mode_dep_parameter)
  call this%convertModeToSlipSystemParameter(mode_dep_parameter*eV2Joule, this%deltaG0_disl)
  deallocate(mode_dep_parameter); nullify(mode_dep_parameter)
  call matf_reader%readVectorParameter("entropy-factor[chi_e]", this%n_slip_modes, mode_dep_parameter)
  call this%convertModeToSlipSystemParameter(mode_dep_parameter, this%chi_e)
  call matf_reader%readParameter("use-jog-screw-model[TRUE/FALSE]", this%use_jog_screw)
  deallocate(mode_dep_parameter); nullify(mode_dep_parameter)
  !
  call matf_reader%skipEmptyLine()
  if (this%n_prec_type>0) then
  allocate(this%deltaG0_prec_ptr(this%n_prec_type, this%n_ss))
  allocate(this%beta_prec_ptr(this%n_prec_type, this%n_ss))
  allocate(this%alpha_prec_ptr(this%n_prec_type, this%n_ss))
  allocate(this%freq_prec_ptr(this%n_prec_type, this%n_ss))
  call matf_reader%readLineAndCheckStringsAreEqual( "-precipitate-related-parameters", dummy_string_array)
  do prec_id = 1,this%n_prec_type
    call matf_reader%readLineAndCheckStringsAreEqual( "precipitate-type-"//trim(adjustl(int2string(prec_id))), dummy_string_array)
    call matf_reader%readVectorParameter("p"//trim(adjustl(int2string(prec_id)))//"-thermal-activation-energy-precipitate[deltaG0_prec]", this%n_slip_modes, mode_dep_parameter)
    call this%convertModeToSlipSystemParameter(mode_dep_parameter*eV2Joule, obstacle_type_variable_temp)
    this%deltaG0_prec_ptr(prec_id,:) = obstacle_type_variable_temp
    deallocate(mode_dep_parameter, obstacle_type_variable_temp); nullify(mode_dep_parameter, obstacle_type_variable_temp)
    !
    call matf_reader%readVectorParameter("p"//trim(adjustl(int2string(prec_id)))//"-precipitate-interspacing-coeff[beta_prec]", this%n_slip_modes, mode_dep_parameter)
    call this%convertModeToSlipSystemParameter(mode_dep_parameter, obstacle_type_variable_temp)
    this%beta_prec_ptr(prec_id,:) = obstacle_type_variable_temp
    deallocate(mode_dep_parameter, obstacle_type_variable_temp); nullify(mode_dep_parameter, obstacle_type_variable_temp)
    !
    call matf_reader%readVectorParameter("p"//trim(adjustl(int2string(prec_id)))//"-hardening-coeff-precipitate[alpha_prec]", this%n_slip_modes, mode_dep_parameter)
    call this%convertModeToSlipSystemParameter(mode_dep_parameter, obstacle_type_variable_temp)
    this%alpha_prec_ptr(prec_id,:) = obstacle_type_variable_temp
    deallocate(mode_dep_parameter, obstacle_type_variable_temp); nullify(mode_dep_parameter, obstacle_type_variable_temp)
    !
    call matf_reader%readVectorParameter("p"//trim(adjustl(int2string(prec_id)))//"-attack-frequency-precipitate[freq_prec]", this%n_slip_modes, mode_dep_parameter)
    call this%convertModeToSlipSystemParameter(mode_dep_parameter, obstacle_type_variable_temp)
    this%freq_prec_ptr(prec_id,:) = obstacle_type_variable_temp
    deallocate(mode_dep_parameter, obstacle_type_variable_temp); nullify(mode_dep_parameter, obstacle_type_variable_temp)
  end do
  call matf_reader%skipEmptyLine()
  endif

  call matf_reader%readLineAndCheckStringsAreEqual( "-strength-parameter-tau0", dummy_string_array)
  call this%tau0_linear_interp%readFromFile(matf_reader, this%n_slip_modes, this%n_ss_per_mode, "tau0", "MPa")
  call matf_reader%skipEmptyLine()

  ! Keep the following lines for a while, so that we can easily go back to linear interpoltation instead of Labusch equation.
  ! call matf_reader%readLineAndCheckStringsAreEqual( "-strength-parameter-tau-Labusch", dummy_string_array)
  ! call this%tau_labusch_linear_interp%readFromFile(file_id, this%n_slip_modes, this%n_ss_per_mode, "tau-labusch", "MPa")
  ! call matf_reader%skipEmptyLine()
  if (this%n_solute_type>0) then
  call matf_reader%readLineAndCheckStringsAreEqual( "-solute-related-parameters", dummy_string_array)
  !TODO need to add a new subroutine to take of array inputs with obstcles and slip-systems
  allocate(this%solute_concentration_ptr(this%n_solute_type))
  allocate(this%K_Labusch_ptr(this%n_solute_type,this%n_ss))
  allocate(this%exponent_Labusch_ptr(this%n_solute_type,this%n_ss))
  allocate(this%DW_ptr(this%n_solute_type,this%n_ss))
  allocate(this%DHc_ptr(this%n_solute_type,this%n_ss))
  allocate(this%m_cd_ptr(this%n_solute_type,this%n_ss))
  allocate(this%freq_cd_ptr(this%n_solute_type,this%n_ss))
  allocate(this%norm_wbar_ptr(this%n_solute_type,this%n_ss))
  allocate(this%disl_core_width_ptr(this%n_solute_type,this%n_ss))
  allocate( this%N_atoms_disl_core_ptr(this%n_solute_type,this%n_ss))
  allocate(this%phi_cd_ptr(this%n_solute_type,this%n_ss))
  allocate(this%alpha_cd_ptr(this%n_solute_type,this%n_ss))
  allocate(this%chi_cd_ptr(this%n_solute_type,this%n_ss))
  do solute_id = 1,this%n_solute_type
    call matf_reader%readLineAndCheckStringsAreEqual( "solute-type-"//trim(adjustl(int2string(solute_id))), dummy_string_array)
    call matf_reader%readParameter( "s"//trim(adjustl(int2string(solute_id)))//"-initial-solute-concentration[c0_cd]", c0_cd_tmp)
    this%solute_concentration_ptr(solute_id) = c0_cd_tmp

    call matf_reader%readVectorParameter("s"//trim(adjustl(int2string(solute_id)))//"-Labusch_coefficient[K_Lab]", this%n_slip_modes, mode_dep_parameter)
    call this%convertModeToSlipSystemParameter(mode_dep_parameter, obstacle_type_variable_temp)
    this%K_Labusch_ptr(solute_id,:) = obstacle_type_variable_temp
    deallocate(mode_dep_parameter, obstacle_type_variable_temp); nullify(mode_dep_parameter, obstacle_type_variable_temp)

    call matf_reader%readVectorParameter("s"//trim(adjustl(int2string(solute_id)))//"-solute_exponent[n_solute]", this%n_slip_modes, mode_dep_parameter)
    call this%convertModeToSlipSystemParameter(mode_dep_parameter, obstacle_type_variable_temp)
    this%exponent_Labusch_ptr(solute_id,:) = obstacle_type_variable_temp
    deallocate(mode_dep_parameter, obstacle_type_variable_temp); nullify(mode_dep_parameter, obstacle_type_variable_temp)

    call matf_reader%readVectorParameter("s"//trim(adjustl(int2string(solute_id)))//"-binding-energy-difference[DW]", this%n_slip_modes, mode_dep_parameter)
    call this%convertModeToSlipSystemParameter(mode_dep_parameter, obstacle_type_variable_temp)
    this%DW_ptr(solute_id,:) = obstacle_type_variable_temp*eV2Joule
    deallocate(mode_dep_parameter, obstacle_type_variable_temp); nullify(mode_dep_parameter, obstacle_type_variable_temp)

    call matf_reader%readVectorParameter("s"//trim(adjustl(int2string(solute_id)))//"-solute-activation-enthalpy[DHc]", this%n_slip_modes, mode_dep_parameter)
    call this%convertModeToSlipSystemParameter(mode_dep_parameter, obstacle_type_variable_temp)
    this%DHc_ptr(solute_id,:) = obstacle_type_variable_temp*eV2Joule
    deallocate(mode_dep_parameter, obstacle_type_variable_temp); nullify(mode_dep_parameter, obstacle_type_variable_temp)

    call matf_reader%readVectorParameter("s"//trim(adjustl(int2string(solute_id)))//"-num_neighbors-core-diffusion[m_cd]", this%n_slip_modes, mode_dep_parameter)
    call this%convertModeToSlipSystemParameter(mode_dep_parameter, obstacle_type_variable_temp)
    this%m_cd_ptr(solute_id,:) = obstacle_type_variable_temp
    deallocate(mode_dep_parameter, obstacle_type_variable_temp); nullify(mode_dep_parameter, obstacle_type_variable_temp)

    call matf_reader%readVectorParameter("s"//trim(adjustl(int2string(solute_id)))//"-frequency-core-diffusion[freq_cd]", this%n_slip_modes, mode_dep_parameter)
    call this%convertModeToSlipSystemParameter(mode_dep_parameter, obstacle_type_variable_temp)
    this%freq_cd_ptr(solute_id,:) = obstacle_type_variable_temp
    deallocate(mode_dep_parameter, obstacle_type_variable_temp); nullify(mode_dep_parameter, obstacle_type_variable_temp)

    call matf_reader%readVectorParameter("s"//trim(adjustl(int2string(solute_id)))//"-normalized-core-width[norm_wbar]", this%n_slip_modes, mode_dep_parameter)
    call this%convertModeToSlipSystemParameter(mode_dep_parameter, obstacle_type_variable_temp)
    this%disl_core_width_ptr(solute_id,:) = this%burgN*obstacle_type_variable_temp
    deallocate(mode_dep_parameter, obstacle_type_variable_temp); nullify(mode_dep_parameter, obstacle_type_variable_temp)

    !! compute the number of atoms within a dislocation core
    this%N_atoms_disl_core_ptr(solute_id,:) = 2._k_real * this%disl_core_width_ptr(solute_id,:)  /(sqrt(3._k_real) * this%burgN**2)

    call matf_reader%readVectorParameter("s"//trim(adjustl(int2string(solute_id)))//"-exponent-core-diffusion[phi_cd]", this%n_slip_modes, mode_dep_parameter)
    call this%convertModeToSlipSystemParameter(mode_dep_parameter, obstacle_type_variable_temp)
    this%phi_cd_ptr(solute_id,:) = obstacle_type_variable_temp
    deallocate(mode_dep_parameter, obstacle_type_variable_temp); nullify(mode_dep_parameter, obstacle_type_variable_temp)

    call matf_reader%readVectorParameter("s"//trim(adjustl(int2string(solute_id)))//"-coeff-energy-variation-along-the-core[alpha_cd]", this%n_slip_modes, mode_dep_parameter)
    call this%convertModeToSlipSystemParameter(mode_dep_parameter, obstacle_type_variable_temp)
    this%alpha_cd_ptr(solute_id,:) = obstacle_type_variable_temp
    deallocate(mode_dep_parameter, obstacle_type_variable_temp); nullify(mode_dep_parameter, obstacle_type_variable_temp)

    call matf_reader%readVectorParameter("s"//trim(adjustl(int2string(solute_id)))//"-prefactor-binding-energy-junction-strength[chi_cd]", this%n_slip_modes, mode_dep_parameter)
    call this%convertModeToSlipSystemParameter(mode_dep_parameter, obstacle_type_variable_temp)
    this%chi_cd_ptr(solute_id,:) = obstacle_type_variable_temp
    deallocate(mode_dep_parameter, obstacle_type_variable_temp); nullify(mode_dep_parameter, obstacle_type_variable_temp)
  end do
  if (ANY (this%solute_concentration_ptr<0._k_real) ) error stop "WKKCTGlide solute_concentration <=0."
  call matf_reader%skipEmptyLine()
  endif


  call matf_reader%readLineAndCheckStringsAreEqual( "-dislocation-evolution-parameters", dummy_string_array)
  call matf_reader%readVectorParameter("coeff-disl-generation[k_disl_gen]", this%n_slip_modes, mode_dep_parameter)
  call this%convertModeToSlipSystemParameter(mode_dep_parameter, this%k_disl_gen)
  deallocate(mode_dep_parameter); nullify(mode_dep_parameter)
  call matf_reader%readVectorParameter("reference-strain-rate-dyn-recovery[edot0_dynrec]", this%n_slip_modes, mode_dep_parameter)
  call this%convertModeToSlipSystemParameter(mode_dep_parameter, this%edot0_dynrec)
  deallocate(mode_dep_parameter); nullify(mode_dep_parameter)
  call matf_reader%readVectorParameter("chi-dyn-recovery[chi_dynrec]", this%n_slip_modes, mode_dep_parameter)
  call this%convertModeToSlipSystemParameter(mode_dep_parameter, this%chi_dynrec)
  deallocate(mode_dep_parameter); nullify(mode_dep_parameter)
  call matf_reader%readVectorParameter("D-dyn-recovery[D_dynrec]", this%n_slip_modes, mode_dep_parameter)
  call this%convertModeToSlipSystemParameter(mode_dep_parameter, this%D_dynrec)
  deallocate(mode_dep_parameter); nullify(mode_dep_parameter)
  call matf_reader%readVectorParameter("g-dyn-recovery[g_dynrec]", this%n_slip_modes, mode_dep_parameter)
  call this%convertModeToSlipSystemParameter(mode_dep_parameter, this%g_dynrec)
  deallocate(mode_dep_parameter); nullify(mode_dep_parameter)

  call matf_reader%readParameter("use-dislocation-static-recovery[TRUE/FALSE]", this%use_static_recovery) 
  if (this%use_static_recovery) then
    call matf_reader%readVectorParameter("k1-static-recovery[bulk,GB]", 2, this%k1_static_recov)
    call matf_reader%readVectorParameter("k2-static-recovery[bulk,GB]", 2, this%k2_static_recov)
  else 
    call matf_reader%readLineAndCheckStringsAreEqual("k1-static-recovery[bulk,GB]", dummy_string_array)
    call matf_reader%readLineAndCheckStringsAreEqual("k2-static-recovery[bulk,GB]", dummy_string_array)
    allocate(this%k1_static_recov(this%n_ss)) 
    this%k1_static_recov = 0._k_real
    allocate(this%k2_static_recov(this%n_ss)) 
    this%k2_static_recov = 0._k_real
  endif

  call matf_reader%readVectorParameter("coeff-trapping-rate-cell-wall[k_disl_trap]", this%n_slip_modes, mode_dep_parameter)
  call this%convertModeToSlipSystemParameter(mode_dep_parameter, this%k_disl_trap)
  deallocate(mode_dep_parameter); nullify(mode_dep_parameter)
  call matf_reader%readVectorParameter("coeff-annyhilation-cell-wall[k_cell_wall_ann]", this%n_slip_modes, mode_dep_parameter)
  call this%convertModeToSlipSystemParameter(mode_dep_parameter, this%k_cell_wall_ann)
  deallocate(mode_dep_parameter); nullify(mode_dep_parameter)
  call matf_reader%skipEmptyLine()

  call matf_reader%readLineAndCheckStringsAreEqual( "-time-march-tolerances", dummy_string_array)
  call matf_reader%readParameter( "max-allowed-dd-percent-increment[unitless]", this%max_allowed_dd_percent_increment)
  call matf_reader%readParameter( "max-allowed-tau-sol-drag-increment[MPa]", this%max_allowed_tau_sol_drag_increment)
  call matf_reader%readParameter( "max-allowed-tau-eline-increment[MPa]", this%max_allowed_tau_eline_increment)
  if (this%max_allowed_dd_percent_increment.le.0._k_real) error stop "max-allowed-dd-percent-increment[unitless] =< 0"
  if (this%max_allowed_tau_sol_drag_increment.le.0._k_real) error stop "max-allowed-tau-sol-drag-increment[MPa] =< 0"
  if (this%max_allowed_tau_eline_increment.le.0._k_real) error stop "max-allowed-tau-eline-increment[MPa] =< 0"

end subroutine

module subroutine initParametersWKKCTGlide(this, phase_id, common_material_parameter_ptr, &
                                            use_damage, n_gauss, n_std_dev_gauss_integration, &
                                            elasticity_obj, crystal_paremeters_ptr)
  use stiffness_base_mod, only : stiffness_base
  use cp_base_mod, only : crystal_paremeters_type 
  use math_constants, only : PI
  use units_conversion_mod, only : MPa2Pa
  implicit none
  class(wkkct_glide), intent(inout) :: this
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

  allocate(this%burgN(this%n_slip_modes))
  this%burgN = crystal_paremeters_ptr%burgVectorL

  this%n_dislocation_loop_type = crystal_paremeters_ptr%number_disloc_loop_type
  if(this%n_dislocation_loop_type.gt.0) then
    allocate(this%loop_crystallography_factor_ptr( crystal_paremeters_ptr%n_ss_total, this%n_dislocation_loop_type))
    this%loop_crystallography_factor_ptr=crystal_paremeters_ptr%loop_crystallography_factor
    write(*,*)'this%loop_crystallography_factor_ptr', this%loop_crystallography_factor_ptr
  endif

end subroutine

module subroutine  addFieldVariablesToGridWKKCTGlide(this)
  use grid_data_var_type_mod
  implicit none
  class(wkkct_glide), intent(inout) :: this

  call this%cp_glide_base%addFieldVariablesToGrid()
  associate (all_grid_data_vars => this%grid_data)

  call all_grid_data_vars%addVar("crss_glide_gauss", ss_generic_vector, & 
    additional_var_dimensions=(/this%n_gauss/), stateful_level=2)
  call all_grid_data_vars%addVar("crss", ss_scalar, stateful_level=2)
  call all_grid_data_vars%addVar("tau_line_tension", ss_scalar, stateful_level=3)
  call all_grid_data_vars%addVar("tau_hall_petch", ss_scalar, stateful_level=2)

  if (this%n_solute_type>0) then
    call all_grid_data_vars%addVar("tau_solute_drag_gauss", ss_generic_vector, &
                                    additional_var_dimensions=(/this%n_gauss/), stateful_level=3)
    call all_grid_data_vars%addVar("tau_m", ss_scalar, stateful_level=2)
    call all_grid_data_vars%addVar("time_solute_aging", ss_generic_vector, &
                                    additional_var_dimensions=(/this%n_gauss/), stateful_level=2)
    call all_grid_data_vars%addVar("time_solute_aging_avg", ss_scalar, stateful_level=2)
    call all_grid_data_vars%addVar("tau_labusch", ss_scalar, stateful_level=2)
  endif

  call all_grid_data_vars%addVar("lambda_s", ss_scalar, stateful_level=2)
  call all_grid_data_vars%addVar("lambda_dislocation", ss_scalar, stateful_level=2)
  call all_grid_data_vars%addVar("lambda_hall_petch", ss_scalar, stateful_level=2)


  call all_grid_data_vars%addVar("prob_d", ss_scalar, stateful_level=2)
  if (this%n_prec_type>0) then
    call all_grid_data_vars%addVar("tau_precipitate_drag", ss_scalar, stateful_level=2)
    call all_grid_data_vars%addVar("prob_p", ss_generic_vector, &
                                 additional_var_dimensions=(/this%n_prec_type/), stateful_level=2)
    call all_grid_data_vars%addVar("diameter_prec", ss_generic_vector, &
                                additional_var_dimensions=(/this%n_prec_type/), stateful_level=2)
    call all_grid_data_vars%addVar("time_prec_climb", ss_generic_vector, &
                                additional_var_dimensions=(/this%n_gauss/), stateful_level=2)
    call all_grid_data_vars%addVar("lambda_precipitate", ss_scalar, stateful_level=2)
  endif
  call all_grid_data_vars%addVar("rho_m", ss_scalar, stateful_level=2)
  call all_grid_data_vars%addVar("rho_m_dot", ss_scalar, stateful_level=2)
  call all_grid_data_vars%addVar("delta_rho_m_generation", ss_scalar, stateful_level=2)
  call all_grid_data_vars%addVar("delta_rho_m_dyn_recovery", ss_scalar, stateful_level=2)

  if (this%use_static_recovery) then
    call all_grid_data_vars%addVar("delta_rho_m_static_recovery", ss_scalar, stateful_level=2)
  endif

  if (this%use_dislocation_cell_wall) then
    call all_grid_data_vars%addVar("tau_cell_wall", ss_scalar, stateful_level=2)
    call all_grid_data_vars%addVar("rho_cw", ss_scalar, stateful_level=2)
    call all_grid_data_vars%addVar("delta_rho_m_trapped_cw", ss_scalar, stateful_level=2)
    call all_grid_data_vars%addVar("delta_rho_cw_annihilation", ss_scalar, stateful_level=2)
    call all_grid_data_vars%addVar("climb_velocity_CellWall", ss_generic_vector, &
    additional_var_dimensions=(/this%n_gauss/), stateful_level=2)
  endif
  call all_grid_data_vars%addVar("sub_grain_size", ss_scalar, stateful_level=2)

  if (this%n_dislocation_loop_type>0) then
    call all_grid_data_vars%addVar("tau_dislocation_loop", ss_scalar, stateful_level=2)
  endif

  ! this model does not calculate the climb velocity. however we request teh grid
  ! to add the "climb_velocity" to the grid so that it is available for the updating
  ! diameter_prec. TODO for optimization we might want to ask for the "climb_velocity" variable
  ! only if needed
  call all_grid_data_vars%addVar("climb_velocity", ss_generic_vector, &
                                 additional_var_dimensions=(/this%n_gauss/), stateful_level=2)

end associate
end subroutine

module subroutine initGridPointersWKKCTGlide(this)
  implicit none
  class(wkkct_glide), intent(inout) :: this

  call this%cp_glide_base%initGridPointers()
  call this%grid_data%getSSGenericVectorDataPointerByName("crss_glide_gauss", this%crss_gauss_grid)
  call this%grid_data%getSSScalarDataPointerByName("crss", this%crss_grid)
  call this%grid_data%getSSScalarDataPointerByName("tau_line_tension", this%tau_line_tension_grid)
  call this%grid_data%getSSScalarDataPointerByNameOld("tau_line_tension", this%tau_line_tension_grid_old)
  call this%grid_data%getSSScalarDataPointerByNameOlder("tau_line_tension", this%tau_line_tension_grid_older)
  call this%grid_data%getSSScalarDataPointerByName("tau_hall_petch", this%tau_hall_petch_grid)
  if (this%n_solute_type>0) then
    call this%grid_data%getSSGenericVectorDataPointerByName("tau_solute_drag_gauss", this%tau_solute_drag_gauss_grid)
    call this%grid_data%getSSGenericVectorDataPointerByNameOld("tau_solute_drag_gauss", this%tau_solute_drag_gauss_grid_old)
    call this%grid_data%getSSGenericVectorDataPointerByNameOlder("tau_solute_drag_gauss", this%tau_solute_drag_gauss_grid_older)
    call this%grid_data%getSSScalarDataPointerByName("tau_m", this%tau_m_grid)
    call this%grid_data%getSSScalarDataPointerByName("tau_labusch", this%tau_labusch_grid)
    call this%grid_data%getSSGenericVectorDataPointerByName("time_solute_aging", this%time_solute_aging_grid)
    call this%grid_data%getSSScalarDataPointerByName("time_solute_aging_avg", this%time_solute_aging_avg_grid)
    call this%grid_data%getSSGenericVectorDataPointerByNameOld("time_solute_aging", this%time_solute_aging_grid_old)
  endif
  call this%grid_data%getSSScalarDataPointerByName("lambda_s", this%lambda_s_grid)
  call this%grid_data%getSSScalarDataPointerByName("lambda_dislocation", this%lambda_dislocation_grid)
  call this%grid_data%getSSScalarDataPointerByName("lambda_hall_petch", this%lambda_hall_petch_grid)

  call this%grid_data%getSSScalarDataPointerByName("prob_d", this%prob_d_grid)
  if (this%n_prec_type>0) then
    call this%grid_data%getSSScalarDataPointerByName("tau_precipitate_drag", this%tau_precipitate_drag_grid)
    call this%grid_data%getSSGenericVectorDataPointerByName("prob_p", this%prob_p_grid)
    call this%grid_data%getSSGenericVectorDataPointerByName("diameter_prec", this%diameter_prec_grid)
    call this%grid_data%getSSGenericVectorDataPointerByName("time_prec_climb", this%time_prec_climb_grid)
    call this%grid_data%getSSGenericVectorDataPointerByNameOld("time_prec_climb", this%time_prec_climb_grid_old)
    call this%grid_data%getSSScalarDataPointerByName("lambda_precipitate", this%lambda_precipitate_grid)
  endif
  call this%grid_data%getSSScalarDataPointerByName("rho_m", this%rho_m_grid)
  call this%grid_data%getSSScalarDataPointerByNameOld("rho_m", this%rho_m_grid_old)
  call this%grid_data%getSSScalarDataPointerByName("rho_m_dot", this%rho_m_dot_grid)
  call this%grid_data%getSSScalarDataPointerByNameOld("rho_m_dot", this%rho_m_dot_grid_old)
  call this%grid_data%getSSScalarDataPointerByName("delta_rho_m_generation", this%delta_rho_m_generation_grid)
  call this%grid_data%getSSScalarDataPointerByName("delta_rho_m_dyn_recovery", this%delta_rho_m_dyn_recovery_grid)

  if (this%use_static_recovery) then
    call this%grid_data%getSSScalarDataPointerByName("delta_rho_m_static_recovery", this%delta_rho_m_static_recovery_grid)
  endif

  if (this%use_dislocation_cell_wall) then
    call this%grid_data%getSSScalarDataPointerByName("tau_cell_wall", this%tau_cell_wall_grid)
    call this%grid_data%getSSScalarDataPointerByName("rho_cw", this%rho_cw_grid)
    call this%grid_data%getSSScalarDataPointerByNameOld("rho_cw", this%rho_cw_grid_old)
    call this%grid_data%getSSScalarDataPointerByName("delta_rho_m_trapped_cw", this%delta_rho_m_trapped_cw_grid)
    call this%grid_data%getSSScalarDataPointerByName("delta_rho_cw_annihilation", this%delta_rho_cw_annihilation_grid)
    call this%grid_data%getSSGenericVectorDataPointerByName("climb_velocity_CellWall", this%climb_velocity_CellWall_grid)
  endif
  if (this%n_dislocation_loop_type>0) then
    call this%grid_data%getSSScalarDataPointerByName("tau_dislocation_loop", this%tau_disl_loop_grid)
  endif
  call this%grid_data%getSSGenericVectorDataPointerByName("climb_velocity", this%climb_velocity_grid)
  call this%grid_data%getSSScalarDataPointerByName("sub_grain_size", this%sub_grain_size_grid)

  ! this is the RVE average equivalent strain rate, and its not a proper grid pointer
  call this%sim_all_macro_data%sim_macro_field_averages%getAverageEqStrainRatePointer(this%avg_eq_strain_rate_ptr)
end subroutine

module subroutine setPointDataWKKCTGlide(this, ix,iy,iz)
  implicit none
  class(wkkct_glide), intent(inout) :: this
  integer, intent(in) :: ix,iy,iz

  associate (nss => this%n_ss, &
             ngauss => this%n_gauss, &
             nprec => this%n_prec_type)
  !calling the parent class to do its job
  call this%cp_glide_base%setPointData(ix,iy,iz)
  this%crss_gauss_ptr => this%crss_gauss_grid(1:ngauss, 1:nss,ix,iy,iz)
  this%crss_ptr => this%crss_grid(1:nss,ix,iy,iz)
  this%tau_line_tension_ptr => this%tau_line_tension_grid(1:nss,ix,iy,iz)
  this%tau_line_tension_ptr_old => this%tau_line_tension_grid_old(1:nss,ix,iy,iz)
  this%tau_hall_petch_ptr => this%tau_hall_petch_grid(1:nss,ix,iy,iz)
  if (this%n_solute_type>0) then
    this%tau_solute_drag_gauss_ptr => this%tau_solute_drag_gauss_grid(1:ngauss,1:nss,ix,iy,iz)
    this%tau_solute_drag_gauss_ptr_old => this%tau_solute_drag_gauss_grid_old(1:ngauss,1:nss,ix,iy,iz)
    this%tau_m_ptr => this%tau_m_grid(1:nss,ix,iy,iz)
    this%tau_labusch_ptr => this%tau_labusch_grid(1:nss,ix,iy,iz)
    this%time_solute_aging_ptr => this%time_solute_aging_grid(1:ngauss,1:nss,ix,iy,iz)
    this%time_solute_aging_avg_ptr => this%time_solute_aging_avg_grid(1:nss,ix,iy,iz)
    this%time_solute_aging_ptr_old => this%time_solute_aging_grid_old(1:ngauss,1:nss,ix,iy,iz)
  else
    if (.not.(associated(this%tau_solute_drag_gauss_ptr))) allocate(this%tau_solute_drag_gauss_ptr(ngauss, nss))
    this%tau_solute_drag_gauss_ptr = 0._k_real
    if (.not.(associated(this%tau_solute_drag_gauss_ptr_old))) allocate(this%tau_solute_drag_gauss_ptr_old(ngauss, nss))
    this%tau_solute_drag_gauss_ptr_old = 0._k_real
    if (.not.(associated(this%tau_labusch_ptr))) allocate(this%tau_labusch_ptr(nss))
    this%tau_labusch_ptr = 0._k_real
    if (.not.(associated(this%time_solute_aging_ptr))) allocate(this%time_solute_aging_ptr(ngauss,nss))
    this%time_solute_aging_ptr = 0._k_real
    if (.not.(associated(this%time_solute_aging_ptr_old))) allocate(this%time_solute_aging_ptr_old(ngauss,nss))
    this%time_solute_aging_ptr_old = 0._k_real
    if (.not.(associated(this%tau_m_ptr))) allocate(this%tau_m_ptr(nss))
    this%tau_m_ptr = 0._k_real
  endif
  this%lambda_s_ptr => this%lambda_s_grid(1:nss,ix,iy,iz)
  this%lambda_dislocation_ptr => this%lambda_dislocation_grid(1:nss,ix,iy,iz)
  this%lambda_hall_petch_ptr => this%lambda_hall_petch_grid(1:nss,ix,iy,iz)
  

  this%prob_d_ptr => this%prob_d_grid(1:nss,ix,iy,iz)
  if (this%n_prec_type>0) then
    this%tau_precipitate_drag_ptr => this%tau_precipitate_drag_grid(1:nss,ix,iy,iz)
    this%prob_p_ptr => this%prob_p_grid(1:nprec,1:nss,ix,iy,iz)
    this%diameter_prec_ptr => this%diameter_prec_grid(1:nprec,1:nss,ix,iy,iz)
    this%time_prec_climb_ptr => this%time_prec_climb_grid(1:ngauss,1:nss,ix,iy,iz)
    this%time_prec_climb_ptr_old => this%time_prec_climb_grid_old(1:ngauss,1:nss,ix,iy,iz)
    this%lambda_precipitate_ptr => this%lambda_precipitate_grid(1:nss,ix,iy,iz)
  else
    if (.not.(associated(this%tau_precipitate_drag_ptr))) allocate(this%tau_precipitate_drag_ptr(nss))
    this%tau_precipitate_drag_ptr = 0._k_real
    if (.not.(associated(this%prob_p_ptr))) allocate(this%prob_p_ptr(1,nss))
    this%prob_p_ptr = 0._k_real
    if (.not.(associated(this%diameter_prec_ptr))) allocate(this%diameter_prec_ptr(1,nss))
    this%diameter_prec_ptr = 0._k_real
    if (.not.(associated(this%time_prec_climb_ptr))) allocate(this%time_prec_climb_ptr(ngauss,nss))
    this%time_prec_climb_ptr = 0._k_real
    if (.not.(associated(this%time_prec_climb_ptr_old))) allocate(this%time_prec_climb_ptr_old(ngauss,nss))
    this%time_prec_climb_ptr_old = 0._k_real
    if (.not.(associated(this%lambda_precipitate_ptr))) allocate(this%lambda_precipitate_ptr(nss))
    this%lambda_precipitate_ptr = 0._k_real
  endif

  if (this%n_dislocation_loop_type > 0) then
    this%tau_disl_loop_ptr => this%tau_disl_loop_grid(1:nss,ix,iy,iz)
  else 
    if (.not.(associated(this%tau_disl_loop_ptr))) allocate(this%tau_disl_loop_ptr(nss))
    this%tau_disl_loop_ptr = 0._k_real
  endif

  this%rho_m_ptr => this%rho_m_grid(1:nss,ix,iy,iz)
  this%rho_m_dot_ptr => this%rho_m_dot_grid(1:nss,ix,iy,iz)
  this%delta_rho_m_generation_ptr => this%delta_rho_m_generation_grid(1:nss,ix,iy,iz)
  this%delta_rho_m_dyn_recovery_ptr => this%delta_rho_m_dyn_recovery_grid(1:nss,ix,iy,iz)

  if (this%use_dislocation_cell_wall) then
    this%tau_cell_wall_ptr =>this%tau_cell_wall_grid(1:nss,ix,iy,iz)
    this%rho_cw_ptr => this%rho_cw_grid(1:nss,ix,iy,iz)
    this%climb_velocity_CellWall_ptr => this%climb_velocity_CellWall_grid(1:ngauss,1:nss,ix,iy,iz)
    this%sub_grain_size_ptr => this%sub_grain_size_grid(1:nss,ix,iy,iz)
    this%delta_rho_m_trapped_cw_ptr => this%delta_rho_m_trapped_cw_grid(1:nss,ix,iy,iz)
    this%delta_rho_cw_annihilation_ptr => this%delta_rho_cw_annihilation_grid(1:nss,ix,iy,iz)
  else 
    if (.not.(associated(this%tau_cell_wall_ptr))) then 
      allocate(this%tau_cell_wall_ptr(nss))
      this%tau_cell_wall_ptr = 0._k_real
    endif
    if (.not.(associated(this%rho_cw_ptr))) then 
      allocate(this%rho_cw_ptr(nss))
      this%rho_cw_ptr = 0._k_real
    endif
    if (.not.(associated(this%sub_grain_size_ptr))) then
       allocate(this%sub_grain_size_ptr(nss))
      this%sub_grain_size_ptr = this%grain_diameter
    endif
    if (.not.(associated(this%climb_velocity_CellWall_ptr))) then
       allocate(this%climb_velocity_CellWall_ptr(ngauss, nss))
       this%climb_velocity_CellWall_ptr = 0._k_real
    endif
    if (.not.(associated(this%delta_rho_m_trapped_cw_ptr))) then
      allocate(this%delta_rho_m_trapped_cw_ptr(nss))
      this%delta_rho_m_trapped_cw_ptr = 0._k_real
    endif
    if (.not.(associated(this%delta_rho_cw_annihilation_ptr))) then
      allocate(this%delta_rho_cw_annihilation_ptr(nss))
      this%delta_rho_cw_annihilation_ptr = 0._k_real
    endif
  endif


  if (this%use_static_recovery) then
    this%delta_rho_m_static_recovery_ptr => this%delta_rho_m_static_recovery_grid(1:nss,ix,iy,iz)
  else 
    if (.not.(associated(this%delta_rho_m_static_recovery_ptr))) then 
      allocate(this%delta_rho_m_static_recovery_ptr(nss))
      this%delta_rho_m_static_recovery_ptr = 0._k_real
    endif
  endif

  this%climb_velocity_ptr => this%climb_velocity_grid(1:ngauss,1:nss,ix,iy,iz)
  
  end associate

  if (.not.associated(this%shear_mod_phase_ptr)) then
    write(*,*) "this%shear_mod_phase_ptr not asssociated"
  endif
  if (.not.associated(this%poisson_phase_ptr)) then
    write(*,*) "this%poisson_phase_ptr not asssociated"
  endif
  if (.not.associated(this%porosity_ptr)) then
    write(*,*) "this%porosity_ptr not asssociated"
  endif

  this%shear_mod = this%shear_mod_phase_ptr
  this%poisson = this%poisson_phase_ptr
end subroutine

module subroutine initStateVariablesAtMaterialPointWKKCTGlide(this)
  use test_utils_mod
  implicit none
  class(wkkct_glide), intent(inout) :: this
  real(k_real), dimension(this%n_ss, this%n_ss) :: bar_alpha_ss_prime
  integer :: prec_id

  this%rho_m_ptr = this%rho0_m
  if (this%use_dislocation_cell_wall) then
    this%rho_cw_ptr = this%rho0_cw
  else 
    this%rho_cw_ptr = 0._k_real
  endif
  this%ss_idx = 1
  this%igauss = 1

  this%time_solute_aging_ptr = this%pre_aging_time ! we need to add the preaging time to the time
  this%time_prec_climb_ptr = this%time
  do prec_id = 1, this%n_prec_type
    this%diameter_prec_ptr(prec_id,:) = this%initial_diameter_prec_ptr(prec_id)
  end do
  
  call this%computeStdDeviation(this%std_dev_ptr)
  call this%computeResolvedStressAndReslovedStressDStress(this%stress6)
  associate (ss_idx => this%ss_idx)
    do ss_idx=1,this%n_ss
      call this%getGaussianValueAndIntegrationWeights(this%rss(ss_idx), this%std_dev_ptr, this%rss_gaussian, this%total_integration_weight_ptr(:,this%ss_idx))
      if (this%n_solute_type> 0) this%time_solute_aging_avg_ptr(ss_idx) = sum(this%time_solute_aging_ptr(:,ss_idx)*this%total_integration_weight_ptr(:,this%ss_idx))
    enddo
  end associate
  call this%computeEffectiveHardeningMatrix(bar_alpha_ss_prime)

  

  call this%computeEffectiveDislocationInterspacingAndProbability(bar_alpha_ss_prime, this%h_alpha_beta, this%lambda_s_ptr, this%prob_d_ptr, this%prob_p_ptr)
  if (this%use_dislocation_cell_wall) then
    call this%computeSubGrainSize(this%h_alpha_beta, this%sub_grain_size_ptr)
  else 
    this%sub_grain_size_ptr = this%grain_diameter
  endif
  call this%computeDeltaTauLineTension(bar_alpha_ss_prime, this%tau_line_tension_ptr)
  this%tau_line_tension_ptr_old = this%tau_line_tension_ptr
  call this%computeStrength(this%crss_gauss_ptr)

  this%tau_solute_drag_gauss_ptr_old = this%tau_solute_drag_gauss_ptr
end subroutine

module subroutine computeTauM(this, solute_id, taum)
  use units_conversion_mod, only : Pa2MPa
  implicit none
  class(wkkct_glide), intent(in) :: this
  real(k_real), intent(out) :: taum
  integer, intent(in) :: solute_id

  associate(alpha_cd => this%alpha_cd_ptr(solute_id,this%ss_idx), &
            burgN => this%burgN(this%ss_idx), &
            w => this%disl_core_width_ptr(solute_id,this%ss_idx))

  taum = alpha_cd * this%computeDEcore(solute_id, this%ss_idx)/(w*burgN)*Pa2MPa
  end associate
end subroutine

module subroutine computeTauEffective(this, rss, tau_effective, dtau_eff_drss)
  use print_utils_mod, only : printToScreen
  implicit none
  class(wkkct_glide), intent(in) :: this
  real(k_real), intent(in) :: rss
  real(k_real), intent(out) :: tau_effective, dtau_eff_drss
  real(k_real) ::  chi

  chi = this%tau_line_tension_ptr(this%ss_idx )
  if (this%solutes_as_backstress) chi = chi + this%tau_m_ptr(this%ss_idx )
  if ((abs(rss) - chi).gt.0._k_real) then
    tau_effective =abs(rss) - chi
    dtau_eff_drss = sign(1._k_real, rss)
  else
    tau_effective = 0._k_real
    dtau_eff_drss = 0._k_real
  endif

  ! !AAK - quadratic blend
  ! tau_back = chi + 1e-16_k_real
  ! tau_ratio = abs(rss)/tau_back
  ! if(tau_ratio.LT.(1._k_real-back_stress_dist)) then
  !   tau_effective = 0._k_real
  !   dtau_eff_drss = 0._k_real
  ! else if(tau_ratio.LT.(1._k_real+back_stress_dist)) then
  !   a0=(1._k_real-back_stress_dist)**2/(4._k_real*back_stress_dist)*tau_back
  !   b0=(back_stress_dist-1._k_real)/(2._k_real*back_stress_dist)
  !   c0=1._k_real/(4._k_real*back_stress_dist*tau_back)
  !   tau_effective = a0+b0*abs(rss)+c0*abs(rss)**2
  !   dtau_eff_drss = b0*sign(1._k_real, rss)+2._k_real*c0*rss
  ! else
  !   tau_effective =abs(rss) - chi
  !   dtau_eff_drss = sign(1._k_real, rss)
  ! end if

end subroutine

module subroutine computeWaitTime(this, rss, crss, frequency, dg0, t_wait, dtwait_drss, use_jog_screw)
  use print_utils_mod, only : printToScreen
  use math_constants, only : kBSI
  use smooth_function_mod, only : symmSmoothClamp, dSymmSmoothClampDY
  class(wkkct_glide), intent(in) :: this
  real(k_real), intent(in) :: rss, crss, frequency, dg0
  real(k_real), intent(out) :: t_wait, dtwait_drss
  real(k_real) :: twait_jog
  logical, intent(in), optional :: use_jog_screw
  real(k_real) :: tau_effective, dtau_eff_drss, Dv
  ! real(k_real) :: normalized_tau_1, normalized_tau_2, normalized_tau, exp_argument, q_argument
  real(k_real) :: normalized_tau, exp_argument, q_argument
  ! real(k_real) :: dnormtau_dtaueff_1, dnormtau_dtaueff_2, dnormtau_dtaueff, dexparg_dqarg, dqarg_dnorm_tau
  real(k_real) :: dnormtau_dtaueff, dexparg_dqarg, dqarg_dnorm_tau
  real(k_real), parameter :: zv_disl = 1._k_real, h=100._k_real

  associate (T => this%temperature, &
              p => this%p_exp(this%ss_idx), &
              q => this%q_exp(this%ss_idx))
  
  call this%computeTauEffective(rss, tau_effective, dtau_eff_drss)
  normalized_tau = tau_effective/crss
  dnormtau_dtaueff = 1._k_real/crss

  if (normalized_tau.le.0._k_real) then
    ! q_argument = 1._k_real
    dqarg_dnorm_tau = 0._k_real
    exp_argument = DG0/(kBSI*T)
    dexparg_dqarg = 0._k_real
  else if (normalized_tau.ge.1._k_real) then
    ! q_argument = 0._k_real
    dqarg_dnorm_tau = 0._k_real
    exp_argument = 0._k_real
    dexparg_dqarg = 0._k_real
  else 
    q_argument = 1._k_real - normalized_tau**p
    dqarg_dnorm_tau = -p*normalized_tau**(p-1._k_real)
    exp_argument = DG0/(kBSI*T)*q_argument**q 
    dexparg_dqarg = DG0/(kBSI*T)*q*q_argument**(q-1._k_real)
  endif

  
  t_wait = exp(exp_argument)/frequency
  dtwait_drss = t_wait * dexparg_dqarg * dqarg_dnorm_tau * dnormtau_dtaueff * dtau_eff_drss
  if (present(use_jog_screw)) then
    if(use_jog_screw) then
      call this%computeSelfDiffusivity(Dv)
      twait_jog = h *this%burgN(this%ss_idx)**3/(this%common_material_parameter_ptr%atomic_volume*zv_disl*this%rho0_m(this%ss_idx)*Dv)
      t_wait = t_wait + twait_jog
    endif
  endif

  end associate
end subroutine

module subroutine computeTotalWaitTime(this, rss, twait_total, dtwait_tot_drss)
  class(wkkct_glide), intent(inout) :: this
  real(k_real), intent(in) :: rss
  real(k_real), intent(out) :: twait_total, dtwait_tot_drss
  real(k_real) :: t_wait_dislocation, t_wait_precipitate, &
                  dt_wait_dislocation_drss, dt_wait_precipitate_drss
  integer :: prec_id
  real(k_real), parameter :: twait_limit=sqrt(huge(1._k_real))
  associate( prob_d => this%prob_d_ptr(this%ss_idx), &
             prob_p => this%prob_p_ptr(:,this%ss_idx), &
             dg0_d => this%deltaG0_disl(this%ss_idx), &
             dg0_p => this%deltaG0_prec_ptr(:,this%ss_idx), &
             crss => this%crss_gauss_ptr(this%igauss, this%ss_idx) )

  call this%computeWaitTime(rss, crss, this%computeFrequencyDislocation(),  dg0_d, t_wait_dislocation, dt_wait_dislocation_drss, this%use_jog_screw)
  twait_total = prob_d*t_wait_dislocation
  dtwait_tot_drss = prob_d*dt_wait_dislocation_drss 

  if (this%n_prec_type>0) then
  do prec_id = 1, this%n_prec_type
    call this%computeWaitTime(rss, crss, this%computeFrequencyPrecipitate(prec_id),  dg0_p(prec_id), t_wait_precipitate, dt_wait_precipitate_drss)
    twait_total = twait_total + prob_p(prec_id)*t_wait_precipitate
    dtwait_tot_drss = dtwait_tot_drss + prob_p(prec_id)*dt_wait_precipitate_drss
  end do
  endif

  if (twait_total.ge.twait_limit) then
    twait_total = twait_limit
    dtwait_tot_drss = 0._k_real
  endif
  end associate

end subroutine

module subroutine computeGammaDotAndDGammaDotDRssSSWKKCTGlide(this, rss, gdot, dgdot_drss)
    use test_utils_mod
    use averages_mod, only : unweightedHarmonicMean
    implicit none
    class(wkkct_glide), intent(inout) :: this
    real(k_real), intent(in) :: rss
    real(k_real), intent(out) :: gdot
    real(k_real), intent(out) :: dgdot_drss

    real(k_real) :: twait_total, dtwait_tot_drss, &
                    vs, dvs_drss, &
                    travel_time, d_travel_time_drss, &
                    shear_velocity
    real(k_real), parameter :: travel_time_slope = 1e-4_k_real
    associate( lambda_s => this%lambda_s_ptr(this%ss_idx), &
               chi_e => this%chi_e, &
               rho_m => this%rho_m_ptr(this%ss_idx), &
               burgN => this%burgN(this%ss_idx), &
               time_prec_climb => this%time_prec_climb_ptr(this%igauss, this%ss_idx), &
               time_solute_aging => this%time_solute_aging_ptr(this%igauss, this%ss_idx), &
               time_solute_aging_old => this%time_solute_aging_ptr_old(this%igauss, this%ss_idx))

    !compute TotalWaitTime
    shear_velocity = this%computeShearVelocity()
    call this%computeTotalWaitTime(rss, twait_total, dtwait_tot_drss)

    time_solute_aging = max(0._k_real, (1._k_real - time_solute_aging_old/(twait_total+1e-16)) *this%dt + time_solute_aging_old)
    !  unweightedHarmonicMean((/this%time+this%pre_aging_time, twait_total/))
    time_prec_climb =  time_solute_aging

    travel_time = lambda_s/shear_velocity*(1._k_real-abs(rss)*travel_time_slope) !m/(m/s) -> s
    d_travel_time_drss = -lambda_s/shear_velocity*sign(1._k_real,rss)*travel_time_slope

    !compute the dislocation velocity
    vs = lambda_s/(twait_total + travel_time)
    dvs_drss = -lambda_s/(twait_total+travel_time)**2 * (dtwait_tot_drss + d_travel_time_drss)

    ! finaly compute gamma dot
    gdot = burgN * rho_m * vs*sign(1._k_real,rss)
    dgdot_drss =  burgN * rho_m *sign(1._k_real,rss) * dvs_drss

#ifdef __DEBUG__
    call test_is_finite(time_solute_aging, "time_solute_aging")
    call test_is_finite(time_prec_climb, "time_prec_climb")
    call test_is_finite(travel_time, "travel_time")
    call test_is_finite(gdot, "gdot")
    call test_is_finite(dgdot_drss, "dgdot_drss")
    call test_is_finite(vs, "vs")
    call test_is_finite(dvs_drss, "dvs_drss")
#endif
    end associate

  end subroutine

module subroutine updateStateVariablesAtMaterialPointInnerLoopWKKCTGlide(this)
  implicit none
  class(wkkct_glide), intent(inout) :: this
  __DECL_CLASS_UNUSED_THIS__
  ! we update everything explicitly after convergence, so nothing to do here
  __SUPPRESS_CLASS_UNUSED_THIS__
end subroutine

module subroutine updateStateVariablesAtMaterialPointOuterLoopWKKCTGlide(this)
  implicit none
  class(wkkct_glide), intent(inout) :: this

  __DECL_CLASS_UNUSED_THIS__
  ! we update everything explicitly after convergence, so nothing to do here
  __SUPPRESS_CLASS_UNUSED_THIS__

end subroutine

module subroutine updateStateVariablesAtMaterialPointStaggeredWKKCTGlide(this)
  implicit none
  class(wkkct_glide), intent(inout) :: this
  real(k_real), dimension(this%n_ss, this%n_ss) :: bar_alpha_ss_prime
  ! TODO WE NEED TO CHECK TOGHETHER THE ORDER IN WHICH VARIABLES ARE UPDATED
  ! I FOLLOWED THE ORDER NATHAN USED, ALSO I

  call this%computeResolvedStressAndReslovedStressDStress(this%stress6)
  associate (ss_idx => this%ss_idx)
  do ss_idx=1,this%n_ss
    call this%getGaussianValueAndIntegrationWeights(this%rss(ss_idx), this%std_dev_ptr, this%rss_gaussian, this%total_integration_weight_ptr(:,this%ss_idx))
    if (this%n_solute_type > 0) this%time_solute_aging_avg_ptr(ss_idx) = sum(this%time_solute_aging_ptr(:,ss_idx)*this%total_integration_weight_ptr(:,this%ss_idx))
  enddo
  end associate

  call this%computeDislocationDensity(this%rho_m_ptr, this%rho_cw_ptr)
  if (this%n_prec_type>0) then
  call this%computeEffectivePrecipitateDiameter(this%diameter_prec_ptr)
  endif
  call this%computeEffectiveHardeningMatrix(bar_alpha_ss_prime)
  call this%computeEffectiveDislocationInterspacingAndProbability(bar_alpha_ss_prime, this%h_alpha_beta, this%lambda_s_ptr, this%prob_d_ptr, this%prob_p_ptr)
  if (this%use_dislocation_cell_wall) then
    call this%computeSubGrainSize(this%h_alpha_beta, this%sub_grain_size_ptr)
  endif
  call this%computeDeltaTauLineTension(bar_alpha_ss_prime, this%tau_line_tension_ptr)

  if (this%evolve_std_dev) call this%computeStdDeviation(this%std_dev_ptr)
  
  call this%computeStrength(this%crss_gauss_ptr)
#ifdef __DEBUG__
  call test_is_finite(this%time_solute_aging_ptr, "time_solute_aging_ptr")
  call test_is_finite(this%time_prec_climb_ptr, "time_prec_climb_ptr")
  call test_is_finite(this%rho_m_ptr, "rho_m_ptr")
  if (this%use_dislocation_cell_wall) then
    call test_is_finite(this%rho_cw_ptr, "rho_cw_ptr")
  endif
  call test_is_finite(this%diameter_prec_ptr, "diameter_prec_ptr")
  call test_is_finite(this%lambda_s_ptr, "lambda_s_ptr")
  call test_is_finite(this%tau_line_tension_ptr, "tau_line_tension_ptr")
  call test_is_finite(this%tau_solute_drag_gauss_ptr, "tau_solute_drag_gauss_ptr")
  call test_is_finite(this%prob_d_ptr, "prob_d_ptr")
  call test_is_finite(this%prob_p_ptr, "prob_p_ptr")
  call test_is_finite(this%std_dev_ptr, "std_dev_ptr")
  call test_is_finite(this%crss_gauss_ptr, "crss_gauss_ptr")
#endif
end subroutine


module subroutine computeEffectivePrecipitateDiameter(this, eff_d_prec)
  implicit none
  class(wkkct_glide), intent(inout) :: this
  real(k_real), dimension(this%n_prec_type, this%n_ss), intent(out) :: eff_d_prec
  real(k_real) :: vclimb_eff_local,twait_eff_local
  real(k_real), parameter :: min_prec_size = 1e-20_k_real
  integer :: prec_id

  associate(initial_diameter_prec => this%initial_diameter_prec_ptr, &
            pdf_local => this%total_integration_weight_ptr, &
            ss_idx => this%ss_idx, &
            igauss => this%igauss, &
            time_prec_climb => this%time_prec_climb_ptr, &
            v_climb => this%climb_velocity_ptr)


  do ss_idx = 1, this%n_ss
    ! compute effective climb velocity and time to by-pass.
    vclimb_eff_local = sum(v_climb(:,ss_idx) * pdf_local(:,ss_idx))
    twait_eff_local = sum(time_prec_climb(:,ss_idx) * pdf_local(:,ss_idx))
    do prec_id = 1, this%n_prec_type
        ! compute effective precipitate size for each type and for each SS.
        eff_d_prec(prec_id,ss_idx) = initial_diameter_prec(prec_id)*exp(-2._k_real * twait_eff_local* abs(vclimb_eff_local) / initial_diameter_prec(prec_id))

    end do
  end do
  end associate
end subroutine

! SUBROUTINE USED FOR UPDATING STATE VARIABLES

module subroutine computeDislocationDensity(this, rho_m, rho_cw)
  use math_constants, only : kBSI
  use units_conversion_mod, only : MPa2Pa
  implicit none
  class(wkkct_glide), intent(in) :: this
  real(k_real), intent(inout), dimension(this%n_ss) :: rho_m, rho_cw
  real(k_real), dimension(this%n_ss)  :: rho_dot

  ! temporary variables
  real(k_real) :: abs_gamma_dot(this%n_ss), &
                  rho_dot_cw(this%n_ss)

  associate(rho_dot_m => this%rho_m_dot_ptr)
  abs_gamma_dot = abs(this%gammadot_ptr)

  rho_dot = this%computeRhoDotGeneration(abs_gamma_dot)
  rho_dot_m = rho_dot 
  this%delta_rho_m_generation_ptr = this%delta_rho_m_generation_ptr + rho_dot*this%dt

  rho_dot = this%computeRhoDotDynamicRecovery(abs_gamma_dot)
  rho_dot_m =  rho_dot_m - rho_dot
  this%delta_rho_m_dyn_recovery_ptr = this%delta_rho_m_dyn_recovery_ptr+rho_dot *this%dt

  if (this%use_dislocation_cell_wall) then
    rho_dot = this%computeRhoDotTrappedCellWall(abs_gamma_dot)
    rho_dot_m =  rho_dot_m - rho_dot
    this%delta_rho_m_trapped_cw_ptr = this%delta_rho_m_trapped_cw_ptr+rho_dot *this%dt
  endif

  if (this%use_static_recovery) then
    rho_dot= this%computeRhoDotStaticRecovery()
    rho_dot_m =  rho_dot_m - rho_dot
    this%delta_rho_m_static_recovery_ptr = this%delta_rho_m_static_recovery_ptr+rho_dot *this%dt
  endif               

  if (this%use_dislocation_cell_wall) then
    rho_dot = this%computeRhoDotTrappedCellWall(abs_gamma_dot)
    rho_dot_cw = rho_dot

    rho_dot = this%computeRhoDotAnnihilationCellWallClimbVelocity()
    rho_dot_cw = rho_dot_cw - rho_dot
    this%delta_rho_cw_annihilation_ptr = this%delta_rho_cw_annihilation_ptr + rho_dot *this%dt

    rho_cw = rho_cw + rho_dot_cw * this%dt
    rho_cw = max(rho_cw, 10._k_real)  
  else 
    rho_dot_cw = 0._k_real
    rho_cw = 0._k_real
  endif
    
  rho_m = rho_m + rho_dot_m * this%dt
  ! min dislocation density value
  rho_m = max(rho_m, 10._k_real)
  
  end associate
end subroutine

module function computeRhoDotGeneration(this, abs_gamma_dot) result(rho_dot)
  implicit none
  class(wkkct_glide), intent(in) :: this
  real(k_real), intent(in) :: abs_gamma_dot(this%n_ss)
  real(k_real) :: rho_dot(this%n_ss)

  rho_dot = this%k_disl_gen/this%lambda_s_ptr * abs_gamma_dot
end function

module function computeK2DynamicRecovery(this, gdot) result(k2)
  use math_constants, only :kBSI
  use units_conversion_mod, only :MPa2Pa
  implicit none
  class(wkkct_glide), intent(in) :: this
  real(k_real), intent(in) :: gdot(this%n_ss)
  real(k_real) :: k2(this%n_ss)

  k2 = (this%chi_dynrec*this%burgN/this%g_dynrec)* &
        (1._k_real - ((kBSI*this%temperature)/(this%D_dynrec*MPa2Pa*this%burgN**3))*log((abs(gdot)+1e-25)/this%edot0_dynrec))* &
    this%k_disl_gen
end function

module function computeRhoDotDynamicRecovery(this, abs_gamma_dot) result(rho_dot)
  implicit none
  class(wkkct_glide), intent(in) :: this
  real(k_real), intent(in) :: abs_gamma_dot(this%n_ss)
  real(k_real) :: rho_dot(this%n_ss)

  rho_dot = this%rho_m_ptr*this%computeK2DynamicRecovery(this%gammadot_ptr)*abs_gamma_dot
end function

module function computeRhoDotStaticRecovery(this) result(rho_dot)
  use math_constants, only :kBSI
  use units_conversion_mod, only :MPa2Pa
  implicit none
  class(wkkct_glide), intent(in) :: this
  real(k_real) :: rho_dot(this%n_ss), Dv, Dp, Deff(this%n_ss)

  if (this%use_static_recovery) then
    call this%computeSelfDiffusivity(Dv)
    call this%computePipeDiffusivity(Dp)

    associate (rho => this%rho_m_ptr(1:this%n_ss), &
               rho_tot => sum(this%rho_m_ptr(1:this%n_ss)), &
               bavg => sum(this%burgN(1:this%n_ss))/this%n_ss, &
               k1 => this%k1_static_recov(this%gb_id), &
               k2 => this%k2_static_recov(this%gb_id), &
               mu => this%shear_mod*MPa2Pa, &
               T => this%temperature, &
               omega => this%common_material_parameter_ptr%atomic_volume, &
               sum_rho_b2 =>sum(this%rho_m_ptr(1:this%n_ss) * this%burgN(1:this%n_ss)**2))

    Deff = Dv + Dp*sum_rho_b2
    
    ! here we use Eq-9 not Eq-10, from Kohnert and Capolungo,
    ! because the latter does not work if rho>rho0
    rho_dot = rho/rho_tot * & ! rescale by slip system
               (k1*Deff/bavg) * rho_tot**(1.5_k_real) * &
               ( exp( (k2*mu*omega*sqrt(sum_rho_b2))/(kBSI*T) ) -1._k_real)

    end associate 
  else 
    rho_dot = 0._k_real
  endif


end function

module function computeRhoDotTrappedCellWall(this, abs_gamma_dot) result(rho_dot)
  implicit none
  class(wkkct_glide), intent(in) :: this
  real(k_real), intent(in) :: abs_gamma_dot(this%n_ss)
  real(k_real) :: rho_dot(this%n_ss)
  if (this%use_dislocation_cell_wall) then
    rho_dot = this%k_disl_trap*abs_gamma_dot/this%sub_grain_size_ptr
  else 
    rho_dot = 0._k_real
  endif
end function

module function computeRhoDotAnnihilationCellWall(this, abs_gamma_dot) result(rho_dot)
  implicit none
  class(wkkct_glide), intent(in) :: this
  real(k_real), intent(in) :: abs_gamma_dot(this%n_ss)
  real(k_real) :: rho_dot(this%n_ss)
  if (this%use_dislocation_cell_wall) then
    rho_dot = this%k_cell_wall_ann * this%rho_cw_ptr * abs_gamma_dot
  else 
    rho_dot = 0._k_real
  endif
end function

module function computeRhoDotAnnihilationCellWallClimbVelocity(this) result(rho_dot)
  implicit none
  class(wkkct_glide), intent(in) :: this
  real(k_real), dimension(this%n_ss) :: rho_dot, vclimb_eff_local
  integer :: ss_idx

  if (this%use_dislocation_cell_wall) then
    associate(v_climb_CW => this%climb_velocity_CellWall_ptr, &
              pdf_local => this%total_integration_weight_ptr)
    vclimb_eff_local = 0._k_Real
    do ss_idx = 1, this%n_ss
      ! compute effective climb velocity in the cell wall.
      vclimb_eff_local(ss_idx) = sum(v_climb_CW(:,ss_idx) * pdf_local(:,ss_idx))
    end do

    rho_dot = this%k_cell_wall_ann * abs(vclimb_eff_local) * ((this%rho_cw_ptr)**(1.5_k_real))

    end associate
  else 
    rho_dot = 0._k_real
  endif
end function

module subroutine computeStrength(this, crss)
  use print_utils_mod, only : printToScreen
  implicit none
  integer :: prec_id, solute_id, loop_id
  class(wkkct_glide), intent(inout) :: this
  real(k_real), dimension(this%n_gauss, this%n_ss), intent(inout) :: crss
  real(k_real), dimension(this%n_ss) :: tau_Labsuch, dtau_disloc_loop, dtau_prec
  real(k_real), dimension(this%n_ss) :: tau0
  call this%tau0_linear_interp%intepolateAllSlipSystems(this%temperature, tau0)
  ! Keep the following lines for a while, so that we can easily go back to linear interpoltation instead of Labusch equation.
  ! call this%tau_labusch_linear_interp%intepolateAllSlipSystems(this%temperature, taulabusch)

  ! hall-petch
  
  associate(burg => this%burgN, &
            n_crss => this%n_crss,  &
            igauss => this%igauss,&
            ss_idx => this%ss_idx)
  
  
  crss = 0._k_real
  this%tau_hall_petch_ptr = this%k_hall_petch*this%shear_mod*this%burgN/sqrt(this%grain_diameter) 
  do igauss=1,this%n_gauss
  crss(igauss, :) = crss(igauss, :) + this%tau_hall_petch_ptr  ** n_crss
  enddo 

  if (this%use_dislocation_cell_wall) then 
    this%tau_cell_wall_ptr = this%shear_mod*this%burgN*sqrt(matmul(this%h_alpha_beta, this%rho_cw_ptr))
    do igauss=1,this%n_gauss
      crss(igauss, :) = crss(igauss, :) + this%tau_cell_wall_ptr  ** n_crss
    enddo 
  endif

  if(this%n_prec_type>0) then
  this%tau_precipitate_drag_ptr  = 0._k_real
  do prec_id = 1, this%n_prec_type
    call this%computeDeltaTauPrecipitateBKS(prec_id, dtau_prec)
    this%tau_precipitate_drag_ptr = this%tau_precipitate_drag_ptr + dtau_prec**n_crss
  end do
  this%tau_precipitate_drag_ptr = (this%tau_precipitate_drag_ptr)**(1._k_real/n_crss)
  do igauss=1,this%n_gauss
    crss(igauss, :) = crss(igauss, :) + this%tau_precipitate_drag_ptr  ** n_crss
  enddo 
  endif
  
  
  if(this%n_solute_type>0) then
    this%tau_solute_drag_gauss_ptr = 0._k_real
    call this%computeDeltaTauSoluteDrag(this%tau_solute_drag_gauss_ptr)
    this%tau_labusch_ptr = 0._k_real
    do solute_id = 1, this%n_solute_type
        call this%computeLabuschStrength(solute_id,tau_Labsuch)
        this%tau_labusch_ptr = this%tau_labusch_ptr + tau_Labsuch**n_crss
    end do
    this%tau_labusch_ptr = (this%tau_labusch_ptr)**(1._k_real/n_crss)

    if (this%solutes_as_backstress) then
      do igauss=1,this%n_gauss
        crss(igauss, :) = crss(igauss, :) + (this%tau_labusch_ptr)  ** n_crss
      enddo 
    else 
      do igauss=1,this%n_gauss
        crss(igauss, :) = crss(igauss, :) + (this%tau_labusch_ptr + this%tau_m_ptr)  ** n_crss
      enddo 
    endif
  endif

  if(this%n_dislocation_loop_type.gt.0) then
      this%tau_disl_loop_ptr = 0._k_real
      do loop_id = 1, this%n_dislocation_loop_type
        call this%computeDeltaTauDislocationLoopsBKS(loop_id, dtau_disloc_loop)
        this%tau_disl_loop_ptr = this%tau_disl_loop_ptr + dtau_disloc_loop**n_crss
      end do
      this%tau_disl_loop_ptr = (this%tau_disl_loop_ptr)**(1._k_real/n_crss)
      do igauss=1,this%n_gauss
        crss(igauss, :) = crss(igauss, :) + this%tau_disl_loop_ptr  ** n_crss
      enddo 
  endif

  do igauss=1,this%n_gauss
    crss(igauss, :) = crss(igauss, :) + tau0**n_crss
  enddo 
  do igauss=1,this%n_gauss
    crss(igauss, :) = crss(igauss, :)**(1._k_real/n_crss)
  enddo 
  
  do ss_idx=1,this%n_ss
    this%crss_ptr(ss_idx) = sum(crss(:, ss_idx)*this%total_integration_weight_ptr(:,ss_idx))
  enddo 
  end associate
end subroutine

module subroutine computeDeltaTauPrecipitate(this, prec_id, dtau_prec)
  implicit none
  integer, intent(in) :: prec_id
  class(wkkct_glide), intent(in) :: this
  real(k_real), dimension(this%n_ss), intent(out) :: dtau_prec

  associate(burg => this%burgN, &
            alpha_prec => this%alpha_prec_ptr(prec_id,1:this%n_ss), &
            diameter_prec => this%diameter_prec_ptr, &
            number_density_prec => this%number_density_prec_ptr(prec_id))

  dtau_prec = this%shear_mod*burg*alpha_prec*sqrt(number_density_prec*diameter_prec(prec_id,:))

  end associate
end subroutine

module subroutine computeDeltaTauPrecipitateBKS(this, prec_id, dtau_prec)
  use math_constants, only : PI
  implicit none
  integer, intent(in) :: prec_id
  class(wkkct_glide), intent(in) :: this
  real(k_real), dimension(this%n_ss), intent(out) :: dtau_prec
  real(k_real), dimension(this%n_ss) :: center_center_spacing, inner_spacing, harmonic_mean_size
  integer :: ss_idx


  associate(burg => this%burgN, &
            poisson => this%poisson, &
            alpha_prec => this%alpha_prec_ptr(prec_id,1:this%n_ss), &
            diameter_prec => this%diameter_prec_ptr, &
            number_density_prec => this%number_density_prec_ptr(prec_id))

  if (this%n_prec_type>0) then
  center_center_spacing = 1._k_real/sqrt(number_density_prec*diameter_prec(prec_id,:))
  inner_spacing = center_center_spacing - diameter_prec(prec_id,:)
  harmonic_mean_size = (1._k_real/diameter_prec(prec_id,:)) + (1._k_real/inner_spacing)
  harmonic_mean_size = 1._k_real/harmonic_mean_size
  do ss_idx = 1, this%n_ss
    if (2._k_real*harmonic_mean_size(ss_idx).le.burg(ss_idx)) harmonic_mean_size(ss_idx) = burg(ss_idx)
    if (inner_spacing(ss_idx).le.burg(ss_idx)) inner_spacing(ss_idx) = burg(ss_idx)
  end do
  dtau_prec = (1._k_real/(2._k_real*PI*sqrt(1._k_real-poisson)))*this%shear_mod*burg/inner_spacing
  dtau_prec = dtau_prec * (log(2._k_real*harmonic_mean_size/burg))**(1.5_k_real)
  dtau_prec = dtau_prec * (log(inner_spacing/burg))**(-0.5_k_real)
  else
    dtau_prec=0._k_real
  endif

  end associate
end subroutine

module subroutine computeDeltaTauDislocationLoopsBKS(this, loop_id, dtau_disloc_loop)
  use math_constants, only : PI
  implicit none
  integer, intent(in) :: loop_id
  class(wkkct_glide), intent(in) :: this
  real(k_real), dimension(this%n_ss), intent(out) :: dtau_disloc_loop
  real(k_real), dimension(this%n_ss) :: center_center_spacing, inner_spacing, harmonic_mean_size
  integer :: ss_idx

  associate(burg => this%burgN, &
            poisson => this%poisson, &
            diameter_disloc_loops => this%initial_diameter_disloc_loops_ptr(loop_id), &
            number_density_disloc_loops => this%number_density_disloc_loops_ptr(loop_id), &
            loop_crystallography_factor_temp => this%loop_crystallography_factor_ptr(:,loop_id) )

  center_center_spacing = 1._k_real/sqrt(number_density_disloc_loops*diameter_disloc_loops)
  inner_spacing = center_center_spacing - diameter_disloc_loops
  harmonic_mean_size = (1._k_real/diameter_disloc_loops) + (1._k_real/inner_spacing)
  harmonic_mean_size = 1._k_real/harmonic_mean_size

  do ss_idx = 1, this%n_ss
    if (2._k_real*harmonic_mean_size(ss_idx).le.burg(ss_idx)) harmonic_mean_size(ss_idx) = burg(ss_idx)
    if (inner_spacing(ss_idx).le.burg(ss_idx)) inner_spacing(ss_idx) = burg(ss_idx)
  end do
  dtau_disloc_loop = (1._k_real/(2._k_real*PI*sqrt(1._k_real-poisson)))*this%shear_mod*burg/inner_spacing
  dtau_disloc_loop = dtau_disloc_loop * (log(2._k_real*harmonic_mean_size/burg))**(1.5_k_real)
  dtau_disloc_loop = dtau_disloc_loop * (log(inner_spacing/burg))**(-0.5_k_real)
  dtau_disloc_loop = dtau_disloc_loop * loop_crystallography_factor_temp
  end associate
end subroutine

module subroutine computeLabuschStrength(this, solute_id, tau_Labsuch)
  implicit none
  integer, intent(in) :: solute_id
  class(wkkct_glide), intent(in) :: this
  real(k_real), dimension(this%n_ss), intent(out) :: tau_Labsuch

  associate(solute_concentration => this%solute_concentration_ptr(solute_id), &
            K_Labusch => this%K_Labusch_ptr(solute_id,:), &
            exponent_Labusch => this%exponent_Labusch_ptr(solute_id,:))
  tau_Labsuch = 0._k_real
  if(this%n_solute_type>0) then
    tau_Labsuch = K_Labusch*this%shear_mod*(solute_concentration**exponent_Labusch)
  endif
  end associate
end subroutine

module subroutine computeDeltaTauLineTension(this, alpha_ss_mobile, tau_line_tension)
  implicit none
  class(wkkct_glide), intent(in) :: this
  real(k_real), dimension(this%n_ss, this%n_ss), intent(in) :: alpha_ss_mobile
  real(k_real), dimension(this%n_ss), intent(out) :: tau_line_tension

  associate(burg => this%burgN)

  tau_line_tension = this%shear_mod*burg*sqrt(matmul(alpha_ss_mobile, this%rho_m_ptr))

#ifdef __DEBUG__
  call test_is_finite(tau_line_tension, "tau_line_tension")
#endif
  end associate
end subroutine

module subroutine computeDeltaTauSoluteDrag(this, tau_sol_drag)
  use mpi_variables_mod, only : mpi_master_rank
  implicit none
  class(wkkct_glide), intent(inout) :: this
  real(k_real), dimension(this%n_gauss, this%n_ss), intent(out) :: tau_sol_drag
  real(k_real) :: taum_local_solute, taum_effective
  integer :: solute_id
  !integer :: igauss, ss_idx

  associate(burg => this%burgN,&
            igauss => this%igauss,&
            ss_idx => this%ss_idx)

  tau_sol_drag = 0._k_real
  if (this%n_solute_type>0) then
  do ss_idx=1,this%n_ss
    do igauss=1,this%n_gauss
    taum_effective = 0._k_real
    do solute_id = 1, this%n_solute_type
      call this%computeTauM(solute_id, taum_local_solute)
        ! tau_sol_drag(ss_idx) = tau_sol_drag(ss_idx) + taum_local * this%total_integration_weight_ptr(igauss, ss_idx)
        taum_effective = taum_effective + taum_local_solute**this%n_crss(ss_idx)
      enddo
      ! write(*,*) 'ss_idx,taum_effective', ss_idx,taum_effective
      tau_sol_drag(igauss, ss_idx) = taum_effective**(1._k_real/this%n_crss(ss_idx))
    end do
  enddo
  endif
  do ss_idx=1,this%n_ss
    this%tau_m_ptr(ss_idx) = sum(tau_sol_drag(:, ss_idx)*this%total_integration_weight_ptr(:,ss_idx))
  enddo 

#ifdef __DEBUG__
  call test_is_finite(tau_sol_drag, "tau_sol_drag")
#endif

  end associate
end subroutine

module subroutine computeEffectiveDislocationInterspacingAndProbability(this, alpha_ss_mobile, alpha_ss_cw, lambda_s, prob_d,prob_p)
  implicit none
  integer :: prec_id
  class(wkkct_glide), intent(in) :: this
  real(k_real), dimension(this%n_ss, this%n_ss), intent(in) :: alpha_ss_mobile, alpha_ss_cw
  real(k_real), dimension(this%n_ss), intent(out) :: lambda_s, prob_d
  real(k_real), dimension(this%n_prec_type, this%n_ss), intent(out) :: prob_p
  real(k_real), dimension(this%n_ss) :: one_over_lambda_disl, one_over_lambda, one_over_lambda_hall_petch
  real(k_real), dimension(this%n_prec_type,this%n_ss) ::one_over_lambda_prec
  

  associate(diameter_prec => this%initial_diameter_prec_ptr)

  one_over_lambda = 0.
  one_over_lambda_prec = 0.
  one_over_lambda_disl = sqrt(matmul(alpha_ss_mobile, this%rho_m_ptr)) + sqrt(matmul(alpha_ss_cw, this%rho_cw_ptr)) !1/m
  this%lambda_dislocation_ptr = 1._k_real/one_over_lambda_disl
  one_over_lambda_hall_petch = this%k_hall_petch/sqrt(this%grain_diameter)
  if (this%k_hall_petch.gt.0._k_real) then 
    this%lambda_hall_petch_ptr = 1._k_real/one_over_lambda_hall_petch
  else 
    this%lambda_hall_petch_ptr = -1._k_real
  endif
  one_over_lambda = one_over_lambda_disl + one_over_lambda_hall_petch

  if (this%n_prec_type>0) then
  do prec_id = 1, this%n_prec_type
    one_over_lambda_prec(prec_id,:) = this%beta_prec_ptr(prec_id,:)*sqrt(this%number_density_prec_ptr(prec_id)*diameter_prec(prec_id)) !1/m
    one_over_lambda = one_over_lambda + one_over_lambda_prec(prec_id,:)
  end do
    if (all(this%beta_prec_ptr>0._k_real)) then
      this%lambda_precipitate_ptr = 1._k_real/sum(one_over_lambda_prec, dim=1)
    else 
      this%lambda_precipitate_ptr = -1._k_real
    endif
  endif

  lambda_s = 1._k_real/(one_over_lambda) !m
  prob_d = (one_over_lambda_disl)*lambda_s
  prob_p = 0._k_real
  if (this%n_prec_type>0) then
  do prec_id = 1, this%n_prec_type
    prob_p(prec_id,:) = (one_over_lambda_prec(prec_id,:))*lambda_s
  end do
  endif
  end associate

end subroutine

module subroutine computeSubGrainSize(this, alpha_ss_cw, sub_grain_size)
  implicit none
  class(wkkct_glide), intent(in) :: this
  real(k_real), dimension(this%n_ss, this%n_ss), intent(in) :: alpha_ss_cw
  real(k_real), dimension(this%n_ss), intent(out) :: sub_grain_size
  real(k_real), dimension(this%n_ss) :: one_over_lambda_cw
  if (this%use_dislocation_cell_wall) then
    associate(grain_diameter => this%common_material_parameter_ptr%grain_diameter, &
              burgN => this%burgN)

    one_over_lambda_cw = sqrt(matmul(alpha_ss_cw, this%rho_cw_ptr))
      sub_grain_size = 1._k_real/(one_over_lambda_cw)
    where (sub_grain_size.ge.grain_diameter) sub_grain_size = grain_diameter
    where (sub_grain_size.le.burgN) sub_grain_size = burgN
    end associate
  else 
    sub_grain_size = this%grain_diameter
  endif
end subroutine

module subroutine computeStdDeviation(this, std_dev)
  implicit none
  class(wkkct_glide), intent(in) :: this
  real(k_real), intent(out) :: std_dev

  associate(etav => this%etav)

  std_dev = etav*this%shear_mod*sqrt(sum(this%rho_m_ptr + this%rho_cw_ptr)) ! this is the standard deviation

  end associate

end subroutine

module function computeDEcoreSaturation(this, solute_id, ss_idx) result(DEcore_sat)
  use math_constants, only : kBSI
  implicit none
  class(wkkct_glide), intent(in) :: this
  integer, intent(in) :: ss_idx, solute_id
  real(k_real) :: DEcore_sat
  DEcore_sat = 0._k_real

  associate( solute_concentration => this%solute_concentration_ptr(solute_id), &
             DW => this%DW_ptr(solute_id,ss_idx), &
             burgN => this%burgN(ss_idx), &
             temperature => this%temperature, &
             N => this%N_atoms_disl_core_ptr(solute_id,ss_idx))
  if (this%n_solute_type>0) then
    DEcore_sat = solute_concentration * N * DW * tanh(DW/(2._k_real*kBSI*temperature)) ! this comes out in newton
  endif
#ifdef __DEBUG__
  call test_is_finite(DEcore_sat, "DEcore_sat ")
#endif

  end associate
end function

module function computeCoreDiffusionTime(this, solute_id, ss_idx) result(td)
  use math_constants, only : kBSI
  implicit none
  class(wkkct_glide), intent(in) :: this
  integer, intent(in) :: ss_idx, solute_id
  real(k_real) :: td

  associate( DHc => this%DHc_ptr(solute_id,ss_idx), &
             DW => this%DW_ptr(solute_id,ss_idx), &
             freq_cd => this%freq_cd_ptr(solute_id,ss_idx), &
             temperature => this%temperature, &
             m_cd => this%m_cd_ptr(solute_id,ss_idx))

    td = exp((DHc-DW/2._k_real)/(kBSI*temperature))/(m_cd*freq_cd)

  end associate
end function

module function computeDEcore(this, solute_id, ss_idx) result(DEcore)
  implicit none
  class(wkkct_glide), intent(in) :: this
  integer, intent(in) :: ss_idx, solute_id
  real(k_real) :: DEcore
  DEcore = 0._k_real
  if (this%n_solute_type>0) then
  associate( phi_cd => this%phi_cd_ptr(solute_id, ss_idx), &
             t_aging => this%time_solute_aging_ptr(this%igauss, ss_idx), &
             td =>  this%computeCoreDiffusionTime(solute_id, ss_idx), &
             DEcoreSat => this%computeDEcoreSaturation(solute_id, ss_idx))

    DEcore = DEcoreSat * (1._k_real - exp(-(t_aging/td)**phi_cd))

#ifdef __DEBUG__
  call test_is_finite(DEcore, "DEcore ")
#endif

  end associate
  endif
end function

module function computeShearVelocity(this) result(shear_velocity)
  use units_conversion_mod, only : MPa2Pa
  implicit none
  class(wkkct_glide), intent(in) :: this
  real(k_real) :: shear_velocity
  shear_velocity = SQRT(this%shear_mod*MPa2Pa/this%rho_mass)
  !  N    m^3    Kg   m^3   m^2
  ! --- * --- = --------- = ---
  ! m^2   Kg    s^2 m Kg    s^2
end function

module function computeBindingEnergyNormalizingFactor(this) result(A)
  use units_conversion_mod, only : MPa2Pa
  use math_constants, only : PI
  implicit none
  class(wkkct_glide), intent(in) :: this
  real(k_real), dimension(this%n_ss) :: A

  A=3._k_real/4._k_real*this%shear_mod*MPa2Pa*this%burgN**2/(PI*(1._k_real-this%poisson))
end function

module function computeFrequencyDislocation(this) result(disl_freq)
  implicit none
  class(wkkct_glide), intent(in) :: this
  real(k_real) :: disl_freq

  associate(chi_e => this%chi_e(this%ss_idx), &
            lambda_s => this%lambda_s_ptr(this%ss_idx))
  disl_freq = chi_e*this%computeShearVelocity()/lambda_s
  ! 1/s = m/s/m
  end associate

end function

module function computeFrequencyPrecipitate(this, prec_id) result(prec_freq)
  implicit none
  class(wkkct_glide), intent(in) :: this
  real(k_real) :: prec_freq
  integer, intent(in) :: prec_id
  prec_freq = this%freq_prec_ptr(prec_id, this%ss_idx)
end function

module subroutine computeEffectiveHardeningMatrix(this, alpha_ss)
  implicit none
  class(wkkct_glide), intent(inout) :: this
  real(k_real), intent(out) :: alpha_ss(this%n_ss, this%n_ss)
  real(k_real) :: A(this%n_ss)
  real(k_real) :: DEcore_effective,chi_times_Decore_by_A
  integer :: i, solute_id
  associate(a_hard => this%h_alpha_beta, &
            pdf_local => this%total_integration_weight_ptr, &
            chi_cd_ptr => this%chi_cd_ptr,&
            igauss => this%igauss, &
            ss_idx => this%ss_idx)


  A = this%computeBindingEnergyNormalizingFactor()
  do i=1, this%n_ss
    do ss_idx=1, this%n_ss
      alpha_ss(i,ss_idx) = 0._k_real
      chi_times_Decore_by_A = 0._k_real
      if (this%n_solute_type>0) then
      do solute_id = 1, this%n_solute_type
        DEcore_effective = 0._k_real
        do igauss = 1,this%n_gauss
          DEcore_effective = DEcore_effective + pdf_local(igauss,ss_idx)*this%computeDEcore(solute_id, ss_idx)
        enddo
        chi_times_Decore_by_A = chi_times_Decore_by_A + chi_cd_ptr(solute_id,ss_idx)*DEcore_effective/A(ss_idx)
      enddo
      endif
      alpha_ss(i,ss_idx) = a_hard(i,ss_idx) * (1._k_real+chi_times_Decore_by_A)**2
    end do
  enddo


  end associate
end subroutine

module subroutine acceptRejectSolutionWKKCTGlide(this, dt_max, accept_solution_flag)
  use mpi_useful_routines_mod, only : MPIComputeQualitySSScalar, MPIMaxIncrementGridSSScalar
  use mpi_variables_mod, only : i_am_mpi_master
  use time_march_mod, only : computeMaxAllowedTimeStepQuality
  use log_file_mod, only : write_detailed_log_to_screen, writeToScreen
  use number_to_string_mod
  implicit none
  class(wkkct_glide), intent(inout) :: this
  real(k_real), intent(out) :: dt_max
  logical, intent(out) :: accept_solution_flag
  logical :: accept_var
  real(k_real) :: max_dt_var, Q
  real(k_real) :: max_rho_m_increment !, max_rho_cw_rel_increment, max_tau_increment
  integer :: npoints, dims(4)

  call writeToScreen("*********************************")
  call writeToScreen("acceptRejectSolutionWKKCTGlide")

  accept_solution_flag = .true.
  dt_max = this%dt_max
  npoints = this%grid_data%getGlobalGridNPoints()
  dims = shape(this%gammadot_grid)
  npoints = npoints * dims(1)

  call MPIComputeQualitySSScalar(this%rho_m_grid,  this%rho_m_grid_old, this%rho_m_dot_grid,  this%rho_m_dot_grid_old, &
  npoints, this%dt, 1e-6_k_real, 1e8_k_real, Q, max_rho_m_increment)
  call computeMaxAllowedTimeStepQuality(Q, this%dt, "rho mobile ",  accept_var, max_rho_m_increment, max_dt_var)
  accept_solution_flag = accept_solution_flag.and.accept_var
  dt_max = min(dt_max, max_dt_var)

  ! ! check cell wall dd increment
  ! if (this%use_dislocation_cell_wall) then
  !   call MPIComputeQualitySSScalar(this%rho_m_grid,  this%rho_m_grid_old, this%rho_m_dot_grid,  this%rho_m_dot_grid_old, &
  !   npoints, this%dt, 1e-3_k_real, 1e6_k_real, Q, max_rho_m_increment)
  !   call computeMaxAllowedTimeStepQuality(Q, this%dt, "rho mobile ",  accept_var, max_rho_m_increment, max_dt_var)
  !   accept_solution_flag = accept_solution_flag.and.accept_var
  !   dt_max = min(dt_max, max_dt_var)
  ! endif

  ! ! check tau solute drag increment
  ! if (this%n_solute_type>0) then
  !     call MPIComputeQualitySSScalar(this%tau_solute_drag_gauss_grid,  this%tau_solute_drag_gauss_grid_old, &
  !                                    (this%tau_solute_drag_gauss_grid-this%tau_solute_drag_gauss_grid_old)/this%dt, &
  !                                    (this%tau_solute_drag_gauss_grid_old-this%tau_solute_drag_gauss_grid_older)/this%dt, &
  !                                    npoints, this%dt, 1e-3_k_real, 1e-3_k_real, Q, max_tau_increment)
  !     call computeMaxAllowedTimeStepQuality(Q, this%dt, "tau_solute ",  accept_var, max_tau_increment, max_dt_var)
  !     accept_solution_flag = accept_solution_flag.and.accept_var
  !     dt_max = min(dt_max, max_dt_var)
  ! endif

  ! ! check line tension solute increment
  ! call MPIComputeQualitySSScalar(this%tau_line_tension_grid,  this%tau_line_tension_grid_old, &
  !                                (this%tau_line_tension_grid-this%tau_line_tension_grid_old)/this%dt, &
  !                                (this%tau_line_tension_grid_old-this%tau_line_tension_grid_older)/this%dt, &
  !                                npoints, this%dt, 1e-3_k_real, 1e-3_k_real, Q, max_tau_increment)
  ! call computeMaxAllowedTimeStepQuality(Q, this%dt, "tau_line_tenstion ",  accept_var, max_tau_increment, max_dt_var)
  ! accept_solution_flag = accept_solution_flag.and.accept_var
  ! dt_max = min(dt_max, max_dt_var)

  if (accept_solution_flag) call writeToScreen("solution might be ACCEPTED")
  if (.not.accept_solution_flag) call writeToScreen("solution will be REJECTED")
  if (i_am_mpi_master.and.write_detailed_log_to_screen) write(*,*) "new allowable dt is ", dt_max
  call writeToScreen("*********************************")
  call writeToScreen("")

end subroutine

module subroutine writeAverageQuantitiesCSVWKKCTGlide(this, csv_writer_obj, write_headers)
  use mpi_useful_routines_mod, only : MPIAverageGridVectorMasked, MPIAverageGridScalarMasked
  use csv_writer_mod, only : csv_writer
  use number_to_string_mod, only : int2string

  implicit none
  class(wkkct_glide), intent(inout) :: this
  class(csv_writer), intent(inout) :: csv_writer_obj
  logical, intent(in) :: write_headers
  real(k_real) :: avg_scalar
  
  call MPIAverageGridScalarMasked(this%std_dev_grid, this%phase_fraction_grid.gt.0._k_real , avg_scalar)
  call csv_writer_obj%AppendScalar(avg_scalar, "stress_std_dev", write_headers)

  call this%writeAvgToCSVByModeAndTotal(csv_writer_obj, "abs_gamma_dot", abs(this%gammadot_grid), write_header=write_headers, write_total=.TRUE.)
  call this%writeAvgToCSVByModeAndTotal(csv_writer_obj, "rho_m", this%rho_m_grid, write_header=write_headers, write_total=.TRUE.)
  call this%writeAvgToCSVByModeAndTotal(csv_writer_obj, "delta_rho_m_generation", this%delta_rho_m_generation_grid, write_header=write_headers, write_total=.TRUE.)
  call this%writeAvgToCSVByModeAndTotal(csv_writer_obj, "delta_rho_m_dyn_recovery", this%delta_rho_m_dyn_recovery_grid, write_header=write_headers, write_total=.TRUE., scaling_factor=-1._k_real)
  call this%writeAvgToCSVByModeAndTotal(csv_writer_obj, "total_crss", this%crss_grid, write_header=write_headers, write_total=.FALSE.)

  if (this%use_static_recovery) then
    call this%writeAvgToCSVByModeAndTotal(csv_writer_obj, "delta_rho_m_static_recovery", this%delta_rho_m_static_recovery_grid, write_header=write_headers, write_total=.TRUE., scaling_factor=-1._k_real)
  endif

  call this%writeAvgToCSVByModeAndTotal(csv_writer_obj, "lambda", this%lambda_s_grid, write_header=write_headers, write_total=.FALSE.)
  call this%writeAvgToCSVByModeAndTotal(csv_writer_obj, "lambda_disl", this%lambda_dislocation_grid, write_header=write_headers, write_total=.FALSE.)
  call this%writeAvgToCSVByModeAndTotal(csv_writer_obj, "lambda_hall_petch", this%lambda_hall_petch_grid, write_header=write_headers, write_total=.FALSE.)
  call this%writeAvgToCSVByModeAndTotal(csv_writer_obj, "tau_line", this%tau_line_tension_grid, write_header=write_headers, write_total=.FALSE.)
  call this%writeAvgToCSVByModeAndTotal(csv_writer_obj, "tau_hall_petch", this%tau_hall_petch_grid, write_header=write_headers, write_total=.FALSE.)

  if (this%use_dislocation_cell_wall) then
    call this%writeAvgToCSVByModeAndTotal(csv_writer_obj, "delta_rho_m_trapped_cw", this%delta_rho_m_trapped_cw_grid, write_header=write_headers, write_total=.TRUE., scaling_factor=-1._k_real)
    call this%writeAvgToCSVByModeAndTotal(csv_writer_obj, "rho_cw", this%rho_cw_grid, write_header=write_headers, write_total=.TRUE.)
    call this%writeAvgToCSVByModeAndTotal(csv_writer_obj, "delta_rho_cw_trapped_m", this%delta_rho_m_trapped_cw_grid, write_header=write_headers, write_total=.TRUE.)
    call this%writeAvgToCSVByModeAndTotal(csv_writer_obj, "delta_rho_cw_annihilation", this%delta_rho_cw_annihilation_grid, write_header=write_headers, write_total=.TRUE., scaling_factor=-1._k_real)
    call this%writeAvgToCSVByModeAndTotal(csv_writer_obj, "lambda_cell_wall", this%sub_grain_size_grid, write_header=write_headers, write_total=.FALSE.)
    call this%writeAvgToCSVByModeAndTotal(csv_writer_obj, "tau_cell_wall", this%tau_cell_wall_grid, write_header=write_headers, write_total=.FALSE.)
  endif

  if (this%n_solute_type> 0) then
    call this%writeAvgToCSVByModeAndTotal(csv_writer_obj, "tau_M", this%tau_m_grid, write_header=write_headers, write_total=.FALSE.)
    call this%writeAvgToCSVByModeAndTotal(csv_writer_obj, "tau_labusch", this%tau_labusch_grid, write_header=write_headers, write_total=.FALSE.)
    call this%writeAvgToCSVByModeAndTotal(csv_writer_obj, "time_solute_aging_avg", this%time_solute_aging_avg_grid, write_header=write_headers, write_total=.FALSE.)
  endif
  if (this%n_prec_type> 0) then
    call this%writeAvgToCSVByModeAndTotal(csv_writer_obj, "tau_precipitate_drag", this%tau_precipitate_drag_grid, write_header=write_headers, write_total=.FALSE.)
    call this%writeAvgToCSVByModeAndTotal(csv_writer_obj, "lambda_prec", this%lambda_precipitate_grid, write_header=write_headers, write_total=.FALSE.)
  endif

  if (this%n_dislocation_loop_type> 0) then
    call this%writeAvgToCSVByModeAndTotal(csv_writer_obj, "tau_dislocation_loop", this%tau_disl_loop_grid, write_header=write_headers, write_total=.FALSE.)
  endif

end subroutine

end submodule
