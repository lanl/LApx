module glide_mod
use kinds
use glide_base_mod, only : cp_glide_base
use all_grid_data_mod, only : all_grid_data
use stiffness_base_mod, only : stiffness_base
use common_material_parameter_mod, only : common_material_parameter
use cp_util_classes_mod, only : slipmode_temperaure_piecewise_linear_parameter

implicit none

type, extends(cp_glide_base) :: hutchinson_glide
  real(k_real), pointer, dimension(:) :: gdot0 => null()
  real(k_real), pointer, dimension(:) :: n_exp => null()
  real(k_real), pointer, dimension(:) :: crss0 => null()
  real(k_real), pointer, dimension(:) :: tau1_hardening => null()
  real(k_real), pointer, dimension(:) :: theta0_hardening => null()
  real(k_real), pointer, dimension(:) :: theta1_hardening => null()
  real(k_real), pointer, dimension(:) :: crss_abs_rel_tol => null()
  real(k_real) :: rss_std_dev = 0.
  real(k_real), pointer, dimension(:,:,:,:) :: crss_grid => null() !-> the critical resolved shear stress
  real(k_real), pointer, dimension(:,:,:,:) :: crss_old_grid => null() !-> the old critical resolved shear stress
  real(k_real), pointer, dimension(:,:,:,:) :: crss_rate_grid => null() !-> the critical resolved shear stress
  real(k_real), pointer, dimension(:,:,:,:) :: crss_rate_old_grid => null() !-> the old critical resolved shear stress
  real(k_real), pointer, dimension(:,:,:) :: gamma_accumulated_grid => null() !-> the accumulated shear
  real(k_real), pointer, dimension(:,:,:) :: gamma_accumulated_old_grid => null() !-> the old accumulated shear
  real(k_real), pointer, dimension(:) :: crss_ptr => null()
  real(k_real), pointer, dimension(:) :: crss_old_ptr => null()
  real(k_real), pointer, dimension(:) :: crss_rate_ptr => null()
  real(k_real), pointer, dimension(:) :: crss_rate_old_ptr => null()
  real(k_real), pointer :: accumulated_gamma_ptr => null()
  real(k_real), pointer :: accumulated_gamma_old_ptr => null()
  logical :: use_hardening = .TRUE.
contains
  procedure :: initParameters => initParametersHutchinsonGlide
  procedure :: addFieldVariablesToGrid => addFieldVariablesToGridHutchinsonGlide
  procedure :: initGridPointers => initGridPointersHutchinsonGlide
  procedure :: initStateVariablesAtMaterialPoint=>initStateVariablesAtMaterialPointHutchinsonGlide
  procedure :: updateStateVariablesAtMaterialPointInnerLoop=>updateStateVariablesAtMaterialPointInnerLoopHutchinsonGlide
  procedure :: updateStateVariablesAtMaterialPointStaggered=>updateStateVariablesAtMaterialPointStaggeredHutchinsonGlide


  procedure :: setPointData=> setPointDataHutchinsonGlide
  procedure :: computeGammaDotAndDGammaDotDRssSS => computeGammaDotAndDGammaDotDRssSSHutchinsonGlide
  procedure :: readMaterialParametersFromFile => readMaterialParametersFromFileHutchinsonGlide

  procedure :: acceptRejectSolution => acceptRejectSolutionHutchinsonGlide
end type

interface
  module subroutine readMaterialParametersFromFileHutchinsonGlide(this, matf_reader)
    use read_from_file_utils, only : file_reader
    implicit none 
    class(hutchinson_glide), intent(inout) :: this
    class(file_reader), intent(inout) :: matf_reader
  end subroutine

  module subroutine initParametersHutchinsonGlide(this, phase_id, common_material_parameter_ptr, &
                                              use_damage, n_gauss, n_std_dev_gauss_integration, &
                                              elasticity_obj, crystal_paremeters_ptr )

    use stiffness_base_mod, only : stiffness_base
    use cp_base_mod, only : crystal_paremeters_type  
    implicit none
    class(hutchinson_glide), intent(inout) :: this
    integer, intent(in) :: phase_id
    class(common_material_parameter), pointer, intent(in) :: common_material_parameter_ptr
    real(k_real) :: ratio_edge_screw
    logical, intent(in) :: use_damage
    integer, intent(in):: n_gauss
    real(k_real), intent(in) :: n_std_dev_gauss_integration
    class(stiffness_base), pointer, intent(in) :: elasticity_obj
    class(crystal_paremeters_type), pointer, intent(in) :: crystal_paremeters_ptr

    end subroutine

    module subroutine initGridPointersHutchinsonGlide(this)
      implicit none
      class(hutchinson_glide), intent(inout) :: this
    end subroutine

    module subroutine initStateVariablesAtMaterialPointHutchinsonGlide(this)
      implicit none
      class(hutchinson_glide), intent(inout) :: this
    end subroutine

    module subroutine setPointDataHutchinsonGlide(this, ix,iy,iz)
      implicit none
      class(hutchinson_glide), intent(inout) :: this
      integer, intent(in) :: ix,iy,iz
    end subroutine

    module subroutine computeGammaDotAndDGammaDotDRssSSHutchinsonGlide(this, rss, gdot, dgdot_drss)
      implicit none
      class(hutchinson_glide), intent(inout) :: this
      real(k_real), intent(in) :: rss
      real(k_real), intent(out) :: gdot
      real(k_real), intent(out) :: dgdot_drss
    end subroutine

    module subroutine   addFieldVariablesToGridHutchinsonGlide(this)
      use grid_data_var_type_mod
      class(hutchinson_glide), intent(inout) :: this
    end subroutine

    module subroutine updateStateVariablesAtMaterialPointInnerLoopHutchinsonGlide(this)
      implicit none
      class(hutchinson_glide), intent(inout) :: this
    end subroutine

    module subroutine updateStateVariablesAtMaterialPointStaggeredHutchinsonGlide(this)
      implicit none
      class(hutchinson_glide), intent(inout) :: this
    end subroutine

    module subroutine acceptRejectSolutionHutchinsonGlide(this, dt_max, accept_solution_flag)
      implicit none
      class(hutchinson_glide), intent(inout) :: this
      real(k_real), intent(out) :: dt_max
      logical, intent(out) :: accept_solution_flag
    end subroutine

end interface

type, extends(cp_glide_base) :: wkkct_glide
  type(slipmode_temperaure_piecewise_linear_parameter) :: tau_labusch_linear_interp
  type(slipmode_temperaure_piecewise_linear_parameter) :: tau0_linear_interp
  ! all required grid pointers
  real(k_real), pointer, dimension(:,:,:,:,:) :: time_prec_climb_grid => null(), &
                                                 time_prec_climb_grid_old => null(), &
                                                 time_solute_aging_grid => null(), &
                                                 time_solute_aging_grid_old => null(), &
                                                 climb_velocity_grid => null(), &
                                                 climb_velocity_CellWall_grid => null(), &
                                                 prob_p_grid => null(), &
                                                 diameter_prec_grid => null(), &
                                                 tau_solute_drag_gauss_grid => null(), &
                                                 tau_solute_drag_gauss_grid_old => null(), &
                                                 tau_solute_drag_gauss_grid_older => null(), &
                                                 crss_gauss_grid => null()

  real(k_real), pointer, dimension(:,:,:,:) :: time_solute_aging_avg_grid => null(), &
                                               tau_line_tension_grid => null(), &
                                               tau_line_tension_grid_old => null(), &
                                               tau_line_tension_grid_older => null(), &
                                               tau_precipitate_drag_grid => null(), &
                                               tau_labusch_grid => null(), &
                                               tau_disl_loop_grid => null(), &
                                               tau_hall_petch_grid => null(), &
                                               tau_cell_wall_grid  => null(), &
                                               lambda_s_grid => null(), &
                                               lambda_dislocation_grid => null(), &
                                               lambda_precipitate_grid => null(), &
                                               lambda_hall_petch_grid => null(), &
                                               crss_grid => null(), &
                                               tau_m_grid =>null(), &
                                               sub_grain_size_grid => null(), &
                                               prob_d_grid => null(), &
                                               rho_m_grid => null(), &
                                               rho_m_grid_old => null(), &
                                               rho_m_dot_grid => null(), &
                                               rho_m_dot_grid_old => null(), &
                                               rho_cw_grid => null(), &
                                               rho_cw_grid_old => null(), &
                                               delta_rho_m_generation_grid => null(), &
                                               delta_rho_m_dyn_recovery_grid => null(), &
                                               delta_rho_m_trapped_cw_grid => null(), &
                                               delta_rho_m_static_recovery_grid => null(), &
                                               delta_rho_cw_annihilation_grid => null()
                                               
  ! all required material point pointers
  real(k_real), pointer, dimension(:,:) :: time_prec_climb_ptr => null(), &
                                           time_prec_climb_ptr_old => null(), &
                                           time_solute_aging_ptr => null(), &
                                           time_solute_aging_ptr_old => null(), &
                                           tau_solute_drag_gauss_ptr  => null(), &
                                           tau_solute_drag_gauss_ptr_old  => null(), &
                                           crss_gauss_ptr => null(), &
                                           climb_velocity_ptr => null(), &
                                           climb_velocity_CellWall_ptr => null(), &
                                           prob_p_ptr => null(), &
                                           deltaG0_prec_ptr => null(), & !-> DETG0P
                                           freq_prec_ptr => null(), & ! -> attfp
                                           beta_prec_ptr => null(), &          ! -> beta_prec
                                           alpha_prec_ptr => null(), &         ! -> alpha_prec
                                           diameter_prec_ptr => null(),&
                                           DW_ptr => null(), &  !-> DW for solute-cross-core term
                                           DHc_ptr => null(), & ! -> DHc for solute-cross-core term
                                           m_cd_ptr => null(), & ! -> Number of neighbors for solute-cross-core term
                                           freq_cd_ptr => null(), &  ! -> attempt frequency for solute-cross-core term
                                           norm_wbar_ptr => null(), &  ! -> normalized core width for solute-cross-core term
                                           N_atoms_disl_core_ptr => null(), &          ! -> Nunmbers of atom per dislocation core
                                           disl_core_width_ptr => null(), & ! w parameter ~7.5
                                           phi_cd_ptr => null() , &         ! -> phi_cd
                                           alpha_cd_ptr => null(), &          ! -> alpha_cd
                                           chi_cd_ptr => null(), &          ! -> chi_cd
                                           K_Labusch_ptr => null(), &          ! -> K_Labusch values for solute strengthening
                                           exponent_Labusch_ptr => null()          ! -> exponent values for Labusch solute strengthening

  real(k_real), pointer, dimension(:) :: time_solute_aging_avg_ptr  => null(), &
                                         tau_line_tension_ptr  => null(), &
                                         tau_line_tension_ptr_old  => null(), &
                                         tau_precipitate_drag_ptr  => null(), &
                                         tau_labusch_ptr  => null(), &
                                         tau_disl_loop_ptr  => null(), &
                                         tau_hall_petch_ptr  => null(), &
                                         tau_cell_wall_ptr  => null(), &
                                         crss_ptr => null(), &
                                         tau_m_ptr =>null(), &
                                         lambda_s_ptr  => null(), &
                                         lambda_dislocation_ptr => null(), &
                                         lambda_precipitate_ptr => null(), &
                                         lambda_hall_petch_ptr => null(), &
                                         sub_grain_size_ptr  => null(), &
                                         prob_d_ptr => null(), &
                                         rho_m_ptr => null(), &
                                         rho_m_dot_ptr => null(), &
                                         rho_cw_ptr => null(), &
                                         delta_rho_m_generation_ptr => null(), &
                                         delta_rho_m_dyn_recovery_ptr => null(), &
                                         delta_rho_m_trapped_cw_ptr => null(), &
                                         delta_rho_m_static_recovery_ptr => null(), &
                                         delta_rho_cw_annihilation_ptr => null()


  ! slip system dependent parameters                       ! OLD NAMES
  real(k_real), pointer, dimension(:) :: tauLabusch_T0 => null(), &     ! -> tauLabusch
                                         tauLabusch_alpha => null(), &     ! -> tauLabusch
                                         p_exp => null(), q_exp => null(), &   ! -> PP_1, PQ
                                         deltaG0_disl => null(),& ! -> DETG0D
                                         burgN => null(), & ! -> BURG_N
                                         tau0_T0 => null(), &        ! -> tau0m
                                         tau0_alpha => null(), &        ! -> tau0m
                                         chi_e => null(), & ! -> chi_e ENTROPY FACTOR
                                         n_crss => null(), &               ! -> n_crss
                                         edot0_dynrec => null(), &
                                         chi_dynrec => null(), &                ! -> chi_dynrec
                                         D_dynrec => null(), g_dynrec => null(), & ! -> D_dynrec, g_dynrec
                                         k_disl_gen => null(), k_disl_trap => null(), k_cell_wall_ann => null(), &   ! -> k_disl_gen, k_disl_trap, ak_c
                                         k1_static_recov => null(), k2_static_recov => null()


  ! precipitate, solute and dislocation loops specific parameters
  real(k_real), pointer, dimension(:) :: number_density_prec_ptr=> null(),& ! ->  precipitate number density
                                         initial_diameter_prec_ptr => null(),& !->  precipitate size (diameter)
                                         solute_concentration_ptr => null(),& !-> initial solute concentration
                                         number_density_disloc_loops_ptr => null(),& !-> Dislocation loops number density
                                         initial_diameter_disloc_loops_ptr => null() !-> Dislcoation loops diameter


  real(k_real) :: etav, &             ! -> eta
                  shear_mod, &
                  poisson, &
                  pre_aging_time, & 
                  grain_diameter, &
                  k_hall_petch

  logical :: use_dislocation_cell_wall = .FALSE.
  logical :: use_static_recovery = .FALSE.
  logical :: solutes_as_backstress = .FALSE.
  logical :: evolve_std_dev = .FALSE.
  logical :: use_jog_screw = .FALSE.
  real(k_real), pointer, dimension(:) :: rho0_m => null(), & 
                                           rho0_cw => null()

  integer :: n_prec_type, n_solute_type, n_dislocation_loop_type
  real(k_real), pointer, dimension(:,:) :: loop_crystallography_factor_ptr => null()

  real(k_real), pointer :: rho_mass => null(), &
                           avg_eq_strain_rate_ptr => null()

  ! time march parameters
  real(k_real) :: max_allowed_dd_percent_increment, &
                  max_allowed_tau_sol_drag_increment, &
                  max_allowed_tau_eline_increment

contains
  procedure :: initParameters => initParametersWKKCTGlide
  procedure :: addFieldVariablesToGrid => addFieldVariablesToGridWKKCTGlide
  procedure :: initGridPointers => initGridPointersWKKCTGlide
  procedure :: initStateVariablesAtMaterialPoint=>initStateVariablesAtMaterialPointWKKCTGlide

  procedure :: setPointData=> setPointDataWKKCTGlide
  procedure :: computeGammaDotAndDGammaDotDRssSS => computeGammaDotAndDGammaDotDRssSSWKKCTGlide
  procedure :: readMaterialParametersFromFile => readMaterialParametersFromFileWKKCTGlide

  ! support procedures
  procedure :: computeTauM
  procedure :: computeTauEffective
  procedure :: computeWaitTime
  procedure :: computeTotalWaitTime

  !procedure used for updating state variables
  procedure :: computeEffectivePrecipitateDiameter
  procedure :: computeDeltaTauPrecipitate
  procedure :: computeDeltaTauPrecipitateBKS
  procedure :: computeDeltaTauDislocationLoopsBKS
  procedure :: computeLabuschStrength
  procedure :: computeStrength
  procedure :: computeDeltaTauLineTension
  procedure :: computeDeltaTauSoluteDrag
  procedure :: computeEffectiveDislocationInterspacingAndProbability
  procedure :: computeSubGrainSize
  procedure :: computeStdDeviation
  procedure :: computeRhoDotGeneration
  procedure :: computeK2DynamicRecovery
  procedure :: computeRhoDotDynamicRecovery
  procedure :: computeRhoDotStaticRecovery
  procedure :: computeRhoDotTrappedCellWall
  procedure :: computeRhoDotAnnihilationCellWall
  procedure :: computeRhoDotAnnihilationCellWallClimbVelocity
  procedure :: computeDislocationDensity
  procedure :: computeShearVelocity
  procedure :: computeBindingEnergyNormalizingFactor
  procedure :: computeFrequencyDislocation
  procedure :: computeFrequencyPrecipitate

  ! the update we wnat to do within the inner loop
  procedure :: updateStateVariablesAtMaterialPointInnerLoop => updateStateVariablesAtMaterialPointInnerLoopWKKCTGlide
  ! the update we wnat to do within the inner loop
  procedure :: updateStateVariablesAtMaterialPointOuterLoop => updateStateVariablesAtMaterialPointOuterLoopWKKCTGlide
  ! the update we wnat to do after convergence
  procedure :: updateStateVariablesAtMaterialPointStaggered => updateStateVariablesAtMaterialPointStaggeredWKKCTGlide

  procedure :: computeDEcoreSaturation
  procedure :: computeCoreDiffusionTime
  procedure :: computeDEcore

  procedure :: computeEffectiveHardeningMatrix
  procedure :: acceptRejectSolution => acceptRejectSolutionWKKCTGlide
  procedure :: writeAverageQuantitiesCSV => writeAverageQuantitiesCSVWKKCTGlide

end type

interface
  module subroutine readMaterialParametersFromFileWKKCTGlide(this, matf_reader)
    use read_from_file_utils, only : file_reader
    implicit none
    class(wkkct_glide), intent(inout) :: this
    class(file_reader), intent(inout) :: matf_reader
  end subroutine

  module subroutine initParametersWKKCTGlide(this, phase_id, common_material_parameter_ptr, &
                                              use_damage, n_gauss, n_std_dev_gauss_integration, &
                                              elasticity_obj, crystal_paremeters_ptr)
    use stiffness_base_mod, only : stiffness_base
    use cp_base_mod, only : crystal_paremeters_type   
    implicit none
    class(wkkct_glide), intent(inout) :: this
    integer, intent(in) :: phase_id
    class(common_material_parameter), pointer, intent(in) :: common_material_parameter_ptr
    logical, intent(in) :: use_damage
    integer, intent(in):: n_gauss
    real(k_real), intent(in) :: n_std_dev_gauss_integration
    class(stiffness_base), pointer, intent(in) :: elasticity_obj
    class(crystal_paremeters_type), pointer, intent(in) :: crystal_paremeters_ptr
    end subroutine

  module subroutine initGridPointersWKKCTGlide(this)
    implicit none
    class(wkkct_glide), intent(inout) :: this
  end subroutine

  module subroutine initStateVariablesAtMaterialPointWKKCTGlide(this)
    implicit none
    class(wkkct_glide), intent(inout) :: this
  end subroutine

  module subroutine setPointDataWKKCTGlide(this, ix,iy,iz)
    implicit none
    class(wkkct_glide), intent(inout) :: this
    integer, intent(in) :: ix,iy,iz
  end subroutine

  module subroutine computeGammaDotAndDGammaDotDRssSSWKKCTGlide(this, rss, gdot, dgdot_drss)
    implicit none
    class(wkkct_glide), intent(inout) :: this
    real(k_real), intent(in) :: rss
    real(k_real), intent(out) :: gdot
    real(k_real), intent(out) :: dgdot_drss
  end subroutine

  ! support procedure
  module subroutine computeTauM(this, solute_id, taum)
    implicit none
    class(wkkct_glide), intent(in) :: this
    real(k_real), intent(out) :: taum
    integer, intent(in) :: solute_id
  end subroutine

  module subroutine computeTauEffective(this, rss, tau_effective, dtau_eff_drss)
    implicit none
    class(wkkct_glide), intent(in) :: this
    real(k_real), intent(in) :: rss
    real(k_real), intent(out) :: tau_effective, dtau_eff_drss
  end subroutine

  module subroutine computeWaitTime(this, rss, crss, frequency, dg0, t_wait, dtwait_drss, use_jog_screw)
    class(wkkct_glide), intent(in) :: this
    real(k_real), intent(in) :: rss, crss, frequency, dg0
    real(k_real), intent(out) :: t_wait, dtwait_drss
    logical, intent(in), optional :: use_jog_screw
  end subroutine

  module subroutine computeTotalWaitTime(this, rss, twait_total, dtwait_tot_drss)
    class(wkkct_glide), intent(inout) :: this
    real(k_real), intent(in) :: rss
    real(k_real), intent(out) :: twait_total, dtwait_tot_drss
  end subroutine

  module subroutine   addFieldVariablesToGridWKKCTGlide(this)
    use grid_data_var_type_mod
    class(wkkct_glide), intent(inout) :: this
  end subroutine

  module subroutine   updateStateVariablesAtMaterialPointInnerLoopWKKCTGlide(this)
    implicit none
    class(wkkct_glide), intent(inout) :: this
  end subroutine

  module subroutine   updateStateVariablesAtMaterialPointOuterLoopWKKCTGlide(this)
    implicit none
    class(wkkct_glide), intent(inout) :: this
  end subroutine

  module subroutine   updateStateVariablesAtMaterialPointStaggeredWKKCTGlide(this)
    class(wkkct_glide), intent(inout) :: this
  end subroutine

  module subroutine computeEffectivePrecipitateDiameter(this, eff_d_prec)
    implicit none
    class(wkkct_glide), intent(inout) :: this
    real(k_real), dimension(this%n_prec_type, this%n_ss), intent(out) :: eff_d_prec
  end subroutine

  module subroutine   computeDislocationDensity(this, rho_m, rho_cw)
    class(wkkct_glide), intent(in) :: this
    real(k_real), intent(inout), dimension(this%n_ss) :: rho_m, rho_cw
  end subroutine

  module function computeRhoDotGeneration(this, abs_gamma_dot) result(rho_dot)
    implicit none
    class(wkkct_glide), intent(in) :: this
    real(k_real), intent(in) :: abs_gamma_dot(this%n_ss)
    real(k_real) :: rho_dot(this%n_ss)
  end function

  module function computeK2DynamicRecovery(this, gdot) result(k2)
    implicit none
    class(wkkct_glide), intent(in) :: this
    real(k_real), intent(in) :: gdot(this%n_ss)
    real(k_real) :: k2(this%n_ss)
  end function

  module function computeRhoDotDynamicRecovery(this, abs_gamma_dot) result(rho_dot)
    implicit none
    class(wkkct_glide), intent(in) :: this
    real(k_real), intent(in) :: abs_gamma_dot(this%n_ss)
    real(k_real) :: rho_dot(this%n_ss)
  end function

  module function computerhodotstaticrecovery(this) result(rho_dot)
    implicit none
    class(wkkct_glide), intent(in) :: this
    real(k_real) :: rho_dot(this%n_ss)
  end function

  module function computeRhoDotTrappedCellWall(this, abs_gamma_dot) result(rho_dot)
    implicit none
    class(wkkct_glide), intent(in) :: this
    real(k_real), intent(in) :: abs_gamma_dot(this%n_ss)
    real(k_real) :: rho_dot(this%n_ss)
  end function

  module function computeRhoDotAnnihilationCellWall(this, abs_gamma_dot) result(rho_dot)
    implicit none
    class(wkkct_glide), intent(in) :: this
    real(k_real), intent(in) :: abs_gamma_dot(this%n_ss)
    real(k_real) :: rho_dot(this%n_ss)
  end function

  module function computeRhoDotAnnihilationCellWallClimbVelocity(this) result(rho_dot)
    implicit none
    class(wkkct_glide), intent(in) :: this
    real(k_real) :: rho_dot(this%n_ss)
  end function

  module subroutine computeStrength(this, crss)
    implicit none
    class(wkkct_glide), intent(inout) :: this
    real(k_real), dimension(this%n_gauss, this%n_ss), intent(inout) :: crss
  end subroutine

  module subroutine   computeDeltaTauPrecipitate(this, prec_id, dtau_prec)
    implicit none
    class(wkkct_glide), intent(in) :: this
    integer, intent(in) :: prec_id
    real(k_real), intent(out),  dimension(this%n_ss) :: dtau_prec
  end subroutine

  module subroutine   computeDeltaTauPrecipitateBKS(this, prec_id, dtau_prec)
    implicit none
    class(wkkct_glide), intent(in) :: this
    integer, intent(in) :: prec_id
    real(k_real), intent(out),  dimension(this%n_ss) :: dtau_prec
  end subroutine

  module subroutine   computeDeltaTauDislocationLoopsBKS(this, loop_id, dtau_disloc_loop)
    implicit none
    class(wkkct_glide), intent(in) :: this
    integer, intent(in) :: loop_id
    real(k_real), intent(out),  dimension(this%n_ss) :: dtau_disloc_loop
  end subroutine

  module subroutine   computeLabuschStrength(this, solute_id, tau_Labsuch)
    implicit none
    class(wkkct_glide), intent(in) :: this
    integer, intent(in) :: solute_id
    real(k_real), intent(out),  dimension(this%n_ss) :: tau_Labsuch
  end subroutine

  module subroutine computeDeltaTauLineTension(this, alpha_ss_mobile, tau_line_tension)
    implicit none
    class(wkkct_glide), intent(in) :: this
    real(k_real), dimension(this%n_ss, this%n_ss), intent(in) :: alpha_ss_mobile
    real(k_real), dimension(this%n_ss), intent(out) :: tau_line_tension
  end subroutine

  module subroutine computeDeltaTauSoluteDrag(this, tau_sol_drag)
    implicit none
    class(wkkct_glide), intent(inout) :: this
    real(k_real), dimension(this%igauss, this%n_ss), intent(out) :: tau_sol_drag
  end subroutine

  module subroutine computeEffectiveDislocationInterspacingAndProbability(this, alpha_ss_mobile, alpha_ss_cw, lambda_s, prob_d, prob_p)
    implicit none
    class(wkkct_glide), intent(in) :: this
    real(k_real), dimension(this%n_ss, this%n_ss), intent(in) :: alpha_ss_mobile, alpha_ss_cw
    real(k_real), dimension(this%n_ss), intent(out) :: lambda_s, prob_d
    real(k_real), dimension(this%n_prec_type, this%n_ss), intent(out) :: prob_p
  end subroutine

  module subroutine computeSubGrainSize(this, alpha_ss_cw, sub_grain_size)
    implicit none
    class(wkkct_glide), intent(in) :: this
    real(k_real), dimension(this%n_ss, this%n_ss), intent(in) ::  alpha_ss_cw
    real(k_real), dimension(this%n_ss), intent(out) :: sub_grain_size
  end subroutine

  module subroutine computeStdDeviation(this, std_dev)
    implicit none
    class(wkkct_glide), intent(in) :: this
    real(k_real), intent(out) :: std_dev
  end subroutine

  module function computeDEcoreSaturation(this, solute_id, ss_idx) result(DEcore_sat)
    implicit none
    class(wkkct_glide), intent(in) :: this
    integer, intent(in) :: ss_idx, solute_id
    real(k_real) :: DEcore_sat
  end function

  module function computeCoreDiffusionTime(this, solute_id, ss_idx) result(td)
    implicit none
    class(wkkct_glide), intent(in) :: this
    integer, intent(in) :: ss_idx, solute_id
    real(k_real) :: td
  end function

  module function computeDEcore(this, solute_id, ss_idx) result(DEcore)
    implicit none
    class(wkkct_glide), intent(in) :: this
    real(k_real) :: DEcore
    integer, intent(in) :: ss_idx, solute_id
  end function

  module function computeShearVelocity(this) result(shear_velocity)
    implicit none
    class(wkkct_glide), intent(in) :: this
    real(k_real) :: shear_velocity
  end function

  module function computeBindingEnergyNormalizingFactor(this) result(A)
    implicit none
    class(wkkct_glide), intent(in) :: this
    real(k_real), dimension(this%n_ss) :: A
  end function

  module function computeFrequencyDislocation(this) result(disl_freq)
    implicit none
    class(wkkct_glide), intent(in) :: this
    real(k_real) :: disl_freq
  end function

  module function computeFrequencyPrecipitate(this, prec_id) result(prec_freq)
    implicit none
    class(wkkct_glide), intent(in) :: this
    integer, intent(in) :: prec_id
    real(k_real) :: prec_freq
  end function

  module subroutine computeEffectiveHardeningMatrix(this, alpha_ss)
    implicit none
    class(wkkct_glide), intent(inout) :: this
    real(k_real), intent(out) :: alpha_ss(this%n_ss, this%n_ss)
  end subroutine

  module subroutine acceptRejectSolutionWKKCTGlide(this, dt_max, accept_solution_flag)
    implicit none
    class(wkkct_glide), intent(inout) :: this
    real(k_real), intent(out) :: dt_max
    logical, intent(out) :: accept_solution_flag
  end subroutine

  module subroutine writeAverageQuantitiesCSVWKKCTGlide(this, csv_writer_obj, write_headers)
    use csv_writer_mod, only : csv_writer
    implicit none
    class(wkkct_glide), intent(inout) :: this
    class(csv_writer), intent(inout) :: csv_writer_obj
    logical, intent(in) :: write_headers
  end subroutine

end interface

contains
subroutine readMaterialParametersFromFileCPGlide(matf_reader, phase_id, & 
            use_damage, n_gauss, n_std_dev_gauss_integration, &
            all_mighty_grid_in,  sim_all_macro_data, the_bc_object, stiffness_ptr, &
            common_material_parameter_ptr, crystal_paremeters_ptr,&
            inelastic_strain_base_ptr)
  use read_from_file_utils, only : file_reader
  use all_mighty_grid_mod, only : all_mighty_grid_type
  use simumaltion_macro_data_mod, only : sim_all_macro_data_obj
  use bc_objects_mod, only : boundary_condition_array_type
  use common_material_parameter_mod, only : common_material_parameter
  use inelastic_strain_mod, only : inelastic_strain_base
  use string_module, only : string_type
  use stiffness_base_mod, only : stiffness_base
  use cp_base_mod, only : crystal_paremeters_type
  implicit none
  type(file_reader), intent(inout) :: matf_reader
  type(string_type) :: material_model_name
  class(inelastic_strain_base), pointer, intent(inout) :: inelastic_strain_base_ptr
  class(common_material_parameter), pointer, intent(inout) :: common_material_parameter_ptr
  class(crystal_paremeters_type), pointer, intent(inout) :: crystal_paremeters_ptr
  integer, intent(in) :: phase_id
  class(all_mighty_grid_type), intent(in), target :: all_mighty_grid_in
  class(sim_all_macro_data_obj), intent(in), target :: sim_all_macro_data
  class(boundary_condition_array_type), intent(in), target :: the_bc_object
  class(stiffness_base), pointer, intent(in) :: stiffness_ptr
  
  class(hutchinson_glide), pointer :: hutchinson_glide_temp => null()
  class(wkkct_glide), pointer :: wkkct_glide_temp => null()
  logical, intent(in) :: use_damage
  integer, intent(in) :: n_gauss
  real(k_real), intent(in) :: n_std_dev_gauss_integration

  call matf_reader%readParameter("--Glide-model", material_model_name)
  select case(material_model_name%getString())
  case ("hutchinson-glide")
    allocate(hutchinson_glide_temp)
    inelastic_strain_base_ptr => hutchinson_glide_temp
    call hutchinson_glide_temp%linkBaseObjects(phase_id, all_mighty_grid_in, sim_all_macro_data, the_bc_object)
    call hutchinson_glide_temp%initParameters(phase_id, common_material_parameter_ptr, &
                                            use_damage, n_gauss, n_std_dev_gauss_integration, &
                                            stiffness_ptr, crystal_paremeters_ptr)
  case ("wkkct-glide")
    allocate(wkkct_glide_temp)
    inelastic_strain_base_ptr => wkkct_glide_temp
    call wkkct_glide_temp%linkBaseObjects(phase_id, all_mighty_grid_in, sim_all_macro_data, the_bc_object)
    call wkkct_glide_temp%initParameters(phase_id, common_material_parameter_ptr, &
                                          use_damage, n_gauss, n_std_dev_gauss_integration, &
                                          stiffness_ptr, crystal_paremeters_ptr)

  case default
   write(*,*) "unrecognized material model for glide", material_model_name%getString()
   write(*,*) "available material models are: "
   write(*,*) "wkkct-glide "
   write(*,*) "hutchinson-glide "
    error stop "abort"
  end select

  call inelastic_strain_base_ptr%readMaterialParametersFromFile(matf_reader)

  select case(material_model_name%getString())
  case ("hutchinson-glide")
    nullify(hutchinson_glide_temp)
  case ("wkkct-glide")
    nullify(wkkct_glide_temp)
  case default
    write(*,*) "somthing is wrong when trying nullifying the pointer for ", material_model_name%getString()
    error stop "abort"
  end select

end subroutine

end module
