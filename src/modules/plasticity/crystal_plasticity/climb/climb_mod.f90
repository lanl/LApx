module climb_mod
use kinds
use climb_base_mod, only : cp_climb_base
use common_material_parameter_mod, only : common_material_parameter
use stiffness_base_mod, only : stiffness_base

implicit none

type, extends(cp_climb_base) :: exponential_climb
  real(k_real), pointer, dimension(:) :: bdot0 => null()
  real(k_real), pointer, dimension(:) :: n_exp => null()
  real(k_real), pointer, dimension(:) :: crss0 => null()
  real(k_real) :: rss_std_dev = 0.
  real(k_real), pointer, dimension(:,:,:,:) :: crss_grid => null() !-> the critical resolved shear stress
  real(k_real), pointer, dimension(:) :: crss_ptr => null()

contains
  procedure :: initParameters => initParametersExponentialClimb
  procedure :: addFieldVariablesToGrid => addFieldVariablesToGridExponentialClimb
  procedure :: initGridPointers => initGridPointersExponentialClimb
  procedure :: initStateVariablesAtMaterialPoint=>initStateVariablesAtMaterialPointExponentialClimb
  procedure :: updateStateVariablesAtMaterialPointInnerLoop=>updateStateVariablesAtMaterialPointInnerLoopExponentialClimb
  procedure :: updateStateVariablesAtMaterialPointStaggered=>updateStateVariablesAtMaterialPointStaggeredExponentialClimb


  procedure :: setPointData=> setPointDataExponentialClimb
  procedure :: computeGammaDotAndDGammaDotDRssSS => computeGammaDotAndDGammaDotDRssSSExponentialClimb
  procedure :: readMaterialParametersFromFile => readMaterialParametersFromFileExponentialClimb

  procedure :: acceptRejectSolution => acceptRejectSolutionExponentialClimb

end type

interface
  module subroutine readMaterialParametersFromFileExponentialClimb(this, matf_reader)
    use read_from_file_utils, only : file_reader
    class(exponential_climb), intent(inout) :: this
    class(file_reader), intent(inout) :: matf_reader
  end subroutine

  module subroutine initParametersExponentialClimb(this, phase_id, common_material_parameter_ptr,&
                                                  use_damage, n_gauss, n_std_dev_gauss_integration, &
                                                   elasticity_obj, crystal_paremeters_ptr)
    use stiffness_base_mod, only : stiffness_base         
    use cp_base_mod, only : crystal_paremeters_type                              
    implicit none
    class(exponential_climb), intent(inout) :: this
    integer, intent(in) :: phase_id
    class(common_material_parameter), pointer, intent(in) :: common_material_parameter_ptr
    logical, intent(in) :: use_damage
    integer, intent(in):: n_gauss
    real(k_real), intent(in) :: n_std_dev_gauss_integration
    class(stiffness_base), pointer, intent(in) :: elasticity_obj
    class(crystal_paremeters_type), pointer, intent(in) :: crystal_paremeters_ptr

    end subroutine

    module subroutine initGridPointersExponentialClimb(this)
      implicit none
      class(exponential_climb), intent(inout) :: this
    end subroutine

    module subroutine initStateVariablesAtMaterialPointExponentialClimb(this)
      implicit none
      class(exponential_climb), intent(inout) :: this
    end subroutine

    module subroutine setPointDataExponentialClimb(this, ix,iy,iz)
      implicit none
      class(exponential_climb), intent(inout) :: this
      integer, intent(in) :: ix,iy,iz
    end subroutine

    module subroutine computeGammaDotAndDGammaDotDRssSSExponentialClimb(this, rss, gdot, dgdot_drss)
      implicit none
      class(exponential_climb), intent(inout) :: this
      real(k_real), intent(in) :: rss
      real(k_real), intent(out) :: gdot
      real(k_real), intent(out) :: dgdot_drss
    end subroutine

    module subroutine   addFieldVariablesToGridExponentialClimb(this)
      use grid_data_var_type_mod
      class(exponential_climb), intent(inout) :: this
    end subroutine

    module subroutine updateStateVariablesAtMaterialPointInnerLoopExponentialClimb(this)
      implicit none
      class(exponential_climb), intent(inout) :: this
    end subroutine

    module subroutine updateStateVariablesAtMaterialPointStaggeredExponentialClimb(this)
      implicit none
      class(exponential_climb), intent(inout) :: this
    end subroutine

    module subroutine acceptRejectSolutionExponentialClimb(this, dt_max, accept_solution_flag)
      implicit none
      class(exponential_climb), intent(inout) :: this
      real(k_real), intent(out) :: dt_max
      logical, intent(out) :: accept_solution_flag
    end subroutine

end interface


type, extends(cp_climb_base) :: transient_climb
  real(k_real), pointer, dimension(:) :: bdot0 => null()
  real(k_real), pointer, dimension(:) :: n_exp => null()
  real(k_real), pointer, dimension(:) :: crss0 => null()
  real(k_real) :: rss_std_dev = 0.
  real(k_real), pointer ::rho_m_grid(:,:,:,:) => null(), &
                          climb_velocity_grid(:,:,:,:,:)  => null(), &
                          climb_velocity_CellWall_grid(:,:,:,:,:)  => null(), &
                          vacancy_conc_grid(:,:,:)  => null(), &
                          interstitial_conc_grid(:,:,:) => null(), & 
                          interstitial_cluster_conc_grid(:,:,:,:) => null(), &
                          vacancy_conc_grid_old(:,:,:)  => null(), &
                          interstitial_conc_grid_old(:,:,:) => null(), & 
                          interstitial_cluster_conc_grid_old(:,:,:,:) => null(), &
                          vacancy_loop_density_grid(:,:,:,:) => null(), &
                          vacancy_loop_radii_grid(:,:,:,:) => null(), &
                          interstitial_loop_density_grid(:,:,:,:) => null(), &
                          interstitial_loop_radii_grid(:,:,:,:) => null(), &
                          dislocation_climb_strain_rate_grid(:,:,:,:) => null(), &
                          dislocation_climb_strain_grid(:,:,:,:) => null(), &
                          dislocation_climb_strain_grid_old(:,:,:,:) => null(), &
                          vacancy_loop_strain_rate_grid(:,:,:,:) => null(), &
                          vacancy_loop_strain_grid(:,:,:,:) => null(), &
                          vacancy_loop_strain_grid_old(:,:,:,:) => null(), &
                          interstitial_loop_strain_rate_grid(:,:,:,:) => null(), &
                          interstitial_loop_strain_grid(:,:,:,:) => null(), &
                          interstitial_loop_strain_grid_old(:,:,:,:) => null(),&
                          loops_climb_strain_rate_grid(:,:,:,:) => null(), &
                          loops_climb_strain_grid(:,:,:,:) => null(), &
                          loops_climb_strain_grid_old(:,:,:,:) => null(), &
                          dose_a_grid(:,:,:) => null(), &
                          dose_c_grid(:,:,:) => null(), &
                          dose_a_grid_old(:,:,:) => null(), &
                          dose_c_grid_old(:,:,:) => null()

  real(k_real), pointer ::rho_m_ptr(:) => null(), &
                          climb_velocity_ptr(:,:)  => null(), &
                          climb_velocity_CellWall_ptr(:,:)  => null(), &
                          vacancy_conc_ptr  => null(), &
                          interstitial_conc_ptr => null(), & 
                          interstitial_cluster_conc_ptr(:) => null(), &
                          vacancy_conc_ptr_old  => null(), &
                          interstitial_conc_ptr_old => null(), & 
                          interstitial_cluster_conc_ptr_old(:) => null(), &
                          vacancy_loop_density_ptr(:) => null(), &
                          vacancy_loop_radii_ptr(:) => null(), &
                          interstitial_loop_density_ptr(:) => null(), &
                          interstitial_loop_radii_ptr(:) => null(), &
                          dislocation_climb_strain_rate_ptr(:) => null(), &
                          dislocation_climb_strain_ptr(:) => null(), &
                          dislocation_climb_strain_ptr_old(:) => null(), &
                          vacancy_loop_strain_rate_ptr(:) => null(), &
                          vacancy_loop_strain_ptr(:) => null(), &
                          vacancy_loop_strain_ptr_old(:) => null(), &
                          interstitial_loop_strain_rate_ptr(:) => null(), &
                          interstitial_loop_strain_ptr(:) => null(), &
                          interstitial_loop_strain_ptr_old(:) => null(), &
                          loops_climb_strain_rate_ptr(:) => null(), &
                          loops_climb_strain_ptr(:) => null(), &
                          loops_climb_strain_ptr_old(:) => null(), &
                          dose_a_ptr => null(), &
                          dose_c_ptr => null(), &
                          dose_a_ptr_old => null(), &
                          dose_c_ptr_old => null()

! material parameters
real(k_real), pointer, dimension(:) :: E_migration_interstitial => null(), &
                                       E_migration_interstitial_cluster => null(), &
                                       interstitial_diffusivity_coefficient => null(), &
                                       interstitial_cluster_diffusivity_coefficient => null()

real(k_real) :: atomic_volume, voxel_volume
! concentration solver variables
real(k_real), pointer :: K2(:,:)=>null(), &
                         c0(:,:)=>null(), &
                         D(:) =>null(), &
                         G(:) => null(), &
                         cdot(:) =>null(), &
                         dcdot_dc(:,:) => null(), &
                         I_conc(:,:) => null(), &
                         Residual(:) =>null(), &
                         dResidualdC(:,:) =>null(), &
                         c_guess(:) =>null(), &
                         delta_c_guess(:) =>null(), &
                         c(:) =>null(), &
                         c_old(:) =>null(), &
                         n_atoms_per_defect(:) => null(), &
                         n_stored_atoms(:) => null(), &
                         n0_stored_atoms(:,:) => null(), &
                         n_dot(:) => null()

contains
  procedure :: initParameters => initParametersTransientClimb
  procedure :: addFieldVariablesToGrid => addFieldVariablesToGridTransientClimb
  procedure :: initGridPointers => initGridPointersTransientClimb
  procedure :: initStateVariablesAtMaterialPoint=>initStateVariablesAtMaterialPointTransientClimb
  procedure :: updateStateVariablesAtMaterialPointOuterLoop=>updateStateVariablesAtMaterialPointOuterLoopTransientClimb
  procedure :: updateStateVariablesAtMaterialPointStaggered=>updateStateVariablesAtMaterialPointStaggeredTransientClimb


  procedure :: setPointData=> setPointDataTransientClimb
  procedure :: computeGammaDotAndDGammaDotDRssSS => computeGammaDotAndDGammaDotDRssSSTransientClimb
  procedure :: readMaterialParametersFromFile => readMaterialParametersFromFileTransientClimb

  procedure :: acceptRejectSolution => acceptRejectSolutionTransientClimb
  procedure :: getStrainRateaAndStressJacobian => getStrainRateaAndStressJacobianTransientClimb !-> the procedure
  procedure :: writeInelastiStrainRateToFile => writeInelastiStrainRateToFileTransientClimb

end type

interface

  module subroutine writeInelastiStrainRateToFileTransientClimb(this, csv_writer_obj, write_headers)
    use csv_writer_mod, only : csv_writer
    implicit none
    class(transient_climb), intent(inout) :: this
    class(csv_writer), intent(inout) :: csv_writer_obj
    logical, intent(in) :: write_headers
  end subroutine

  module subroutine readMaterialParametersFromFileTransientClimb(this, matf_reader)
    use read_from_file_utils, only : file_reader
    implicit none
    class(transient_climb), intent(inout) :: this
    class(file_reader), intent(inout) :: matf_reader
  end subroutine

  module subroutine initParametersTransientClimb(this, phase_id, common_material_parameter_ptr,&
                                                use_damage, n_gauss, n_std_dev_gauss_integration, &
                                                elasticity_obj, crystal_paremeters_ptr)
    use stiffness_base_mod, only : stiffness_base         
    use cp_base_mod, only : crystal_paremeters_type                              
    implicit none
    class(transient_climb), intent(inout) :: this
    integer, intent(in) :: phase_id
    class(common_material_parameter), pointer, intent(in) :: common_material_parameter_ptr
    logical, intent(in) :: use_damage
    integer, intent(in):: n_gauss
    real(k_real), intent(in) :: n_std_dev_gauss_integration
    class(stiffness_base), pointer, intent(in) :: elasticity_obj
    class(crystal_paremeters_type), pointer, intent(in) :: crystal_paremeters_ptr
    
  end subroutine

    module subroutine initGridPointersTransientClimb(this)
      implicit none
      class(transient_climb), intent(inout) :: this
    end subroutine

    module subroutine initStateVariablesAtMaterialPointTransientClimb(this)
      implicit none
      class(transient_climb), intent(inout) :: this
    end subroutine

    module subroutine setPointDataTransientClimb(this, ix,iy,iz)
      implicit none
      class(transient_climb), intent(inout) :: this
      integer, intent(in) :: ix,iy,iz
    end subroutine

    module subroutine computeGammaDotAndDGammaDotDRssSSTransientClimb(this, rss, gdot, dgdot_drss)
      implicit none
      class(transient_climb), intent(inout) :: this
      real(k_real), intent(in) :: rss
      real(k_real), intent(out) :: gdot
      real(k_real), intent(out) :: dgdot_drss
    end subroutine

    module subroutine   addFieldVariablesToGridTransientClimb(this)
      use grid_data_var_type_mod
      class(transient_climb), intent(inout) :: this
    end subroutine

    module subroutine updateStateVariablesAtMaterialPointOuterLoopTransientClimb(this)
      implicit none
      class(transient_climb), intent(inout) :: this
    end subroutine

    module subroutine updateStateVariablesAtMaterialPointStaggeredTransientClimb(this)
      implicit none
      class(transient_climb), intent(inout) :: this
    end subroutine

    module subroutine acceptRejectSolutionTransientClimb(this, dt_max, accept_solution_flag)
      implicit none
      class(transient_climb), intent(inout) :: this
      real(k_real), intent(out) :: dt_max
      logical, intent(out) :: accept_solution_flag
    end subroutine

    module subroutine getStrainRateaAndStressJacobianTransientClimb(this, stress6, epsilon_dot, depsilon_dot_dstress, ix, iy, iz)
      implicit none
      class(transient_climb), intent(inout) :: this
      real(k_real), target, intent(in) ::stress6(6)
      real(k_real), dimension(6), intent(out) :: epsilon_dot
      real(k_real), dimension(6,6), intent(out) :: depsilon_dot_dstress
      integer, intent(in) ::  ix, iy, iz
    end subroutine

end interface

!******************************************************************************!
! WKKCT climb base type. This model is based on Wen et al. 2020
! This class is abstract, and cannot be used. The climbvelcotity function is empty
! When extending this class one should also remember to add the appropiate parameters,
! field variables and pointers
!******************************************************************************!
type, extends(cp_climb_base) :: wkkct_climb_base
  real(k_real), pointer, dimension(:) :: burgN => null()

  ! girddata pointers
  real(k_real), pointer, dimension(:,:,:,:,:) :: climb_velocity_grid => null()
  real(k_real), pointer, dimension(:,:,:,:,:) :: climb_velocity_CellWall_grid => null()
  real(k_real), pointer, dimension(:,:,:,:) :: rho_m_grid => null()
  real(k_real), pointer, dimension(:,:,:,:) :: rho_cw_grid => null()
  real(k_real), pointer, dimension(:,:,:,:) :: sub_grain_size_grid => null()

  real(k_real), pointer, dimension(:,:) :: climb_velocity_ptr => null()
  real(k_real), pointer, dimension(:,:) :: climb_velocity_CellWall_ptr => null()
  real(k_real), pointer, dimension(:) :: rho_m_ptr => null()
  real(k_real), pointer, dimension(:) :: rho_cw_ptr => null()
  real(k_real), pointer, dimension(:) :: sub_grain_size_ptr  => null()
  real(k_real) :: climb_reduce_ratio

  !! paramters
  real(k_real) :: zi0, Di0, zv0, E_migration_interstitial, shear_mod

contains
  ! basic procedures we need to override
  procedure :: initParametersWKKCTClimbBase
  procedure :: addFieldVariablesToGrid => addFieldVariablesToGridWKKCTClimbBase
  procedure :: initGridPointers => initGridPointersWKKCTClimbBase
  procedure :: setPointData => setPointDataWKKCTClimbBase
  procedure :: initStateVariablesAtMaterialPoint => initStateVariablesAtMaterialPointWKKCTClimbBase
  procedure :: updateStateVariablesAtMaterialPointInnerLoop => updateStateVariablesAtMaterialPointInnerLoopWKKCTClimbBase
  procedure :: updateStateVariablesAtMaterialPointStaggered => updateStateVariablesAtMaterialPointStaggeredWKKCTClimbBase

  procedure :: computeClimbVelocity => computeClimbVelocityWKKCTClimbBase
  procedure :: computeClimbVelocityCellWall => computeClimbVelocityCellWallWKKCTClimbBase
  procedure :: computeZinterstitial => computeZinterstitialWKKCTClimbBase
  procedure :: computeZVacancy => computeZVacancyWKKCTClimbBase
  procedure :: computeGammaDotAndDGammaDotDRssSS=> computeGammaDotAndDGammaDotDRssSSWKKCTClimbBase
  procedure :: computeConcentrationInterstiatialBulk => computeConcentrationInterstiatialBulkWKKCTBase
  procedure :: computeConcentrationVacancyBulk => computeConcentrationVacancyBulkWKKCTBase
  procedure :: computeConcentrationVacancyCellWall => computeConcentrationVacancyCellWallWKKCTBase
  procedure :: computeConcentrationVacancyCore => computeConcentrationVacancyCoreWKKCTBase
  procedure :: computeInterstialDiffusivity
  procedure :: resolved_stress_equilibrium_concentration
  procedure :: cellwall_local_concentration
  procedure :: hydrostatic_stress_equilibrium_concentration
end type

interface

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
      

  end subroutine

  module subroutine  addFieldVariablesToGridWKKCTClimbBase(this)
    implicit none
    class(wkkct_climb_base), intent(inout) :: this
  end subroutine

  module subroutine initGridPointersWKKCTClimbBase(this)
    implicit none
    class(wkkct_climb_base), intent(inout) :: this
  end subroutine

  module subroutine setPointDataWKKCTClimbBase(this, ix,iy,iz)
    implicit none
    class(wkkct_climb_base), intent(inout) :: this
    integer, intent(in) :: ix,iy,iz
  end subroutine

  module subroutine initStateVariablesAtMaterialPointWKKCTClimbBase(this)
    implicit none
    class(wkkct_climb_base), intent(inout) :: this
  end subroutine

  module subroutine updateStateVariablesAtMaterialPointInnerLoopWKKCTClimbBase(this)
    implicit none
    class(wkkct_climb_base), intent(inout) :: this
  end subroutine

  module subroutine updateStateVariablesAtMaterialPointStaggeredWKKCTClimbBase(this)
    implicit none
    class(wkkct_climb_base), intent(inout) :: this
  end subroutine

  module subroutine computeClimbVelocityWKKCTClimbBase(this, rss, vclimb, dvclimb_dtauclimb)
    implicit none
    class(wkkct_climb_base), intent(in) :: this
    real(k_real), intent(in) :: rss
    real(k_real), intent(out) :: vclimb, dvclimb_dtauclimb
  end subroutine

  module subroutine computeClimbVelocityCellWallWKKCTClimbBase(this, rss, vclimb_CW, dvclimb_dtauclimb_CW)
    implicit none
    class(wkkct_climb_base), intent(in) :: this
    real(k_real), intent(in) :: rss
    real(k_real), intent(out) :: vclimb_CW, dvclimb_dtauclimb_CW
  end subroutine

  module subroutine computeZinterstitialWKKCTClimbBase(this, tau_climb, zi, dzi_dtauclimb)
    implicit none
    class(wkkct_climb_base), intent(in) :: this
    real(k_real), intent(in) :: tau_climb
    real(k_real), intent(out) :: zi, dzi_dtauclimb
  end subroutine

  module subroutine computeZVacancyWKKCTClimbBase(this, tau_climb, zv, dzv_dtauclimb)
    implicit none
    class(wkkct_climb_base), intent(in) :: this
    real(k_real), intent(in) :: tau_climb
    real(k_real), intent(out) :: zv, dzv_dtauclimb
  end subroutine

  module subroutine computeGammaDotAndDGammaDotDRssSSWKKCTClimbBase(this, rss, gdot, dgdot_drss)
    implicit none
    class(wkkct_climb_base), intent(inout) :: this
    real(k_real), intent(in) :: rss
    real(k_real), intent(out) :: gdot
    real(k_real), intent(out) :: dgdot_drss
  end subroutine

  module subroutine computeConcentrationInterstiatialBulkWKKCTBase(this, tau_climb, ci_bulk, dci_bulk_dtauclimb)
    implicit none
    class(wkkct_climb_base), intent(in) :: this
    real(k_real), intent(in) :: tau_climb
    real(k_real), intent(out) :: ci_bulk, dci_bulk_dtauclimb
  end subroutine

  module subroutine computeConcentrationVacancyBulkWKKCTBase(this, tau_climb, cv_bulk, dcv_bulk_dtauclimb)
    implicit none
    class(wkkct_climb_base), intent(in) :: this
    real(k_real), intent(in) :: tau_climb
    real(k_real), intent(out) :: cv_bulk, dcv_bulk_dtauclimb
  end subroutine

  module subroutine computeConcentrationVacancyCellWallWKKCTBase(this, rho_cw, cv_CW, dcv_CW_dtauclimb)
    implicit none
    class(wkkct_climb_base), intent(in) :: this
    real(k_real), intent(in) :: rho_cw
    real(k_real), intent(out) :: cv_CW, dcv_CW_dtauclimb
  end subroutine

  module subroutine computeConcentrationVacancyCoreWKKCTBase(this, tau_climb, cv_core, dcv_core_dtauclimb)
    implicit none
    class(wkkct_climb_base), intent(in) :: this
    real(k_real), intent(in) :: tau_climb
    real(k_real), intent(out) :: cv_core, dcv_core_dtauclimb
  end subroutine

  module subroutine computeInterstialDiffusivity(this, Di)
    class(wkkct_climb_base), intent(in) :: this
    real(k_real), intent(out) :: Di
  end subroutine

  module subroutine resolved_stress_equilibrium_concentration(this, tau_climb, ctau, dctau_dtauclimb)
    implicit none
    class(wkkct_climb_base), intent(in) :: this
    real(k_real), intent(in) :: tau_climb
    real(k_real), intent(out) :: ctau, dctau_dtauclimb
  end subroutine

  module subroutine cellwall_local_concentration(this, rho_cw, c_CW_rho, dcCW_dtauclimb)
    implicit none
    class(wkkct_climb_base), intent(in) :: this
    real(k_real), intent(in) :: rho_cw
    real(k_real), intent(out) :: c_CW_rho, dcCW_dtauclimb
  end subroutine

  module subroutine hydrostatic_stress_equilibrium_concentration(this, csigma_h, dcsigma_h_dtauclimb)
    implicit none
    class(wkkct_climb_base), intent(in) :: this
    real(k_real), intent(out) :: csigma_h, dcsigma_h_dtauclimb
  end subroutine

end interface


!********************************************************************************!
! The climb model NOT considering irradiation described in Wen et al. 2020 EQ. 20
!********************************************************************************!
type, extends(wkkct_climb_base) :: wkkct_climb

contains
  procedure :: initParameters=>initParametersWKKCTClimb
  procedure :: initGridPointers => initGridPointersWKKCTClimb
  procedure :: setPointData => setPointDataWKKCTClimb
  procedure :: computeConcentrationInterstiatialBulk => computeConcentrationInterstiatialBulkWKKCTClimb
  procedure :: computeConcentrationVacancyBulk => computeConcentrationVacancyBulkWKKCTClimb
  procedure :: computeConcentrationVacancyCellWall => computeConcentrationVacancyCellWallWKKCTClimb
  procedure :: computeConcentrationVacancyCore => computeConcentrationVacancyCoreWKKCTClimb
  procedure :: readMaterialParametersFromFile => readMaterialParametersFromFileWKKCTClimb
  procedure :: computeZinterstitial => computeZinterstitialWKKCTClimb
  procedure :: computeZVacancy => computeZVacancyWKKCTClimb
  procedure :: acceptRejectSolution => acceptRejectSolutionWKKCTClimb
end type

interface
  module subroutine readMaterialParametersFromFileWKKCTClimb(this, matf_reader)
    use read_from_file_utils, only : file_reader
    implicit none
    class(wkkct_climb), intent(inout) :: this
    class(file_reader), intent(inout) :: matf_reader
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
  end subroutine

  module subroutine initGridPointersWKKCTClimb(this)
    implicit none
    class(wkkct_climb), intent(inout) :: this
  end subroutine

  module subroutine setPointDataWKKCTClimb(this, ix,iy,iz)
    implicit none
    class(wkkct_climb), intent(inout) :: this
    integer, intent(in) :: ix,iy,iz
  end subroutine

  module subroutine computeZinterstitialWKKCTClimb(this, tau_climb, zi, dzi_dtauclimb)
    implicit none
    class(wkkct_climb), intent(in) :: this
    real(k_real), intent(in) :: tau_climb
    real(k_real), intent(out) :: zi, dzi_dtauclimb
  end subroutine

  module subroutine computeZVacancyWKKCTClimb(this, tau_climb, zv, dzv_dtauclimb)
    implicit none
    class(wkkct_climb), intent(in) :: this
    real(k_real), intent(in) :: tau_climb
    real(k_real), intent(out) :: zv, dzv_dtauclimb
  end subroutine

  module subroutine computeConcentrationInterstiatialBulkWKKCTClimb(this, tau_climb, ci_bulk, dci_bulk_dtauclimb)
    implicit none
    class(wkkct_climb), intent(in) :: this
    real(k_real), intent(in) :: tau_climb
    real(k_real), intent(out) :: ci_bulk, dci_bulk_dtauclimb
  end subroutine

  module subroutine computeConcentrationVacancyBulkWKKCTClimb(this, tau_climb, cv_bulk, dcv_bulk_dtauclimb)
    implicit none
    class(wkkct_climb), intent(in) :: this
    real(k_real), intent(in) :: tau_climb
    real(k_real), intent(out) :: cv_bulk, dcv_bulk_dtauclimb
  end subroutine

  module subroutine computeConcentrationVacancyCellWallWKKCTClimb(this, rho_cw, cv_CW, dcv_CW_dtauclimb)
    implicit none
    class(wkkct_climb), intent(in) :: this
    real(k_real), intent(in) :: rho_cw
    real(k_real), intent(out) :: cv_CW, dcv_CW_dtauclimb
  end subroutine

  module subroutine computeConcentrationVacancyCoreWKKCTClimb(this, tau_climb, cv_core, dcv_core_dtauclimb)
    implicit none
    class(wkkct_climb), intent(in) :: this
    real(k_real), intent(in) :: tau_climb
    real(k_real), intent(out) :: cv_core, dcv_core_dtauclimb
  end subroutine

  module subroutine acceptRejectSolutionWKKCTClimb(this, dt_max, accept_solution_flag)
    implicit none
    class(wkkct_climb), intent(inout) :: this
    real(k_real), intent(out) :: dt_max
    logical, intent(out) :: accept_solution_flag
  end subroutine

end interface

!******************************************************************************!
! An extension of the climb model proposed by Wen et al. 2020 that includes
! irradiation damage
!******************************************************************************!
type, extends(wkkct_climb_base) :: wkkct_climb_irradiation
  real(k_real), pointer, dimension(:,:,:) :: conc_interstitial_grid => null(), &
                                             conc_vacancy_grid => null()
  real(k_real), pointer :: conc_interstitial_ptr => null(), &
                           conc_vacancy_ptr => null()

  ! parameters
  real(k_real) :: Zsi, pK0

contains
  procedure :: initParameters=>initParametersWKKCTClimbIrradiation
  procedure :: addFieldVariablesToGrid => addFieldVariablesToGridWKKCTClimbIrradiation
  procedure :: initGridPointers => initGridPointersWKKCTClimbIrradiation
  procedure :: setPointData => setPointDataWKKCTClimbIrradiation
  procedure :: initStateVariablesAtMaterialPoint => initStateVarsAtMaterialPointWKKCTClimbIrradiation
  procedure :: updateStateVariablesAtMaterialPointInnerLoop => updateStateVarsAtMaterialPointInnerLoopWKKCTClimbIrradiation
  procedure :: updateStateVariablesAtMaterialPointStaggered => updateStateVarsAtMaterialPointStaggeredWKKCTClimbIrradiation

  procedure :: computeSaturationConcentrationWKKCTClimbIrradiation
  procedure :: computeConcentrationInterstiatialBulk => computeConcentrationInterstiatialBulkWKKCTClimbIrradiation
  procedure :: computeConcentrationVacancyBulk => computeConcentrationVacancyBulkWKKCTClimbIrradiation
  procedure :: computeConcentrationVacancyCore => computeConcentrationVacancyCoreWKKCTClimbIrradiation
  procedure :: readMaterialParametersFromFile => readMaterialParametersFromFileWKKCTClimbIrradiation
  procedure :: computeZinterstitial => computeZinterstitialWKKCTClimbIrradiation
  procedure :: computeZVacancy => computeZVacancyWKKCTClimbIrradiation
end type

interface
  module subroutine readMaterialParametersFromFileWKKCTClimbIrradiation(this, matf_reader)
    use read_from_file_utils, only : file_reader
    implicit none
    class(wkkct_climb_irradiation), intent(inout) :: this
    class(file_reader), intent(inout) :: matf_reader
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
  end subroutine

  module subroutine  addFieldVariablesToGridWKKCTClimbIrradiation(this)
    implicit none
    class(wkkct_climb_irradiation), intent(inout) :: this
  end subroutine

  module subroutine initGridPointersWKKCTClimbIrradiation(this)
    implicit none
    class(wkkct_climb_irradiation), intent(inout) :: this
  end subroutine

  module subroutine setPointDataWKKCTClimbIrradiation(this, ix,iy,iz)
    implicit none
    class(wkkct_climb_irradiation), intent(inout) :: this
    integer, intent(in) :: ix,iy,iz
  end subroutine

  module subroutine initStateVarsAtMaterialPointWKKCTClimbIrradiation(this)
    implicit none
    class(wkkct_climb_irradiation), intent(inout) :: this
  end subroutine

  module subroutine updateStateVarsAtMaterialPointInnerLoopWKKCTClimbIrradiation(this)
    implicit none
    class(wkkct_climb_irradiation), intent(inout) :: this
  end subroutine

  module subroutine updateStateVarsAtMaterialPointStaggeredWKKCTClimbIrradiation(this)
    implicit none
    class(wkkct_climb_irradiation), intent(inout) :: this
  end subroutine

  module subroutine computeZinterstitialWKKCTClimbIrradiation(this, tau_climb, zi, dzi_dtauclimb)
    implicit none
    class(wkkct_climb_irradiation), intent(in) :: this
    real(k_real), intent(in) :: tau_climb
    real(k_real), intent(out) :: zi, dzi_dtauclimb
  end subroutine

  module subroutine computeZVacancyWKKCTClimbIrradiation(this, tau_climb, zv, dzv_dtauclimb)
    implicit none
    class(wkkct_climb_irradiation), intent(in) :: this
    real(k_real), intent(in) :: tau_climb
    real(k_real), intent(out) :: zv, dzv_dtauclimb
  end subroutine

  module subroutine computeConcentrationInterstiatialBulkWKKCTClimbIrradiation(this, tau_climb, ci_bulk, dci_bulk_dtauclimb)
    implicit none
    class(wkkct_climb_irradiation), intent(in) :: this
    real(k_real), intent(in) :: tau_climb
    real(k_real), intent(out) :: ci_bulk, dci_bulk_dtauclimb
  end subroutine

  module subroutine computeConcentrationVacancyBulkWKKCTClimbIrradiation(this, tau_climb, cv_bulk, dcv_bulk_dtauclimb)
    implicit none
    class(wkkct_climb_irradiation), intent(in) :: this
    real(k_real), intent(in) :: tau_climb
    real(k_real), intent(out) :: cv_bulk, dcv_bulk_dtauclimb
  end subroutine


  module subroutine computeConcentrationVacancyCoreWKKCTClimbIrradiation(this, tau_climb, cv_core, dcv_core_dtauclimb)
    implicit none
    class(wkkct_climb_irradiation), intent(in) :: this
    real(k_real), intent(in) :: tau_climb
    real(k_real), intent(out) :: cv_core, dcv_core_dtauclimb
  end subroutine

  module subroutine  computeSaturationConcentrationWKKCTClimbIrradiation(this)
    implicit none
    class(wkkct_climb_irradiation), intent(inout) :: this
  end subroutine

end interface
!
! !******************************************************************************!
! ! The chemo-mechanical climb model
! !******************************************************************************!
! type, extends(wkkct_climb_base) :: wkkct_climb_chemo_mech
!   real(k_real), pointer, dimension(:,:,:) :: concentration_grid => null()
!   real(k_real), pointer :: concentration_ptr => null()
!
!   ! parameters
!   real(k_real) :: z_vacancy0, diff_vacancy, conc_vacancy_0
!
! contains
!   procedure :: initParameters=>initParametersWKKCTClimbChemoMechanical
!   procedure :: computeClimbVelocity => computeClimbVelocityWKKCTClimbChemoMechanical
!   procedure :: initGridPointers => initGridPointersWKKCTClimbChemoMechanical
!   procedure :: setPointData => setPointDataWKKCTClimbChemoMechanical
! end type
!
! interface
!   module subroutine initParametersWKKCTClimbChemoMechanical(this, phase_id, n_ss_per_mode, use_damage, n_gauss, &
!                               burgVectorL, climb_reduce_ratio, ratio_edge, zv0, diff_v, c_v0)
!     implicit none
!     class(wkkct_climb_chemo_mech), intent(inout) :: this
!     integer, intent(in) :: phase_id, n_ss_per_mode(:)
!     real(k_real), target, intent(in), dimension(:) :: burgVectorL
!     logical, intent(in) :: use_damage
!     integer, intent(in):: n_gauss
!     real(k_real), target, intent(in) :: climb_reduce_ratio, ratio_edge, &           !base wkkct parameters
!                                         zv0, diff_v, c_v0                     !model specific parameters
!   end subroutine
!
!   module subroutine computeClimbVelocityWKKCTClimbChemoMechanical(this, rss)
!     implicit none
!     class(wkkct_climb_chemo_mech), intent(inout) :: this
!     real(k_real), intent(in) :: rss
!   end subroutine
!
!   module subroutine initGridPointersWKKCTClimbChemoMechanical(this)
!     implicit none
!     class(wkkct_climb_chemo_mech), intent(inout) :: this
!   end subroutine
!
!   module subroutine setPointDataWKKCTClimbChemoMechanical(this, ix,iy,iz)
!     implicit none
!     class(wkkct_climb_chemo_mech), intent(inout) :: this
!     integer, intent(in) :: ix,iy,iz
!   end subroutine
!
! end interface
contains
subroutine readMaterialParametersFromFileCPClimb(matf_reader, phase_id, &
                                                  use_damage, n_gauss, n_std_dev_gauss_integration, &
                                                  all_mighty_grid_in, sim_all_macro_data, the_bc_object, stiffness_ptr, &
                                                  common_material_parameter_ptr, crystal_paremeters_ptr, inelastic_strain_base_ptr)
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

!CLIMB MODELS
class(transient_climb), pointer :: transient_climb_temp  => null()
class(exponential_climb), pointer :: exponential_climb_temp  => null()
class(wkkct_climb), pointer :: wkkct_climb_temp  => null()
class(wkkct_climb_irradiation), pointer :: wkkct_climb_irradiation_temp  => null()
logical, intent(in) :: use_damage
integer, intent(in) :: n_gauss
real(k_real), intent(in) :: n_std_dev_gauss_integration


call matf_reader%readParameter("--Climb-model", material_model_name)
select case(material_model_name%getString())
case ("wkkct-climb")
  allocate(wkkct_climb_temp)
  inelastic_strain_base_ptr => wkkct_climb_temp
  call wkkct_climb_temp%linkBaseObjects(phase_id, all_mighty_grid_in, sim_all_macro_data, the_bc_object)
  call wkkct_climb_temp%initParameters(phase_id, common_material_parameter_ptr,&
                                      use_damage, n_gauss, n_std_dev_gauss_integration, stiffness_ptr, crystal_paremeters_ptr)
case ("wkkct-climb-irradiation")
  allocate(wkkct_climb_irradiation_temp)
  inelastic_strain_base_ptr => wkkct_climb_irradiation_temp
  call wkkct_climb_irradiation_temp%linkBaseObjects(phase_id, all_mighty_grid_in, sim_all_macro_data, the_bc_object)
  call wkkct_climb_irradiation_temp%initParameters(phase_id, common_material_parameter_ptr,&
                                                  use_damage, n_gauss, n_std_dev_gauss_integration, stiffness_ptr, crystal_paremeters_ptr)
case ("exponential-climb")
  allocate(exponential_climb_temp)
  inelastic_strain_base_ptr => exponential_climb_temp
  call exponential_climb_temp%linkBaseObjects(phase_id, all_mighty_grid_in, sim_all_macro_data, the_bc_object)
  call exponential_climb_temp%initParameters(phase_id, common_material_parameter_ptr,&
                                             use_damage, n_gauss, n_std_dev_gauss_integration, stiffness_ptr, crystal_paremeters_ptr)

case ("transient-climb")
  allocate(transient_climb_temp)
  inelastic_strain_base_ptr => transient_climb_temp
  call transient_climb_temp%linkBaseObjects(phase_id, all_mighty_grid_in, sim_all_macro_data, the_bc_object)
  call transient_climb_temp%initParameters(phase_id, common_material_parameter_ptr,&
                                        use_damage, n_gauss, n_std_dev_gauss_integration, stiffness_ptr, crystal_paremeters_ptr)
case default
  write(*,*) "unrecognized material model for climb", material_model_name%getString()
  write(*,*) "available material models are: "
  write(*,*) "transient-climb, exponential-climb, wkkct-climb, and  wkkct-climb-irradiation"
  
  error stop "abort"
end select

call inelastic_strain_base_ptr%readMaterialParametersFromFile(matf_reader)
select case(material_model_name%getString())
case ("wkkct-climb")
  nullify(wkkct_climb_temp)
case ("wkkct-climb-irradiation")
  nullify(wkkct_climb_irradiation_temp)
case ("exponential-climb")
  nullify(exponential_climb_temp)
case ("transient-climb")
  nullify(transient_climb_temp)
case default
  write(*,*) "somthing is worng when trying nullifying the pointer for ", material_model_name%getString()
  error stop "abort"
end select

end subroutine

end module
