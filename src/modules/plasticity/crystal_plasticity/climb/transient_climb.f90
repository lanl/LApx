#include "macro_debug.fpp"

submodule(climb_mod) transient_climb_smod
implicit none
contains
!
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

  call this%initParametersClimbBase(phase_id, common_material_parameter_ptr,&
                                    use_damage, n_gauss, n_std_dev_gauss_integration, &
                                    elasticity_obj, crystal_paremeters_ptr)

  this%elasticity_obj => elasticity_obj
end subroutine

module subroutine readMaterialParametersFromFileTransientClimb(this, matf_reader)
  use read_from_file_utils, only : file_reader
  implicit none
  class(transient_climb), intent(inout) :: this
  class(file_reader), intent(inout) :: matf_reader
  ! real(k_real), dimension(:), pointer :: mode_dep_parameter => null()

  ! call readVectorParameter(file_id, "beta-dot-0[unitless]", this%n_slip_modes, mode_dep_parameter)
  ! call this%convertModeToSlipSystemParameter(mode_dep_parameter, this%bdot0)
  ! deallocate(mode_dep_parameter); nullify(mode_dep_parameter)
  ! call readVectorParameter(file_id, "n-exponent[unitless]", this%n_slip_modes, mode_dep_parameter)
  ! call this%convertModeToSlipSystemParameter(mode_dep_parameter, this%n_exp)
  ! deallocate(mode_dep_parameter); nullify(mode_dep_parameter)
  call matf_reader%readParameter("n-exponent-damage[unitless]", this%n_exp_damage)
  ! call readVectorParameter(file_id, "initial-critical-resolved-shear-stress[MPa]", this%n_slip_modes, mode_dep_parameter)
  ! call this%convertModeToSlipSystemParameter(mode_dep_parameter, this%crss0)
  ! deallocate(mode_dep_parameter); nullify(mode_dep_parameter)
  ! call readParameter(file_id, "rss-std-deviation[MPa]", this%rss_std_dev)
  this%atomic_volume = this%common_material_parameter_ptr%atomic_volume
  this%vacancy_diffusivity_coefficient => this%common_material_parameter_ptr%vacancy_diffusivity_coefficient ! point to the bulk value
  this%E_formation_vacancy => this%common_material_parameter_ptr%vacancy_formation_energy ! point to the bulk value
  this%E_migration_vacancy => this%common_material_parameter_ptr%vacancy_migration_energy ! point to the bulk value
  this%E_migration_interstitial => this%common_material_parameter_ptr%vacancy_migration_energy ! point to the bulk value
  this%E_migration_interstitial_cluster => this%common_material_parameter_ptr%vacancy_migration_energy ! point to the bulk value
  this%interstitial_diffusivity_coefficient => this%common_material_parameter_ptr%vacancy_diffusivity_coefficient ! point to the bulk value
  this%interstitial_cluster_diffusivity_coefficient => this%common_material_parameter_ptr%vacancy_diffusivity_coefficient ! point to the bulk value
  this%voxel_volume = 1._k_real

end subroutine

module subroutine  addFieldVariablesToGridTransientClimb(this)
  use grid_data_var_type_mod
  implicit none
  class(transient_climb), intent(inout) :: this

  call this%cp_climb_base%addFieldVariablesToGrid()
  associate (all_grid_data_vars => this%grid_data)

  call all_grid_data_vars%addVar("rho_m", ss_scalar, stateful_level=2)
  call all_grid_data_vars%addVar("climb_velocity", ss_generic_vector, &
                                 additional_var_dimensions=(/this%n_gauss/), stateful_level=2)
  call all_grid_data_vars%addVar("climb_velocity_cellwall", ss_generic_vector, &
                                 additional_var_dimensions=(/this%n_gauss/), stateful_level=2)
  call all_grid_data_vars%addVar("vacancy_concentration", scalar, stateful_level=2)
  call all_grid_data_vars%addVar("interstitial_concentration", scalar, stateful_level=2)
  call all_grid_data_vars%addVar("interstitial_cluster_concentration", generic_vector, &
                                 additional_var_dimensions=(/this%n_cd/),stateful_level=2)
  call all_grid_data_vars%addVar("vacancy_loop_density", generic_vector, &
                                 additional_var_dimensions=(/this%n_cd/),stateful_level=2)
  call all_grid_data_vars%addVar("vacancy_loop_radii", generic_vector, &
                                 additional_var_dimensions=(/this%n_cd/),stateful_level=2)
  call all_grid_data_vars%addVar("interstitial_loop_density", generic_vector, &
                                 additional_var_dimensions=(/this%n_cd/),stateful_level=2)
  call all_grid_data_vars%addVar("interstitial_loop_radii", generic_vector, &
                                 additional_var_dimensions=(/this%n_cd/),stateful_level=2)

  call all_grid_data_vars%addVar("dislocation_climb_strain_rate", vector6)
  call all_grid_data_vars%addVar("dislocation_climb_strain", vector6,stateful_level=2)
  
  call all_grid_data_vars%addVar("vacancy_loop_strain_rate", vector6)
  call all_grid_data_vars%addVar("vacancy_loop_strain", vector6,stateful_level=2)

  call all_grid_data_vars%addVar("interstitial_loop_strain_rate", vector6)
  call all_grid_data_vars%addVar("interstitial_loop_strain", vector6,stateful_level=2)

  call all_grid_data_vars%addVar("loops_climb_strain_rate", vector6)
  call all_grid_data_vars%addVar("loops_climb_strain", vector6,stateful_level=2)

  call all_grid_data_vars%addVar("dose_a", scalar,stateful_level=2)
  call all_grid_data_vars%addVar("dose_c", scalar,stateful_level=2)

  end associate
end subroutine

module subroutine initGridPointersTransientClimb(this)
  implicit none
  class(transient_climb), intent(inout) :: this

  call this%cp_climb_base%initGridPointers()

  call this%grid_data%getSSScalarDataPointerByName("rho_m", this%rho_m_grid)
  call this%grid_data%getSSGenericVectorDataPointerByName("climb_velocity", this%climb_velocity_grid)
  call this%grid_data%getSSGenericVectorDataPointerByName("climb_velocity_CellWall", this%climb_velocity_CellWall_grid)
  call this%grid_data%getScalarDataPointerByName("vacancy_concentration", this%vacancy_conc_grid)
  call this%grid_data%getScalarDataPointerByName("interstitial_concentration", this%interstitial_conc_grid)
  call this%grid_data%getGenericVectorDataPointerByName("interstitial_cluster_concentration", this%interstitial_cluster_conc_grid)
  call this%grid_data%getScalarDataPointerByNameOld("vacancy_concentration", this%vacancy_conc_grid_old)
  call this%grid_data%getScalarDataPointerByNameOld("interstitial_concentration", this%interstitial_conc_grid_old)
  call this%grid_data%getGenericVectorDataPointerByNameOld("interstitial_cluster_concentration", this%interstitial_cluster_conc_grid_old)

  call this%grid_data%getGenericVectorDataPointerByName("vacancy_loop_density", this%vacancy_loop_density_grid)
  call this%grid_data%getGenericVectorDataPointerByName("vacancy_loop_radii", this%vacancy_loop_radii_grid)

  call this%grid_data%getGenericVectorDataPointerByName("interstitial_loop_density", this%interstitial_loop_density_grid)
  call this%grid_data%getGenericVectorDataPointerByName("interstitial_loop_radii", this%interstitial_loop_radii_grid)
  
  call this%grid_data%getVector6DataPointerByName("dislocation_climb_strain_rate",  this%dislocation_climb_strain_rate_grid)
  call this%grid_data%getVector6DataPointerByName("dislocation_climb_strain", this%dislocation_climb_strain_grid)
  call this%grid_data%getVector6DataPointerByNameOld("dislocation_climb_strain", this%dislocation_climb_strain_grid_old)
  
  call this%grid_data%getVector6DataPointerByName("vacancy_loop_strain_rate",  this%vacancy_loop_strain_rate_grid)
  call this%grid_data%getVector6DataPointerByName("vacancy_loop_strain", this%vacancy_loop_strain_grid)
  call this%grid_data%getVector6DataPointerByNameOld("vacancy_loop_strain", this%vacancy_loop_strain_grid_old)

  call this%grid_data%getVector6DataPointerByName("interstitial_loop_strain_rate",  this%interstitial_loop_strain_rate_grid)
  call this%grid_data%getVector6DataPointerByName("interstitial_loop_strain", this%interstitial_loop_strain_grid)
  call this%grid_data%getVector6DataPointerByNameOld("interstitial_loop_strain", this%interstitial_loop_strain_grid_old)

  call this%grid_data%getVector6DataPointerByName("loops_climb_strain_rate",  this%loops_climb_strain_rate_grid)
  call this%grid_data%getVector6DataPointerByName("loops_climb_strain", this%loops_climb_strain_grid)
  call this%grid_data%getVector6DataPointerByNameOld("loops_climb_strain", this%loops_climb_strain_grid_old)

  call this%grid_data%getScalarDataPointerByName("dose_a", this%dose_a_grid)
  call this%grid_data%getScalarDataPointerByName("dose_c", this%dose_c_grid)
  call this%grid_data%getScalarDataPointerByNameOld("dose_a", this%dose_a_grid_old)
  call this%grid_data%getScalarDataPointerByNameOld("dose_c", this%dose_c_grid_old)

end subroutine

module subroutine setPointDataTransientClimb(this, ix,iy,iz)
  implicit none
  class(transient_climb), intent(inout) :: this
  integer, intent(in) :: ix,iy,iz

  !calling the parent class to do its job
  call this%cp_climb_base%setPointData(ix,iy,iz)
  this%rho_m_ptr => this%rho_m_grid(:,ix,iy,iz)
  this%climb_velocity_ptr => this%climb_velocity_grid(:,:,ix,iy,iz)
  this%climb_velocity_CellWall_ptr => this%climb_velocity_CellWall_grid(:,:,ix,iy,iz)
  this%vacancy_conc_ptr => this%vacancy_conc_grid(ix,iy,iz)
  this%interstitial_conc_ptr => this%interstitial_conc_grid(ix,iy,iz)
  this%interstitial_cluster_conc_ptr => this%interstitial_cluster_conc_grid(:,ix,iy,iz)
  this%vacancy_conc_ptr_old => this%vacancy_conc_grid_old(ix,iy,iz)
  this%interstitial_conc_ptr_old => this%interstitial_conc_grid_old(ix,iy,iz)
  this%interstitial_cluster_conc_ptr_old => this%interstitial_cluster_conc_grid_old(:,ix,iy,iz)

  this%vacancy_loop_density_ptr => this%vacancy_loop_density_grid(:,ix,iy,iz)
  this%vacancy_loop_radii_ptr => this%vacancy_loop_radii_grid(:,ix,iy,iz)

  this%interstitial_loop_density_ptr => this%interstitial_loop_density_grid(:,ix,iy,iz)
  this%interstitial_loop_radii_ptr => this%interstitial_loop_radii_grid(:,ix,iy,iz)

  this%dislocation_climb_strain_rate_ptr => this%dislocation_climb_strain_rate_grid(:,ix,iy,iz)
  this%dislocation_climb_strain_ptr => this%dislocation_climb_strain_grid(:,ix,iy,iz)
  this%dislocation_climb_strain_ptr_old => this%dislocation_climb_strain_grid_old(:,ix,iy,iz)

  this%vacancy_loop_strain_rate_ptr => this%vacancy_loop_strain_rate_grid(:,ix,iy,iz)
  this%vacancy_loop_strain_ptr => this%vacancy_loop_strain_grid(:,ix,iy,iz)
  this%vacancy_loop_strain_ptr_old => this%vacancy_loop_strain_grid_old(:,ix,iy,iz)
  
  this%interstitial_loop_strain_rate_ptr => this%interstitial_loop_strain_rate_grid(:,ix,iy,iz)
  this%interstitial_loop_strain_ptr => this%interstitial_loop_strain_grid(:,ix,iy,iz)
  this%interstitial_loop_strain_ptr_old => this%interstitial_loop_strain_grid_old(:,ix,iy,iz) 

  this%loops_climb_strain_rate_ptr => this%loops_climb_strain_rate_grid(:,ix,iy,iz)
  this%loops_climb_strain_ptr => this%loops_climb_strain_grid(:,ix,iy,iz)
  this%loops_climb_strain_ptr_old => this%loops_climb_strain_grid_old(:,ix,iy,iz) 

  this%dose_a_ptr => this%dose_a_grid(ix,iy,iz) 
  this%dose_c_ptr => this%dose_c_grid(ix,iy,iz) 

  this%dose_a_ptr_old => this%dose_a_grid_old(ix,iy,iz) 
  this%dose_c_ptr_old => this%dose_c_grid_old(ix,iy,iz) 
end subroutine
!
module subroutine computeGammaDotAndDGammaDotDRssSSTransientClimb(this, rss, gdot, dgdot_drss)
  implicit none
  class(transient_climb), intent(inout) :: this
  real(k_real), intent(in) :: rss
  real(k_real), intent(out) :: gdot
  real(k_real), intent(out) :: dgdot_drss
  __DECL_CLASS_UNUSED_THIS__
  __DECL_UNUSED_REAL__



  __SUPPRESS_CLASS_UNUSED_THIS__
  __SUPPRESS_UNUSED_REAL__(rss)

  ! associate (bdot0 => this%bdot0, &
  !            crss => this%crss_ptr, &
  !            n => this%n_exp, &
  !            ss_idx => this%ss_idx )
  gdot = 0._k_real
  dgdot_drss = 0._k_real
  !   gdot = bdot0(ss_idx)*(abs(rss/crss(ss_idx)))**n(ss_idx) * &
  !                           sign(1._k_real,rss)

  !   ! we are assuming dcrss/drss = 0
  !   dgdot_drss = n(ss_idx) * bdot0(ss_idx)*&
  !                       (abs(rss/crss(ss_idx)))**(n(ss_idx)-1) / crss(ss_idx)

  ! end associate
end subroutine

module subroutine initStateVariablesAtMaterialPointTransientClimb(this)
  use test_utils_mod
  implicit none
  class(transient_climb), intent(inout) :: this

  this%vacancy_conc_ptr = 0.001_k_real
  this%interstitial_conc_ptr = 0.001_k_real
  this%interstitial_cluster_conc_ptr(1:3) = 0.0001_k_real
  this%interstitial_loop_radii_ptr = 0._k_real
  this%interstitial_loop_radii_ptr(1:3) = 1e-9_k_real*107._k_real
  this%interstitial_loop_density_ptr = 0._k_real
  this%interstitial_loop_density_ptr(1:3) = 1e18_k_real

  this%vacancy_loop_radii_ptr = 0._k_real
  this%vacancy_loop_radii_ptr(4) = 1e-9_k_real*107._k_real
  this%vacancy_loop_density_ptr = 0._k_real
  this%vacancy_loop_density_ptr(4) = 1e18_k_real*3._k_real*10._k_real
end subroutine

module subroutine updateStateVariablesAtMaterialPointOuterLoopTransientClimb(this)
  use math_constants, only : PI, kBSI
  use units_conversion_mod, only : MPa2Pa
  use tensor_math_mod, only : getIdentity
  use solve_linear_system, only : solve_Ax_eq_b
  implicit none
  class(transient_climb), intent(inout) :: this
  integer :: n_eq_v,   idx_v_start, idx_v_end,&
             n_eq_i,   idx_i_start, idx_i_end,&
             n_eq_ic,  idx_ic_start, idx_ic_end,&
             n_eq_rho, idx_rho_start, idx_rho_end,&
             n_eq_vl,  idx_vl_start, idx_vl_end,&
             n_eq_il,  idx_il_start, idx_il_end,&
             n_eq_gb,  idx_gb_start, idx_gb_end,&
             n_eq_total, n_conc_eq, &
             idx, iter, i, j
  real(k_real), pointer, dimension(:) :: Kvl_v => null() , &
                                         Kvl_i => null() , & 
                                         Kvl_ic => null() , &
                                         Kil_v => null(), &
                                         Kil_i => null(), &
                                         Kil_ic => null(), &
                                         Krho_v => null(), &
                                         Krho_i => null(), &
                                         Kgb_v => null(), &
                                         Kgb_i => null()
  real(k_real) :: rc_il, rc_vl, c0th_v, D_i, D_v, D_ic, f_ic, f_r, n_ic, gnrt
  real(k_real), dimension(this%n_ss) :: rss_without_damage
  real(k_real), dimension(6,this%n_ss) :: drss_without_damage_dstress

  call this%computeResolvedStressAndReslovedStressDStress(this%stress_ptr, with_damage_in=.FALSE.)
  rss_without_damage = this%rss
  drss_without_damage_dstress = this%drss_dstress

  call this%computeResolvedStressAndReslovedStressDStress(this%stress_ptr)

  call this%computeClimbDirectionStress(this%stress_ptr) 
 

  ! first thing we need to assemble the K^2 matrix which will have
  ! 1 row vacancy
  ! 1 row interstitial
  ! n_cd rows interstial clusters
  ! n_ss rows for dislocations
  ! n_cd rows for vacancy loops
  ! n_cd rows for interstitial loops

  f_ic=0.13_k_real
  f_r=0.97_k_real
  n_ic = 10._k_real
  gnrt = 1._k_real

  n_eq_total = 0
  n_eq_v=1;          idx_v_start=1;                         idx_v_end=idx_v_start+n_eq_v-1;
  n_eq_i=1;          idx_i_start=idx_v_start+n_eq_v;        idx_i_end=idx_i_start+n_eq_i-1;
  n_eq_ic=this%n_cd; idx_ic_start=idx_i_start+n_eq_i;       idx_ic_end=idx_ic_start+n_eq_ic-1;
  n_eq_rho=this%n_ss;idx_rho_start=idx_ic_start+n_eq_ic;    idx_rho_end=idx_rho_start+n_eq_rho-1;
  n_eq_vl=this%n_cd; idx_vl_start=idx_rho_start+n_eq_rho;   idx_vl_end=idx_vl_start+n_eq_vl-1;
  n_eq_il=this%n_cd; idx_il_start=idx_vl_start+n_eq_vl;     idx_il_end=idx_il_start+n_eq_il-1;
  n_eq_GB=this%n_cd; idx_gb_start=idx_il_start+n_eq_il;     idx_gb_end=idx_gb_start+n_eq_gb-1;


  ! write(*,*) "n_eq_v, idx_v_start, idx_v_end ", n_eq_v, idx_v_start, idx_v_end
  ! write(*,*) "n_eq_i, idx_i_start, idx_i_end ", n_eq_i, idx_i_start, idx_i_end
  ! write(*,*) "n_eq_ic, idx_ic_start, idx_ic_end ", n_eq_ic, idx_ic_start, idx_ic_end
  ! write(*,*) "n_eq_rho, idx_rho_start, idx_rho_end ", n_eq_rho, idx_rho_start, idx_rho_end
  ! write(*,*) "n_eq_vl, idx_vl_start, idx_vl_end ", n_eq_vl, idx_vl_start, idx_vl_end
  ! write(*,*) "n_eq_il, idx_il_start, idx_il_end ", n_eq_il, idx_il_start, idx_il_end
  ! write(*,*) "n_eq_GB, idx_GB_start, idx_GB_end ", n_eq_GB, idx_GB_start, idx_GB_end

  n_eq_total = n_eq_v+n_eq_i+n_eq_ic+n_eq_rho+n_eq_vl+n_eq_il+n_eq_GB
  n_conc_eq = n_eq_v+n_eq_i+n_eq_ic
  if (.not.(associated(this%K2))) allocate(this%K2(n_eq_total,n_eq_total))
  if (.not.(associated(this%c0))) allocate(this%c0(n_eq_total, n_eq_total))
  if (.not.(associated(this%G))) allocate(this%G(n_eq_total))
  if (.not.(associated(this%D))) allocate(this%D(n_eq_total))
  if (.not.(associated(this%cdot))) allocate(this%cdot(n_conc_eq))
  if (.not.(associated(this%dcdot_dc))) allocate(this%dcdot_dc(n_conc_eq, n_conc_eq))
  if (.not.(associated(this%I_conc))) allocate(this%I_conc(n_conc_eq, n_conc_eq))
  call getIdentity(this%I_conc)
  if (.not.(associated(this%Residual))) allocate(this%Residual(n_conc_eq))
  if (.not.(associated(this%dResidualdC))) allocate(this%dResidualdC(n_conc_eq, n_conc_eq))
  if (.not.(associated(this%c_guess))) allocate(this%c_guess(n_conc_eq))
  if (.not.(associated(this%delta_c_guess))) allocate(this%delta_c_guess(n_conc_eq))
  if (.not.(associated(this%c))) allocate(this%c(n_conc_eq))
  if (.not.(associated(this%c_old))) allocate(this%c_old(n_conc_eq))
  if (.not.(associated(this%n_atoms_per_defect))) allocate(this%n_atoms_per_defect(n_eq_total))
  if (.not.(associated(this%n_stored_atoms))) allocate(this%n_stored_atoms(n_eq_total))
  if (.not.(associated(this%n0_stored_atoms))) allocate(this%n0_stored_atoms(n_eq_total, n_eq_total))
  if (.not.(associated(this%n_dot))) allocate(this%n_dot(n_eq_total))

  if (.not.(associated(Kvl_v))) allocate(Kvl_v(n_eq_vl))
  if (.not.(associated(Kvl_i))) allocate(Kvl_i(n_eq_vl))
  if (.not.(associated(Kvl_ic))) allocate(Kvl_ic(n_eq_vl))
 
  if (.not.(associated(Kil_v))) allocate(Kil_v(n_eq_il))
  if (.not.(associated(Kil_i))) allocate(Kil_i(n_eq_il))
  if (.not.(associated(Kil_ic))) allocate(Kil_ic(n_eq_il))

  if (.not.(associated(Krho_v))) allocate(Krho_v(n_eq_rho))
  if (.not.(associated(Krho_i))) allocate(Krho_i(n_eq_rho))

  if (.not.(associated(Kgb_v))) allocate(Kgb_v(n_eq_gb))
  if (.not.(associated(Kgb_i))) allocate(Kgb_i(n_eq_gb))
  
  
  associate(Nvl => this%vacancy_loop_density_ptr, &
            rvl => this%vacancy_loop_radii_ptr, &
            Nil => this%interstitial_loop_density_ptr, &
            ril => this%interstitial_loop_radii_ptr, &
            rho => this%rho_m_ptr, &
            Diff_v => this%vacancy_diffusivity_coefficient, &
            Diff_i => this%interstitial_diffusivity_coefficient, &
            Diff_ic => this%interstitial_cluster_diffusivity_coefficient, &
            T => this%temperature, &
            dt => this%dt, &
            Em_v=> this%E_migration_vacancy, &
            Ef_v=> this%E_formation_vacancy, &
            Em_i=> this%E_migration_interstitial, &
            Em_ic=> this%E_migration_interstitial_cluster, &
            omega => this%atomic_volume, &
            K2 => this%K2, &
            D => this%D, &
            C0 => this%C0, &
            G => this%G, &
            cdot => this%cdot, &
            dcdot_dc => this%dcdot_dc, &
            R => this%Residual, &
            dR_dC => this%dResidualdC, &
            I_conc => this%I_conc, &
            c_guess => this%c_guess, &
            delta_c_guess => this%delta_c_guess, &
            c => this%c, &
            c_old=> this%c_old, &
            V => this%voxel_volume, &
            n_atoms_per_defect => this%n_atoms_per_defect, &
            n_stored_atoms => this%n_stored_atoms, &
            n0_stored_atoms => this%n0_stored_atoms, &
            n_dot => this%n_dot &
            )
  

  ! Nvl = 1._k_real
  ! rvl = 1._k_real
  ! Nil = 1._k_real
  ! ril = 1._k_real
  rc_vl = 1e-5_k_real
  rc_il = 1e-5_k_real

  Kvl_v = 2._k_real*PI*rvl*Nvl
  Kvl_i = 2._k_real*PI*rvl*Nvl
  Kvl_ic =  (PI*rc_vl)**2 * (PI*rvl*Nvl) * 2._k_real * PI*rvl*Nvl

  Kil_v = 2._k_real*PI*ril*Nil
  Kil_i = 2._k_real*PI*ril*Nil
  Kil_ic = (PI*rc_il)**2 * (PI*ril*Nil) * 2._k_real * PI*ril*Nil

  Krho_v = rho
  Krho_i = rho

  K2=0._k_real
  ! filling the matrix
  do idx=idx_v_start, idx_v_end
    K2(idx, idx_vl_start:idx_vl_end) = Kvl_v 
    K2(idx, idx_il_start:idx_il_end)  = Kil_v 
    K2(idx, idx_rho_start:idx_rho_end)   = Krho_v 
    ! K2(idx, idx_gb_start:idx_gb_end)   = Kgb_v 
  enddo
  do idx=idx_i_start, idx_i_end
    K2(idx, idx_vl_start:idx_vl_end) = Kvl_i
    K2(idx, idx_il_start:idx_il_end)   = Kil_i
    K2(idx, idx_rho_start:idx_rho_end)   = Krho_i
    ! K2(idx, idx_gb_start:idx_gb_end)   = Kgb_i
  enddo

  ! this are only diagonal (i.e. ic interact only with il and vl with the same burger vector)
  i = 0
  do idx=idx_ic_start, idx_ic_end-1
    i = i + 1
    K2(idx, idx_il_start+i-1) = Kil_ic(i)/100._k_real
    K2(idx, idx_vl_start+i-1) = Kvl_ic(i)/100._k_real
  enddo
  
  ! if (this%ix.eq.1.and.this%iy.eq.1.and.this%iz.eq.1) then
  ! do idx=1,n_eq_total
  !   write(*,*) "idx K2 ", idx, K2(idx,:)
  ! end do
  ! endif
  ! fill the equilibrium concentration matrix (c0=cth*csigma)
  c0th_v = exp(-Ef_v(this%gb_id)/(kBSI*T))
  C0 = 0._k_real
  do idx=idx_v_start, idx_v_end
    C0(idx, idx_vl_start:idx_vl_end) = c0th_v*exp(this%cds*MPa2Pa*omega/(kBSI*T)) 
    C0(idx, idx_il_start:idx_il_end) = c0th_v*exp(this%cds*MPa2Pa*omega/(kBSI*T)) 
    C0(idx, idx_rho_start:idx_rho_end)   = c0th_v*exp(this%rss*MPa2Pa*omega/(kBSI*T)) 
    ! write(*,*) "C0(idx, idx_vl_start:idx_vl_end) ", exp(this%cds*MPa2Pa*omega/(kBSI*T)) 
    ! write(*,*) "C0(idx, idx_il_start:idx_il_end)  ", exp(this%cds*MPa2Pa*omega/(kBSI*T)) 
    ! write(*,*) "C0(idx, idx_rho_start:idx_rho_end) ", exp(this%rss*MPa2Pa*omega/(kBSI*T)) 
    ! write(*,*) "this%cds", this%cds
    ! write(*,*) "this%rss", this%rss
  enddo
  ! if (this%ix.eq.1.and.this%iy.eq.1.and.this%iz.eq.1) then
  !   do idx=1,n_eq_total
  !     write(*,*) "idx C0 ", idx, C0(idx,:)
  !   end do
  !   endif

  ! C0 = 0._k_real
  ! fill the diffusivity matrix (D=d0*exp(-E_migration/(KBSI*T)))
  ! write(*,*) "Diff_v ", Diff_v
  ! write(*,*) "Em_v ", Em_v
  ! write(*,*) "this%gb_id ", this%gb_id
  D_v = Diff_v(this%gb_id) * exp(-Em_v(this%gb_id)/(kBSI*T))
  D = 0._k_real
  do idx=idx_v_start, idx_v_end
    D(idx) = D_v
  enddo

  D_i = Diff_i(this%gb_id) * exp(-Em_i(this%gb_id)/(kBSI*T))
  do idx=idx_i_start, idx_i_end
    D(idx) = D_i
  enddo

  D_ic = Diff_ic(this%gb_id) * exp(-Em_ic(this%gb_id)/(kBSI*T))
  do idx=idx_ic_start, idx_ic_end-1 ! no cluster in the c direction
    D(idx) = D_ic
  enddo
  D(idx_ic_end) = 0._k_real 
  
  G =0._k_real
  do idx=idx_v_start, idx_v_end
    G(idx) = gnrt*(1._k_real-f_r)
  enddo

  do idx=idx_i_start, idx_i_end
    G(idx) = gnrt*(1._k_real-f_r)*(1._k_real-f_ic)
  enddo

  do idx=idx_ic_start, idx_ic_end-1
    G(idx) = gnrt*(1._k_real-f_r)*(f_ic)/(3._k_real*n_ic)
  enddo


  c(idx_v_start) = this%vacancy_conc_ptr_old
  c(idx_i_start) = this%interstitial_conc_ptr_old
  c(idx_ic_start:idx_ic_end) = this%interstitial_cluster_conc_ptr_old
  c_old(idx_v_start) = this%vacancy_conc_ptr_old
  c_old(idx_i_start) = this%interstitial_conc_ptr_old
  c_old(idx_ic_start:idx_ic_end) = this%interstitial_cluster_conc_ptr_old

  c_guess = c
  ! build cdot 
  cdot = 0._k_real
  dcdot_dc = 0._k_real
  do idx=1,n_conc_eq
    cdot(idx) = G(idx) - sum(K2(idx,:)*(c_guess(idx)-C0(idx,:)))*D(idx)
    dcdot_dc(idx,idx) = -sum(K2(idx,:))*D(idx)
  enddo

  ! get the new concentration
  iter = 0
  ! write(*,*) "shape(R) ", shape(R)
  ! write(*,*) "shape(dR_dc) ", shape(dR_dc)
  ! write(*,*) "shape(cdot) ", shape(cdot)
  ! write(*,*) "shape(dcdot_dc) ", shape(dcdot_dc)
  ! write(*,*) "shape(K2) ", shape(K2)
  ! write(*,*) "shape(c_guess) ", shape(c_guess)
  ! write(*,*) "shape(C0) ", shape(C0)
  ! write(*,*) "shape(D) ", shape(D)
  ! write(*,*) "shape(delta_c_guess) ", shape(delta_c_guess)
  ! write(*,*) "cdot ", cdot
  ! write(*,*) "c_old ", c_old
  ! write(*,*) "c_guess ", c_guess
  ! write(*,*) "this%rss ", this%rss
  ! write(*,*) "this%cds ", this%cds

  R = c_guess - (cdot*dt+c_old)
 
  do while (any(abs(R)>1e-10_k_real))
    ! write(*,*) "iter ", iter
    ! write(*,*) "dR_dc ", dR_dc
    dR_dc = I_conc - dcdot_dc*dt
    iter = iter+ 1
    call solve_Ax_eq_b(dR_dc, R, delta_c_guess)
    c_guess = c_guess - delta_c_guess
    do idx=1,n_conc_eq
      cdot(idx) = G(idx) - sum(K2(idx,:)*(c_guess(idx)-C0(idx,:)))*D(idx)
      dcdot_dc(idx,idx) = -sum(K2(idx,:))*D(idx)
    enddo
    R = c_guess - (cdot*dt+c_old)
    ! write(*,*) "R iter", iter, R
  enddo
  if (iter > 3 ) then 
    write(*,*) "iter", iter
    write(*,*) "R ", R
    error stop "this should converge in one iteration"
  endif
  this%c =  c_guess
  this%vacancy_conc_ptr = this%c(idx_v_start)
  this%interstitial_conc_ptr = this%c(idx_i_start)
  this%interstitial_cluster_conc_ptr = this%c(idx_ic_start:idx_ic_end)


  ! compute the number of defects
  n_atoms_per_defect = 0._k_real 
  n_atoms_per_defect(idx_v_start:idx_v_end) = -1._k_real
  n_atoms_per_defect(idx_i_start:idx_i_end) = 1._k_real
  n_atoms_per_defect(idx_ic_start:idx_ic_end) = n_ic

  n_stored_atoms = 0._k_real
  n0_stored_atoms = 0._k_real
  n_stored_atoms(1:n_conc_eq) = this%c * V  * n_atoms_per_defect(1:n_conc_eq)
  do idx = 1, n_eq_total
    n0_stored_atoms(idx,:) = c0(idx,:) * V  * n_atoms_per_defect(idx)
  enddo

  do j = 1, n_eq_total
     n_dot(j) = 0._k_real
  do i = 1, n_eq_total
      n_dot(j) =n_dot(j) + K2(i,j)*(n_stored_atoms(i)-n0_stored_atoms(i,j))*D(i) &
                  -  K2(j,i)*(n_stored_atoms(j)-n0_stored_atoms(j,i))*D(j)
  enddo
  enddo
  
  this%inelastic_strain_rate_ptr= 0._k_real
  this%dislocation_climb_strain_rate_ptr = 0._k_real
  do idx=1,this%n_ss
    this%inelastic_strain_rate_ptr = this%inelastic_strain_rate_ptr + &
                                     this%schmid_ptr(:,idx) * n_dot(idx_rho_start+idx-1)
    this%dislocation_climb_strain_rate_ptr = this%dislocation_climb_strain_rate_ptr + &
                                     this%schmid_ptr(:,idx) * n_dot(idx_rho_start+idx-1)                            
  enddo 
  this%dislocation_climb_strain_ptr = this%dislocation_climb_strain_ptr_old + &
                                      this%dislocation_climb_strain_rate_ptr*this%dt

  this%vacancy_loop_strain_rate_ptr= 0._k_real     
  this%interstitial_loop_strain_rate_ptr= 0._k_real     
  this%loops_climb_strain_rate_ptr= 0._k_real                                 
  do idx=1,this%n_cd
    this%inelastic_strain_rate_ptr = this%inelastic_strain_rate_ptr + &
                                     this%climb_direction_tensor_ptr(:,idx) * n_dot(idx_vl_start+idx-1)
    this%inelastic_strain_rate_ptr = this%inelastic_strain_rate_ptr + &
                                     this%climb_direction_tensor_ptr(:,idx) * n_dot(idx_il_start+idx-1)
    this%inelastic_strain_rate_ptr = this%inelastic_strain_rate_ptr + &
                                     this%climb_direction_tensor_ptr(:,idx) * n_dot(idx_gb_start+idx-1)
    
    this%vacancy_loop_strain_rate_ptr = this%vacancy_loop_strain_rate_ptr + &
                                     this%climb_direction_tensor_ptr(:,idx) * n_dot(idx_vl_start+idx-1)
    this%interstitial_loop_strain_rate_ptr = this%interstitial_loop_strain_rate_ptr + &
                                     this%climb_direction_tensor_ptr(:,idx) * n_dot(idx_il_start+idx-1)   
  enddo 
  this%loops_climb_strain_rate_ptr = this%vacancy_loop_strain_rate_ptr + this%interstitial_loop_strain_rate_ptr

  this%vacancy_loop_strain_ptr = this%vacancy_loop_strain_ptr_old + &
                                      this%vacancy_loop_strain_rate_ptr*this%dt
  this%interstitial_loop_strain_ptr = this%interstitial_loop_strain_ptr_old + &
                                      this%interstitial_loop_strain_rate_ptr*this%dt
  this%loops_climb_strain_ptr = this%loops_climb_strain_ptr_old + &
                                      this%loops_climb_strain_rate_ptr*this%dt
  end associate

  if (this%ix.eq.1.and.this%iy.eq.1.and.this%iz.eq.1) then
  write(*,*) "c ", this%c
  write(*,*) "c_old ", this%c_old
  write(*,*) "c_dot ", (this%c-this%c_old)/this%dt
  write(*,*) "n_dot(V) ", this%n_dot(idx_V_start:idx_V_end)
  write(*,*) "sum(n_dot(V)) ", sum(this%n_dot(idx_V_start:idx_V_end))
  write(*,*) "n_dot(I) ", this%n_dot(idx_I_start:idx_I_end)
  write(*,*) "sum(n_dot(I)) ", sum(this%n_dot(idx_I_start:idx_I_end))
  write(*,*) "n_dot(IC) ", this%n_dot(idx_IC_start:idx_IC_end)
  write(*,*) "sum(n_dot(IC)) ", sum(this%n_dot(idx_IC_start:idx_IC_end))
  write(*,*) "n_dot(rho) ", this%n_dot(idx_rho_start:idx_rho_end)
  write(*,*) "sum(n_dot(rho)) ", sum(this%n_dot(idx_rho_start:idx_rho_end))
  write(*,*) "n_dot(IL) ", this%n_dot(idx_il_start:idx_il_end)
  write(*,*) "sum(n_dot(IL)) ", sum(this%n_dot(idx_il_start:idx_il_end))
  write(*,*) "n_dot(VL) ", this%n_dot(idx_vl_start:idx_vl_end)
  write(*,*) "sum(n_dot(VL)) ", sum(this%n_dot(idx_vl_start:idx_vl_end))
  write(*,*) "n_dot(GB) ", this%n_dot(idx_gb_start:idx_gb_end)
  write(*,*) "sum(n_dot(GB)) ", sum(this%n_dot(idx_gb_start:idx_gb_end))
  write(*,*) "sum(n_dot) ", sum(this%n_dot)*this%atomic_volume
  endif
  ! clean up
  deallocate(this%D); nullify(this%D)
  deallocate(this%C0); nullify(this%C0)
  deallocate(this%K2); nullify(this%K2)
  deallocate(this%G); nullify(this%G)
  
  deallocate(Kvl_v); nullify(Kvl_v)
  deallocate(Kvl_i); nullify(Kvl_i)
  deallocate(Kvl_ic); nullify(Kvl_ic)

  deallocate(Kil_v); nullify(Kil_v)
  deallocate(Kil_i); nullify(Kil_i)
  deallocate(Kil_ic); nullify(Kil_ic)
 
  deallocate(Krho_v); nullify(Krho_v)
  deallocate(Krho_i); nullify(Krho_i)

  deallocate(Kgb_v); nullify(Kgb_v)
  deallocate(Kgb_i); nullify(Kgb_i)


  

end subroutine

module subroutine updateStateVariablesAtMaterialPointStaggeredTransientClimb(this)
  implicit none
  class(transient_climb), intent(inout) :: this
  __DECL_CLASS_UNUSED_THIS__
  ! real(k_real) :: fluence, dose, dose_rate, &
  !                 drho_ddose_a, drho_ddose_c, &
  !                 dddose_a_dt, dddose_c_dt, &
  !                 drhoa_dt, drhoc_dt, &
  !                 dNc_ddose(4), dNa_ddose(4), Ncmax


  __SUPPRESS_CLASS_UNUSED_THIS__
  ! fluence = 1._k_real
  ! dose = 15._k_real/7.2e25_k_real*fluence
  ! dose_rate = dose/this%dt
  ! this%dose_a_ptr = min(this%dose_a_ptr_old + dose_rate*this%dt, 20.83_k_real)
  ! this%dose_c_ptr = min(this%dose_c_ptr_old + dose_rate*this%dt, 20.83_k_real)

  ! dddose_a_dt = (this%dose_a_ptr - this%dose_a_ptr_old)/this%dt
  ! dddose_c_dt = (this%dose_c_ptr - this%dose_c_ptr_old)/this%dt

  ! drho_ddose_a = 0._k_real
  ! if (this%dose_a_ptr_old.lt.20.83_k_real) &
  !   drho_ddose_a = 1.98e11_k_real*dose**2 - 6.23e12*dose + 4.76e13_k_real

  ! drho_ddose_c = 0._k_real
  ! if (this%dose_c_ptr_old.lt.20.83_k_real) &
  !   drho_ddose_c = 1.87e11_k_real*dose**2 - 4.58e12*dose + 3.26e13_k_real

  ! drhoa_dt = drho_ddose_a*dddose_a_dt
  ! drhoc_dt = drho_ddose_c*dddose_c_dt

  ! if (this%vacancy_loop_density_ptr_old(4).lt.NVmax) then
  !   dNc_ddose = Ncmax*A*exp( A*(this%dose_c_ptr - dose_critialc_c) )
  ! endif
  
end subroutine

module subroutine acceptRejectSolutionTransientClimb(this, dt_max, accept_solution_flag)
  use mpi_useful_routines_mod, only : MPIMaxRelativeIncrementGridSSScalar, MPIMaxIncrementGridSSScalar
  use mpi_variables_mod, only : i_am_mpi_master
  use time_march_mod, only : computeMaxAllowedTimeStepLinear
  implicit none
  class(transient_climb), intent(inout) :: this
  real(k_real), intent(out) :: dt_max
  logical, intent(out) :: accept_solution_flag

  if (i_am_mpi_master) write(*,*) "*********************************"
  if (i_am_mpi_master) write(*,*) "acceptRejectSolutionTransientClimb"

  accept_solution_flag = .true.
  dt_max = this%dt_max

  if (i_am_mpi_master) write(*,*) "nothing to check"

  if (accept_solution_flag.and.i_am_mpi_master) write(*,*) "solution might be ACCEPTED"
  if (.not.accept_solution_flag.and.i_am_mpi_master) write(*,*) "solution will be REJECTED"
  if (i_am_mpi_master) write(*,*) "new allowable dt is ", dt_max
  if (i_am_mpi_master) write(*,*) "*********************************"
  if (i_am_mpi_master) write(*,*) ""

end subroutine

module subroutine getStrainRateaAndStressJacobianTransientClimb(this, stress6, epsilon_dot, depsilon_dot_dstress, ix, iy, iz)
  implicit none
  class(transient_climb), intent(inout) :: this
  real(k_real), target, intent(in) ::stress6(6)
  real(k_real), dimension(6), intent(out) :: epsilon_dot
  real(k_real), dimension(6,6), intent(out) :: depsilon_dot_dstress
  integer, intent(in) ::  ix, iy, iz
  __DECL_UNUSED_VECTOR_PTR__

  __SUPPRESS_UNUSED_VECTOR_WARNING__(stress6)
  ! first we set all grid pointers to the local material point
  call this%setPointData(ix, iy, iz)
  epsilon_dot = this%inelastic_strain_rate_ptr 
  depsilon_dot_dstress = 0._k_real
  ! call this%computeEpsilonDotAndDepsilonDotDStress(stress6, epsilon_dot, depsilon_dot_dstress)

  this%inelastic_strain_rate_ptr(:) = epsilon_dot
end subroutine

module subroutine writeInelastiStrainRateToFileTransientClimb(this, csv_writer_obj, write_headers)
  use write_to_file_utils_mod  
  use csv_writer_mod, only : csv_writer
  use mpi_useful_routines_mod, only : MPIAverageGridVectorMasked, MPIAverageGridScalarMasked, MPIAverageGridVectorMasked
  use change_tensor_basis, only : chg_basis_vector6_to_tensor2

  implicit none 
  class(transient_climb), intent(inout) :: this
  class(csv_writer), intent(inout) :: csv_writer_obj
  ! real(k_real) :: strain_rate_tensor(3,3), strain_rate_vector(6), scalar_avg, vector4_avg(4), !& 
  
  logical, intent(in) :: write_headers


  call this%cp_climb_base%writeInelastiStrainRateToFile( csv_writer_obj, write_headers)


  ! if (present(add_trailing_coma)) then
  !   add_trailing_coma_ = add_trailing_coma
  ! else
  !   add_trailing_coma_=.TRUE.
  ! end if

  ! ! call this%isInitialized()


  ! if(write_headers) then
  !   call csv_writer_obj%AppendTensorHeader("dislocation_climb_strain_rate")
  !   call csv_writer_obj%AppendTensorHeader("loops_climb_strain_rate")
  !   call csv_writer_obj%AppendTensorHeader("vacancy_loop_strain_rate")
  !   call csv_writer_obj%AppendTensorHeader("interstitial_loop_strain_rate")
  !   call csv_writer_obj%AppendTensorHeader("dislocation_climb_strain")
  !   call csv_writer_obj%AppendTensorHeader("loops_climb_strain")
  !   call csv_writer_obj%AppendTensorHeader("vacancy_loop_strain")
  !   call csv_writer_obj%AppendTensorHeader("interstitial_loop_strain")
  !   call csv_writer_obj%AppendScalarHeader("vacancy_concentration")
  !   call csv_writer_obj%AppendScalarHeader("interstitial_concentration")
  !   call csv_writer_obj%AppendVectorHeader(4, "interstitial_cluster_concentration")
  ! else 

  ! endif

  !   call MPIAverageGridVectorMasked(this%dislocation_climb_strain_rate_grid, this%phase_grid==this%phase_id, strain_rate_vector)
  !   call chg_basis_vector6_to_tensor2(strain_rate_vector,  strain_rate_tensor)
  !   if (i_am_mpi_master) call write2DArrayToFileInline(file_id, strain_rate_tensor, "strain", add_trailing_coma=.TRUE.)

  !   call MPIAverageGridVectorMasked(this%loops_climb_strain_rate_grid, this%phase_grid==this%phase_id, strain_rate_vector)
  !   call chg_basis_vector6_to_tensor2(strain_rate_vector,  strain_rate_tensor)
  !   if (i_am_mpi_master) call write2DArrayToFileInline(file_id, strain_rate_tensor, "strain", add_trailing_coma=.TRUE.)

  !   call MPIAverageGridVectorMasked(this%vacancy_loop_strain_rate_grid, this%phase_grid==this%phase_id, strain_rate_vector)
  !   call chg_basis_vector6_to_tensor2(strain_rate_vector,  strain_rate_tensor)
  !   if (i_am_mpi_master) call write2DArrayToFileInline(file_id, strain_rate_tensor, "strain", add_trailing_coma=.TRUE.)

  !   call MPIAverageGridVectorMasked(this%interstitial_loop_strain_rate_grid, this%phase_grid==this%phase_id, strain_rate_vector)
  !   call chg_basis_vector6_to_tensor2(strain_rate_vector,  strain_rate_tensor)
  !   if (i_am_mpi_master) call write2DArrayToFileInline(file_id, strain_rate_tensor, "strain", add_trailing_coma=.TRUE.)

  !   call MPIAverageGridVectorMasked(this%dislocation_climb_strain_grid, this%phase_grid==this%phase_id, strain_rate_vector)
  !   call chg_basis_vector6_to_tensor2(strain_rate_vector,  strain_rate_tensor)
  !   if (i_am_mpi_master) call write2DArrayToFileInline(file_id, strain_rate_tensor, "strain", add_trailing_coma=.TRUE.)

  !   call MPIAverageGridVectorMasked(this%loops_climb_strain_grid, this%phase_grid==this%phase_id, strain_rate_vector)
  !   call chg_basis_vector6_to_tensor2(strain_rate_vector,  strain_rate_tensor)
  !   if (i_am_mpi_master) call write2DArrayToFileInline(file_id, strain_rate_tensor, "strain", add_trailing_coma=.TRUE.)

  !   call MPIAverageGridVectorMasked(this%vacancy_loop_strain_grid, this%phase_grid==this%phase_id, strain_rate_vector)
  !   call chg_basis_vector6_to_tensor2(strain_rate_vector,  strain_rate_tensor)
  !   if (i_am_mpi_master) call write2DArrayToFileInline(file_id, strain_rate_tensor, "strain", add_trailing_coma=.TRUE.)

  !   call MPIAverageGridVectorMasked(this%interstitial_loop_strain_grid, this%phase_grid==this%phase_id, strain_rate_vector)
  !   call chg_basis_vector6_to_tensor2(strain_rate_vector,  strain_rate_tensor)
  !   if (i_am_mpi_master) call write2DArrayToFileInline(file_id, strain_rate_tensor, "strain", add_trailing_coma=.TRUE.)

  !   call MPIAverageGridScalarMasked(this%vacancy_conc_grid, this%phase_grid==this%phase_id, scalar_avg)
  !   if (i_am_mpi_master) call writeScalarToFileInline(file_id, scalar_avg, add_trailing_coma=.TRUE.)

  !   call MPIAverageGridScalarMasked(this%interstitial_conc_grid, this%phase_grid==this%phase_id, scalar_avg)
  !   if (i_am_mpi_master) call writeScalarToFileInline(file_id, scalar_avg, add_trailing_coma=.TRUE.)

  !   ! call MPIAverageGridVectorMasked(this%interstitial_cluster_conc_grid, this%phase_grid==this%phase_id, vector4_avg)
  !   ! if (i_am_mpi_master) call writeVectorToFileInline(file_id, scalar_avg, add_trailing_coma=add_trailing_coma_)
  ! endif
  
end subroutine

end submodule
