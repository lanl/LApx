MODULE GLOBAL
use kinds
use all_grid_data_mod
use fft_grid_data_mod, only : fft_grid_data_type
use simumaltion_macro_data_mod, only : sim_all_macro_data_obj
use gauss_legendre_integration_mod
use glide_mod, only : wkkct_glide
use climb_mod, only : wkkct_climb, wkkct_climb_irradiation
use isotropic_diffusion_mod, only : nabarro_herring_plus_coble_diffusion
use phase_material_mod, only : phase_material_type
use phase_material_array_mod, only : phase_material_array_type
use all_mighty_grid_mod, only : all_mighty_grid_type
use inelastic_strain_mod, only : inelastic_strain_base
use stiffness_base_mod, only : stiffness_base
use porosity_base_mod, only : porosity_base
USE, intrinsic :: iso_c_binding
use string_module, only : string_type
use bc_objects_mod, only : boundary_condition_array_type
use microstructure_info_mod, only : microstructure_info_type
implicit none


type(sim_all_macro_data_obj), pointer :: sim_all_macro_data => null()
class(inelastic_strain_base), pointer :: inelastic_strain_base_ptr => null()
class(stiffness_base), pointer :: elasticity_base_ptr => null()
class(porosity_base), pointer :: porosity_base_ptr => null()
class(phase_material_type),pointer :: phase_material_ptr
type(phase_material_array_type) :: phase_material_array
type(microstructure_info_type) :: microstrucutre_info

! specific offsets for the the third dimension
integer(C_INTPTR_T) :: npts3_rank, npts3_offset, npts3_start, npts3_end
integer :: n_gb_voxels

! the grid data container
class(all_grid_data), pointer :: common_grid_data_global => null()
class(fft_grid_data_type), pointer :: fft_grid_data_global => null()
type(all_mighty_grid_type), pointer :: all_mighty_grid_global => null()
type(string_type) :: microstructure_file_name, bc_file_name
class(boundary_condition_array_type), pointer :: the_bc_object => null()
! common griddata used by all modules
integer(k_int) , dimension(:,:,:), pointer :: jgrain => null()
integer(k_int) , dimension(:,:,:), pointer :: jphase => null()
integer(k_int) , dimension(:,:,:), pointer :: idgb => null()
real(k_real) , dimension(:,:,:), pointer :: tridistance => null()
real(k_real) , dimension(:,:,:,:), pointer :: gbnormals => null()
real(k_real) , dimension(:,:,:,:), pointer :: neibs => null()
real(k_real) :: B(3,3,6)
! tolerances
real(k_real) :: material_stress_abs_tol = 1e-7
real(k_real) :: imposed_stressBC_abs_tol = 1e-7
real(k_real) :: nl_strain_rate_rel_tol = 1e-5
real(k_real) :: nl_stress_rel_tol = 1e-5
real(k_real) :: tmarch_critical_strain_rate = 0.9
real(k_real) :: strain_rate_value_for_using_abs_tol = 1e-12
real(k_real) :: nl_strain_rate_abs_tol = 1e-15
real(k_real), parameter :: min_meaningful_strain_rate_to_check_convergence = 1e-16
integer :: max_nl_iterations = 0
integer :: max_material_iterations = 0

include "mechanics_variables.inc"
! include "vacancy_transport_variables.inc"

real(k_real) :: voxel_size(3), voxel_volume
real(k_Real), pointer :: tdot =>null(), &
                         tdot_min =>null(), &
                         tdot_max  =>null(), &
                         tdot_pre  =>null(), &
                         tdot_ref  =>null(), &
                         time => null()

real(k_real) :: vel_grad_macro(3,3),vel_grad_correction(3,3),rve_stress(3,3), rve_stress_t(3,3), IVAR
real(k_real) :: epav(3,3),edotpav(3,3),etotav(3,3), edottav(3,3)
real(k_real) :: eel(3,3),eeldot(3,3),eel_avg(3,3),eeldot_avg(3,3)
real(k_Real) :: avg_rve_climb_tensor(3,3),avg_rve_glide_tensor(3,3)
real(k_real) :: dis_grad_macro(3,3)

! ***** TO REMOVE BLOCK ******
INTEGER :: IINNERFAIL
INTEGER :: NPH
real(k_real) :: WGT



INTEGER ::  kperiodic_image, neighbors
! PPC
real(k_real) :: rve_stress_vm, &
                rve_stress_vm_old, &
                rve_equivalent_total_strain_rate, &
                rve_equivalent_total_strain_rate_old, &
                rve_equivalent_total_strain_rate_older, &
                ERRS_REL,ERRE_REL, ERRS_ABS, ERRE_ABS, ERRSBC, &
                ERRS_REL_MAX, ERRE_REL_MAX, ERRS_ABS_MAX, ERRE_ABS_MAX
logical :: update_texture=.true.,phase_hardening=.true.
! ***** TO REMOVE BLOCK ******
logical, pointer, dimension(:) :: IGAS => null()
! ***** end TO REMOVE BLOCK ******
!New output EL
logical :: output_err,output_fields,output_tex,output_ssc,printtimestamp,output_str
integer :: ssc_moviestep,err_moviestep,fld_moviestep,tex_moviestep,str_moviestep

! ELAS
real(k_real) :: c0(3,3,3,3),s0(3,3,3,3),c066(6,6), s066(6,6)
real(k_real) :: XLSEC_new(6,6), XLSEC(6,6), XLSEC_old(6,6),&
                XMSEC_new(6,6), XMSEC(6,6), XMSEC_old(6,6),&
                XLSEC_new_3333(3,3,3,3), XLSEC_3333(3,3,3,3), XLSEC_old_3333(3,3,3,3), &
                XMSEC_new_3333(3,3,3,3), XMSEC_3333(3,3,3,3), XMSEC_old_3333(3,3,3,3)

real(k_real) :: rve_equivalent_total_strain, evmp, dvmp,evmtot,evme,dvmtot,dvme
INTEGER :: IMICRO
integer :: ii,jj,kk,ll,outer_loop_iteration
real(k_real):: avg_rve_glide,avg_rve_climb

!!!! RESTART VARIABLES:
logical :: restart_from_dump=.false.
type(string_type) :: restart_file_name, dump_file_base_name, field_file_base_name
integer :: num_restart_file_to_keep, write_restart_file_every_n_steps, &
           write_field_files_every_n_steps

!!!!! debug variables
logical :: write_non_converged_mp_info_to_file=.false.

contains
  subroutine initGridData()
    use grid_data_types
    use mechanics_module_init
    use all_grid_data_mod
    implicit none


    integer :: nss
    ! this variables are here for now
    logical :: use_mechanics, use_vecancy_transport
    nss = 1
    use_mechanics = .TRUE.
    use_vecancy_transport = .TRUE.


    call common_grid_data_global%addVar("phase_fraction", generic_vector, additional_var_dimensions=(/microstrucutre_info%getNumPhases()/))
    call common_grid_data_global%addVar("grain_id", scalar_integer)
    call common_grid_data_global%addVar("gb_id", scalar_integer)
    call common_grid_data_global%addVar("triple_line_distance", scalar)
    call common_grid_data_global%addVar("gb_normals", real_vector)
    call common_grid_data_global%addVar("grain_neighbors", generic_vector , additional_var_dimensions=(/7/))

    if (use_mechanics) call getMechanicsGridVariables(common_grid_data_global)

    call all_mighty_grid_global%AMGAllocateGridVariables()

    call assingGridDataVariables()
  end subroutine

  subroutine assingGridDataVariables()
    implicit none
    logical :: use_mechanics, use_vecancy_transport
    use_mechanics = .TRUE.
    use_vecancy_transport = .TRUE.

    call common_grid_data_global%getScalarIntegerDataPointerByName("grain_id", jgrain)
    call common_grid_data_global%getScalarIntegerDataPointerByName("gb_id", idgb)
    call common_grid_data_global%getScalarDataPointerByName("triple_line_distance", tridistance)
    call common_grid_data_global%getRealVectorDataPointerByName("gb_normals", gbnormals)
    call common_grid_data_global%getGenericVectorDataPointerByName("grain_neighbors", neibs)

    if (use_mechanics) call assignMechanicsPointers()
    ! if (use_vecancy_transport) call assignVacancyTransportPointers()

  end subroutine

  subroutine assignMechanicsPointers()
    implicit none
    call common_grid_data_global%getTensor2DataPointerByName("velocity_gradient", vel_grad)
    call common_grid_data_global%getTensor2DataPointerByNameOld("velocity_gradient", vel_grad_old)
    call common_grid_data_global%getTensor2DataPointerByName("total_displacement_gradient", disgradtot)
    call common_grid_data_global%getTensor2DataPointerByNameOld("total_displacement_gradient", disgradtot_old)
    call common_grid_data_global%getTensor2DataPointerByName("cauchy_stress", grid_stress)
    call common_grid_data_global%getTensor2DataPointerByNameOld("cauchy_stress", grid_stress_old)
    call common_grid_data_global%getTensor2DataPointerByName("cauchy_stress_rate", grid_stress_rate)
    call common_grid_data_global%getTensor2DataPointerByNameOld("cauchy_stress_rate", grid_stress_rate_old)

    call common_grid_data_global%getTensor2DataPointerByName("elastic_strain", eelfield)
    call common_grid_data_global%getTensor2DataPointerByNameOld("elastic_strain", eelfield_old)
    call common_grid_data_global%getTensor2DataPointerByName("elastic_strain_rate", eelfield_rate)

    call common_grid_data_global%getTensor2DataPointerByName("plastic_strain", ept)
    call common_grid_data_global%getTensor2DataPointerByNameOld("plastic_strain", ept_old)
    call common_grid_data_global%getTensor2DataPointerByName("plastic_strain_rate", edotp)
    call common_grid_data_global%getTensor2DataPointerByNameOld("plastic_strain_rate", edotp_old)

    call common_grid_data_global%getTensor2DataPointerByName("total_strain", total_strain) !->WE NEED TO CHECK THIS
    call common_grid_data_global%getTensor2DataPointerByNameOld("total_strain", total_strain_old)
    call common_grid_data_global%getTensor2DataPointerByName("total_strain_rate", total_strain_rate) !->WE CAN just store the angles instead of the matrices if its more convenient
    call common_grid_data_global%getTensor2DataPointerByNameOld("total_strain_rate", total_strain_rate_old) !->WE CAN just store the angles instead of the matrices if its more convenient


    call common_grid_data_global%getMatrix66DataPointerByName("stiffness", stiffness66)
    call common_grid_data_global%getMatrix66DataPointerByNameOld("stiffness", stiffness66_old)

  end subroutine

END MODULE GLOBAL
