module inner_loop
use, intrinsic :: IEEE_ARITHMETIC
use kinds
use tensor_math_mod
use test_utils_mod
use print_utils_mod
use mpi_variables_mod, only : i_am_mpi_master
IMPLICIT NONE

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine evpal
  USE global, only: all_mighty_grid_global, IINNERFAIL, tdot, &
                    erre_rel, errs_rel, erre_abs, errs_abs,&
                    erre_rel_max, errs_rel_max, erre_abs_max, errs_abs_max,&
                    material_stress_abs_tol, max_material_iterations,&
                    grid_stress, grid_stress_rate, grid_stress_old, wgt, edotp, stiffness66, vel_grad, &
                    XLSEC_new, XLSEC, XLSEC_old, &
                    XLSEC_new_3333,  XLSEC_3333,  XLSEC_old_3333, &
                    XMSEC_new, XMSEC, XMSEC_old, &
                    XMSEC_new_3333, XMSEC_3333, XMSEC_old_3333, &
                    phase_material_array
  use change_tensor_basis, only : chg_basis_tensor2_to_vector6, chg_basis_vector6_to_tensor2, chg_basis_matrix66_to_tensor4
  use matrix_inversion, only : lapackInverseSymmetric, lapackInverseNonSymmetric
  use mpi_useful_routines_mod, only : MPISumMatrix, MPISumScalar, MPISumScalarInteger, MPIMAXScalar
  use tensor_math_mod, only : vectorNorm
  use log_file_mod, only : write_detailed_log_to_screen
  Implicit none
  real(k_real) :: R_new, R_old, alpha
  real(k_real), dimension(6) :: stress6, equilibrated_stress_old, &
                                equilibrated_stress6, equilibrated_strain_rate6, &
                                constituive_strain_rate6, &
                                edot_el, edot_pl, &
                                R_stress, R_strain, &
                                initial_residual, &
                                delta_stress_guess, &
                                stress_begin_NR, &
                                stress_end_NR


  real(k_real), dimension(6,6) :: xjacobinv, dconstituive_strain_rate6_dsigma, &
                                  dRstrain_dstress, dRstress_dstress, &
                                  S, & !-> local complaince matrix, &
                                  L_local, M_local, &
                                  XLSEC_proc
  real(k_real) :: erre_proc_rel, errs_proc_rel, erre_proc_abs, errs_proc_abs, &
                  erre_proc_rel_max, errs_proc_rel_max, erre_proc_abs_max, errs_proc_abs_max
  integer :: iter_stress, ix, iy, iz, N_VOXEL_FAIL_PROC, inv_info, fix_point_inner_loop_iter
  logical :: al_converged, run_diagnostic, fix_point_inner_loop_converged
  integer, parameter :: fix_point_inner_loop_max_iter=100
  integer ::  x_start_rank, x_end_rank, y_start_rank, y_end_rank, z_start_rank, z_end_rank
  
  
  call all_mighty_grid_global%AMGGetLoopLimitsRank(x_start_rank, x_end_rank, y_start_rank, y_end_rank, z_start_rank, z_end_rank)

  XLSEC_proc = 0._k_real
  N_VOXEL_FAIL_PROC = 0
  IINNERFAIL = 0
  erre_proc_rel = 0._k_real
  errs_proc_rel = 0._k_real
  erre_proc_abs = 0._k_real
  errs_proc_abs = 0._k_real
  erre_proc_rel_max = 0._k_real
  errs_proc_rel_max = 0._k_real
  erre_proc_abs_max = 0._k_real
  errs_proc_abs_max = 0._k_real

  DO iz=z_start_rank, z_end_rank; DO iy=y_start_rank, y_end_rank; DO  ix=x_start_rank, x_end_rank

      run_diagnostic =.FALSE.
      ! new stress update
666   call chg_basis_tensor2_to_vector6(grid_stress(:,:,ix,iy,iz), equilibrated_stress6)
      call chg_basis_tensor2_to_vector6(getSymmetricPart(vel_grad(:,:,ix,iy,iz)), equilibrated_strain_rate6)
      call chg_basis_tensor2_to_vector6(grid_stress_old(:,:,ix,iy,iz), equilibrated_stress_old)

      stress6 = equilibrated_stress6 ! -> init first stress guess
      constituive_strain_rate6 = 0._k_real
      call lapackInverseSymmetric(stiffness66(:,:,ix,iy,iz) ,S)
      fix_point_inner_loop_iter = 0
      fix_point_inner_loop_converged = .false.
      do while (fix_point_inner_loop_iter.lt.fix_point_inner_loop_max_iter.and.(.not.(fix_point_inner_loop_converged)))
      ! compute stress and strain norm to compute global errors
      iter_stress = 0
      al_converged= .false.
      alpha = 1._k_real
      fix_point_inner_loop_iter = fix_point_inner_loop_iter + 1
      stress_begin_NR = stress6
      call phase_material_array%updateStateVariablesInnerLoop(ix,iy,iz)
      do while(iter_stress.lt.max_material_iterations.and.(.not.(al_converged)))
        iter_stress = iter_stress+1
        
        ! if (ix*iy*iz.eq.1.and.i_am_mpi_master) write(*,*) "iter_stress", iter_stress
        !get the plastic strain rate and its derivatives
        call phase_material_array%getConstitutiveStrainRate(stress6, &
                                                            constituive_strain_rate6, &
                                                            dconstituive_strain_rate6_dsigma, &
                                                            edot_el, edot_pl, &
                                                            ix,iy,iz)

        ! if (ix*iy*iz.eq.1.and.i_am_mpi_master)then
        !   write(*,*) "bla"
        ! endif
        R_stress = ( stress6-equilibrated_stress6 )
        R_strain =  matmul(XLSEC, constituive_strain_rate6 - equilibrated_strain_rate6)

        R_New = vectorNorm(R_stress+R_strain)
        if ((iter_stress>1).and.(R_New.ge.R_old)) then
          alpha = alpha*0.5_k_real
          stress6 = stress6 - delta_stress_guess*0.5_k_real
          call phase_material_array%getConstitutiveStrainRate(stress6, &
                                                            constituive_strain_rate6, &
                                                            dconstituive_strain_rate6_dsigma, &
                                                            edot_el, edot_pl, &
                                                            ix,iy,iz)

          R_stress = ( stress6-equilibrated_stress6 )
          R_strain =  matmul(XLSEC, constituive_strain_rate6 - equilibrated_strain_rate6)
          R_New = vectorNorm(R_stress+R_strain)

        endif
          R_old = R_New

        if (run_diagnostic.and.write_detailed_log_to_screen) then
          write(*,*) "Innerloop iter_stress", iter_stress, " new stress ", stress6
          write(*,*) "Innerloop iter_stress", iter_stress, " Residual ", R_stress+R_strain, " Residual norm ", R_New
        endif
        if (iter_stress==1) initial_residual = abs(R_stress+R_strain)

        if ((R_New.le.material_stress_abs_tol) &
          .and.(iter_stress>1)) then ! we converged
          al_converged=.true.
          stress_end_NR = stress6
        else ! we didn't converged yet
          ! compute jacobian and next stress guess
          dRstress_dstress = I66
          dRstrain_dstress = matmul(XLSEC, dconstituive_strain_rate6_dsigma)
          call lapackInverseNonSymmetric(dRstress_dstress+dRstrain_dstress, xjacobinv, inv_info)
          if (inv_info.ne.0) exit
          ! compute new stress guess
          delta_stress_guess = -matmul(xjacobinv,R_stress+R_strain)*alpha
          if (run_diagnostic.and.write_detailed_log_to_screen) then
            write(*,*) "Innerloop iter_stress", iter_stress, " alpha ", alpha
            write(*,*) "Innerloop iter_stress", iter_stress, "dRstrain_dstress", dRstrain_dstress
            write(*,*) "Innerloop iter_stress", iter_stress, "XLSEC", XLSEC
            write(*,*) "Innerloop iter_stress", iter_stress, "S/tdot", S/tdot
            write(*,*) "Innerloop iter_stress", iter_stress, "dconstituive_strain_rate6_dsigma", dconstituive_strain_rate6_dsigma
            write(*,*) "Innerloop iter_stress", iter_stress, "jacob", dRstress_dstress+dRstrain_dstress
            write(*,*) "Innerloop iter_stress", iter_stress, "xjacobinv", xjacobinv
            write(*,*) "Innerloop iter_stress", iter_stress, " delta_stress_guess ", delta_stress_guess, " delta stress norm ", vectorNorm(delta_stress_guess)
          endif
          if ((vectorNorm(delta_stress_guess)<material_stress_abs_tol).and.(iter_stress>1).and.( R_New.lt.(material_stress_abs_tol*100._k_real))) then
            al_converged=.true.
          else
            stress6 = stress6 + delta_stress_guess
          endif
        endif ! END Newton raphson convergence checks
      enddo ! END Newton raphson
        if (al_converged) then
          ! R_stress = ( stress_begin_NR-stress_end_NR)
          ! R_New = vectorNorm(R_stress)
          ! if (run_diagnostic.and.write_detailed_log_to_screen) then
          !   write(*,*) "Innerloop fix point iteration", fix_point_inner_loop_max_iter, " Residual ", R_stress, " Residual norm ", R_New
          !   write(*,*) "Innerloop fix point iteration", fix_point_inner_loop_max_iter, " stress_begin_NR ", stress_begin_NR
          !   write(*,*) "Innerloop fix point iteration", fix_point_inner_loop_max_iter, " stress_end_NR ", stress_end_NR
          ! endif
          fix_point_inner_loop_converged = .TRUE.
          !  ((R_New.lt.material_stress_abs_tol).and.fix_point_inner_loop_iter.ge.2)
        else 
          exit
        endif
      enddo ! fix point NR

      if(fix_point_inner_loop_converged) then
        ! if (i_am_mpi_master.and.(ix*iy*iz==1)) write(*,*) "fix point iteration converged in ", fix_point_inner_loop_iter,  "iterations"
        call updateALErrors(stress6, constituive_strain_rate6, equilibrated_stress6, equilibrated_strain_rate6, &
            errs_proc_abs, erre_proc_abs, errs_proc_rel, erre_proc_rel, &
            errs_proc_abs_max, erre_proc_abs_max, errs_proc_rel_max, erre_proc_rel_max)
        call chg_basis_vector6_to_tensor2(stress6, grid_stress(:,:,ix,iy,iz))



        call chg_basis_vector6_to_tensor2(edot_pl, edotp(:,:,ix,iy,iz))
        grid_stress_rate(:,:,ix,iy,iz) = (grid_stress(:,:,ix,iy,iz) - grid_stress_old(:,:,ix,iy,iz))/tdot
        M_local = dconstituive_strain_rate6_dsigma
        call lapackInverseSymmetric(M_local, L_local, inv_info )
        if (inv_info.ne.0) al_converged = .false.
        XLSEC_proc = XLSEC_proc + L_local*wgt
      endif

      if (.not.fix_point_inner_loop_converged)then
        if (.not.run_diagnostic) then
          run_diagnostic = .TRUE.
          GOTO 666
        else
          run_diagnostic=.FALSE.
        endif
        N_VOXEL_FAIL_PROC = N_VOXEL_FAIL_PROC+1
        if (run_diagnostic.and.write_detailed_log_to_screen) then
          WRITE(*,*) "Warning: innerloop failed for voxel ", ix, iy, iz
          WRITE(*,*) "Initial Residual ", initial_residual
          WRITE(*,*) "Final Residual ", abs(R_stress+R_strain)
        endif
      endif

  END DO; END DO; END DO

  call MPISumScalarInteger(N_VOXEL_FAIL_PROC, IINNERFAIL)
  if (IINNERFAIL.eq.0) then
    call MPISumMatrix(XLSEC_proc, XLSEC_new)
    call MPISumScalar(erre_proc_rel, erre_rel)
    call MPISumScalar(errs_proc_rel, errs_rel)
    call MPISumScalar(erre_proc_abs, erre_abs)
    call MPISumScalar(errs_proc_abs, errs_abs)
    call MPIMAXScalar(erre_proc_rel_max, erre_rel_max)
    call MPIMAXScalar(errs_proc_rel_max, errs_rel_max)
    call MPIMAXScalar(erre_proc_abs_max, erre_abs_max)
    call MPIMAXScalar(errs_proc_abs_max, errs_abs_max)

    ! compute XLSEC_new
    call lapackInverseSymmetric(XLSEC_new, XMSEC_new)
    call chg_basis_matrix66_to_tensor4(XLSEC_new, XLSEC_new_3333)
    call chg_basis_matrix66_to_tensor4(XMSEC_new, XMSEC_new_3333)

    ! XLSEC = XLSEC_new
    ! XMSEC = XMSEC_new
    ! XMSEC_3333 = XMSEC_new_3333
    ! XLSEC_3333 = XLSEC_new_3333

    call phase_material_array%updateStateVariablesOuterLoop()
  else
    XLSEC = XLSEC_old
    XMSEC = XMSEC_old
    XMSEC_3333 = XMSEC_old_3333
    XLSEC_3333 = XLSEC_old_3333
  endif

end subroutine evpal

subroutine updateALErrors(stress_sol, strain_rate_sol, stress_initial, strain_rate_inital, &
                           errs_proc_abs, erre_proc_abs, errs_proc_rel, erre_proc_rel, &
                           errs_proc_abs_max, erre_proc_abs_max, errs_proc_rel_max, erre_proc_rel_max)
  use global, only : WGT
  implicit none
  real(k_real), intent(in), dimension(6) :: stress_sol, strain_rate_sol, stress_initial, strain_rate_inital
  real(k_real), intent(inout) :: errs_proc_abs, erre_proc_abs, errs_proc_rel, erre_proc_rel, &
                                 errs_proc_abs_max, erre_proc_abs_max, errs_proc_rel_max, erre_proc_rel_max
  real(k_real) :: norm_stress, norm_strain_rate
  norm_stress = vectorNorm(stress_sol-stress_initial)
  norm_strain_rate = vectorNorm(strain_rate_sol - strain_rate_inital)
  !and update the global error
  errs_proc_abs = errs_proc_abs + norm_stress*WGT
  erre_proc_abs = erre_proc_abs + norm_strain_rate*WGT
  errs_proc_rel = errs_proc_rel + norm_stress/vectorNorm(0.5_k_real*(stress_sol+stress_initial))*WGT
  erre_proc_rel = erre_proc_rel + norm_strain_rate/vectorNorm(0.5_k_real*(strain_rate_sol+strain_rate_inital))*WGT
  errs_proc_abs_max = max(errs_proc_abs_max, norm_stress )
  erre_proc_abs_max = max(erre_proc_abs_max, norm_strain_rate )
  errs_proc_rel_max = max(errs_proc_rel_max, norm_stress/vectorNorm(0.5_k_real*(stress_sol+stress_initial)) )
  erre_proc_rel_max = max(erre_proc_rel_max, norm_strain_rate/vectorNorm(0.5_k_real*(strain_rate_sol+strain_rate_inital)) )
end subroutine

subroutine get_smacro
  USE global, only: grid_stress, vel_grad_correction, XMSEC_new_3333,&
                    rve_stress, rve_stress_vm, &
                    ERRSBC, the_bc_object, vel_grad_macro
  use mpi_useful_routines_mod, only : MPIAverageGridMatrix
  use solve_linear_system, only : solve_mixed_linear_system_evpbc
  use tensor_math_mod, only : getSymmetricPart
  use voigt_indicial_conversion_mod, only : Tensor2ToVectorVoigt
  implicit none
    real(k_real), dimension(3,3) :: stress_bc, vel_grad_bc, strain_rate_solution, stress_solution
    integer, dimension(3,3) ::  istress_rate_components, i_vel_grad_comp
    integer, dimension(6) :: i_vel_grad_comp_6, istress_rate_components_6
    call MPIAverageGridMatrix(grid_stress, rve_stress)


    call the_bc_object%computeCurrentBC(stress=stress_bc, istress_rate_components=istress_rate_components, &
                                        vel_grad=vel_grad_bc, ivel_grad_components=i_vel_grad_comp)
    ERRSBC=tensor2Norm(istress_rate_components*(stress_bc-rve_stress))

    i_vel_grad_comp_6 =Tensor2ToVectorVoigt(i_vel_grad_comp)
    istress_rate_components_6 =Tensor2ToVectorVoigt(istress_rate_components)

    call solve_mixed_linear_system_EVPBC(getSymmetricPart(vel_grad_macro), rve_stress, &
                                i_vel_grad_comp_6, istress_rate_components_6, &
                                getSymmetricPart(vel_grad_bc), stress_bc, &
                                XMSEC_new_3333, strain_rate_solution, stress_solution)

    vel_grad_correction = strain_rate_solution - getSymmetricPart(vel_grad_macro)

    rve_stress_vm = computeVMEquivalentStress(rve_stress)

end subroutine get_smacro


! subroutine writeFailedIntegrationPointTofile(ix,iy,iz, equilibrated_stress6, equilibrated_stress_old, equilibrated_strain_rate6, S)
!   use mpi_variables_mod, only : mpi_rank
!   use global, only : grid_data, time, tdot, XLSEC
!   implicit none
!   integer, intent(in) :: ix,iy,iz
!   integer :: i
!   real(k_real), dimension(6), intent(in) :: equilibrated_stress6, equilibrated_stress_old, equilibrated_strain_rate6
!   real(k_real), dimension(6,6), intent(in) :: S
!   integer :: file_id
!   logical:: ex
!   character(len=100) :: f_string, fname

!   write(f_string,*) mpi_rank
!   file_id = 167 + mpi_rank
!   fname = "dump_points_rank_"//trim(adjustl(f_string))//".txt"
!   INQUIRE (FILE=trim(adjustl(fname)), EXIST=ex)
!   if (ex) then
!     open(unit = file_id, file = trim(adjustl(fname)), status="old", position="append", action="write")
!   else
!      open(unit = file_id, file = trim(adjustl(fname)), status="new", action="write")
!   endif
!   write(file_id, *) "mpi_rank: ", mpi_rank
!   write(file_id, *) "location: ", ix,iy,iz
!   write(file_id, *) "Inital Values"
!   write(file_id, *) "time"
!   write(file_id, *) time
!   write(file_id, *) "tdot"
!   write(file_id, *) tdot
!   write(file_id, *) "XLSEC"
!   do i=1,6
!   write(file_id, *) XLSEC(:,i)
!   enddo
!   write(file_id, *) "initial_stress_guess"
!   write(file_id, *) equilibrated_stress6
!   write(file_id, *) "converged_stress_previous_time_step"
!   write(file_id, *) equilibrated_stress_old
!   write(file_id, *) "equilibrated_strain_rate6"
!   write(file_id, *) equilibrated_strain_rate6
!   write(file_id, *) "compliance"
!   do i=1,6
!   write(file_id, *) S(:,i)
!   enddo
!   write(file_id, *) ""
!   write(file_id, *) "GridData"
!   call grid_data%DumpMaterialPointValuesToTextFile(ix, iy, iz, file_id)
!   write(file_id, *) ""
!   close(unit=file_id)
! end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module inner_loop
