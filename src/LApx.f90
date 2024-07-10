PROGRAM EVPFFT
  use mpi, only : MPI_Barrier, MPI_COMM_WORLD
  use mpi_variables_mod, only : init_mpi, finalize_mpi, mpi_rank, my_mpi_err
  use mpi_useful_routines_mod, only : MPILogicalOR, MPIAverageGridMatrix, MPISumMatrix, AmIMPIMaster
  USE global, only : all_mighty_grid_global, sim_all_macro_data, the_bc_object, phase_material_array, restart_file_name, field_file_base_name, &
  fft_grid_data_global, restart_from_dump, rve_equivalent_total_strain_rate, rve_equivalent_total_strain_rate_old, rve_equivalent_total_strain_rate_older, &
  rve_stress, c066, s066, c0, s0, rve_stress_t, grid_stress, erre_abs, erre_abs_max, erre_rel, erre_rel_max, errsbc, &
  XLSEC_new, XLSEC, XLSEC_old,&
  XMSEC_new, XMSEC, XMSEC_old,&
  XLSEC_new_3333, XLSEC_3333, XLSEC_old_3333, &
  XMSEC_new_3333, XMSEC_3333, XMSEC_old_3333, &
  vel_grad, errs_abs, errs_abs_max, errs_rel_max, IINNERFAIL, imicro, &
  imposed_stressBC_abs_tol, material_stress_abs_tol, max_material_iterations, max_nl_iterations, min_meaningful_strain_rate_to_check_convergence, &
  nl_strain_rate_abs_tol, nl_strain_rate_rel_tol, nl_stress_rel_tol, outer_loop_iteration, output_err,&
  phase_hardening, rve_stress_vm, rve_stress_vm_old, strain_rate_value_for_using_abs_tol, tdot, time, tmarch_critical_strain_rate, total_strain_rate, update_texture, &
  errs_rel, vel_grad_correction, vel_grad_macro, voxel_size, voxel_volume,wgt, write_field_files_every_n_steps, dump_file_base_name, write_restart_file_every_n_steps, num_restart_file_to_keep


  USE fileIO
  USE inner_loop
  ! USE vacancy_transport
  USE, intrinsic :: iso_c_binding
  use matrix_inversion
  use change_tensor_basis
  use kinds
  use tensor_math_mod
  use elasticity_mod
  use porosity_base_mod
  ! use texture_mod
  use gb_normals_mod, only : initGBproperties
  use read_microstructure_file_mod, only : read_microstructure_from_file
  ! use texture_mod, only : update_rotation_matrix, update_stiffness
  use print_utils_mod
  use compute_averages_mod, only : UpdateRateFields, StressGuessNextOuterLoop
  use time_march_mod, only : acceptRejectSolution
  use bc_mod, only : computeInitialGuess
  use hdf5, only : h5open_f, h5close_f
  use log_file_mod
  use number_to_string_mod
  use termination_criterion_mod, only : checkSimulationTermination, terminate_simulation_global
  use bc_objects_mod, only : i_vel_grad
  use fft_grid_data_mod, only : fft_frequency_vector_norm, &
                                fft_frequency_tensor_K
  use csv_writer_mod, only : csv_writer

  implicit none

real(k_real) :: a(3,3),g1(3,3,3,3)
complex(k_real) :: complex_tensor2(3,3)

real(k_real) :: errsbc_abs,&  ! aboslute applied stress error
                errsbc_rel    ! relative applied stress error

integer, parameter :: dump_restart_nsteps = 2
integer, parameter :: n_dump_to_keep = 3


!! AC
integer :: last_converged_iterations
!!IAL2
integer :: IAL=1, ierr

integer :: ix, iy, iz, ix_global,iy_global,iz_global


logical :: converged = .false.

! time march variables
real(k_real) :: accept_reject_dt, time_read, dt_read
logical :: accept_solution = .true.
logical :: new_bc_block = .true.
logical :: bc_block_completed = .false.
! temp varaibles for checking tolerances
real(k_real) :: avg_tot_strain_rate(3,3), &
                avg_tot_strain_rate_norm

integer :: x_start_rank, x_end_rank, y_start_rank, y_end_rank, z_start_rank, z_end_rank
integer :: npts1, npts2, npts3


call init_mpi()
call h5open_f(ierr)
allocate(all_mighty_grid_global)
allocate(sim_all_macro_data)
allocate(the_bc_object)
call sim_all_macro_data%init()
! call the_bc_object%init()

call input ! while reading we initialize griddata and fftw arrays
call all_mighty_grid_global%AMGSetDumpForRestartOptions(num_restart_file_to_keep, write_restart_file_every_n_steps, dump_file_base_name%getString())
call all_mighty_grid_global%AMGSetWriteFieldOptions(write_field_files_every_n_steps, field_file_base_name%getString())
call the_bc_object%initBC()


call all_mighty_grid_global%AMGGetLoopLimitsRank(x_start_rank, x_end_rank, y_start_rank, y_end_rank, z_start_rank, z_end_rank)
call all_mighty_grid_global%AMGGetGlobalGridDimension(npts1, npts2, npts3)
call phase_material_array%updatePhaseParameters()
call phase_material_array%updateRotationMatrix()
call phase_material_array%updateStiffness()
call phase_material_array%updateAverageStiffness(c066, s066, c0, s0)
call phase_material_array%initTexture()
call phase_material_array%initStateVariables()

if(restart_from_dump) then
  call all_mighty_grid_global%AMGReloadFromDump(restart_file_name%getString(), time_read, dt_read, imicro, XLSEC, XLSEC_old)
  write(*,*) mpi_rank, "time_read, dt_read", time_read, dt_read
  write(*,*) mpi_rank, "imicro ", imicro
  call MPI_Barrier(MPI_COMM_WORLD, my_mpi_err)
  call sim_all_macro_data%sim_time%setTimeFromReload(time_read, dt_read)
  ! compute XLSEC

  call lapackInverseSymmetric(XLSEC, XMSEC)
  call lapackInverseSymmetric(XLSEC_OLD, XMSEC_OLD)
  call chg_basis_matrix66_to_tensor4(XLSEC, XLSEC_3333)
  call chg_basis_matrix66_to_tensor4(XLSEC_old, XLSEC_old_3333)
  call chg_basis_matrix66_to_tensor4(XMSEC, XMSEC_3333)
  call chg_basis_matrix66_to_tensor4(XMSEC_OLD, XMSEC_old_3333)
endif

call all_mighty_grid_global%AMGupdateAllStatefulVariables() ! saving state variables old values so we have teh references we need
if (.not.(sim_all_macro_data%sim_time%initialized)) error stop "sim_time not initialized, Abort!"
call sim_all_macro_data%sim_macro_field_averages%updateAverageQuantities()

! call writeDDaverages(dd_csv_writer, rho_m, rho_w, time, .TRUE.)
call writeToLogFile(" ", .TRUE.)
call phase_material_array%writeInelastiStrainRateToFile()
call phase_material_array%writeAverageQuantitiesCSV(write_headers=.TRUE.)
call phase_material_array%writeAverageQuantitiesCSV(write_headers=.FALSE.)
call sim_all_macro_data%sim_time%resetDtToDt0()
call sim_all_macro_data%sim_time%advanceTime()



!check tolerances
if(material_stress_abs_tol.le.0._k_real) error stop "the material_stress_abs_tol must be a small positive number. Check input file"
if(nl_strain_rate_rel_tol.le.0._k_real) error stop "the non-linear strain relative tolerance must be a small positive number. Check input file"
if(nl_stress_rel_tol.le.0._k_real) error stop "the non-linear stress relative tolerance must be a small positive number. Check input file"
if(imposed_stressBC_abs_tol.le.0._k_real) error stop "the imposed stress absolute tolerance must be a small positive number. Check input file"

! check iterations
if(max_material_iterations.le.0) error stop "the max number of material iterations should be a positive value. Check input file"
if(max_nl_iterations.le.0._k_real) error stop "the max number of non linear iterations should be a positive value. Check input file"
if (tmarch_critical_strain_rate.lt.1e-6_k_real) error stop "the critical strian rate for TIME_MARCH must be gt 1e-6. Check input file"



!here e start initializing variables TODO ! this should be between the reloading and first step
vel_grad_macro = i_vel_grad! in rate form - MAK,9/15/2020

if(restart_from_dump) call StressGuessNextOuterLoop
! HERE BEGINS THE START OVER THE LOOP
!stop
!stop
!!! HERE IS THE LOOP OVER THE BOUNDARY CONDITIONS (IN THIS CASE STEPS)
! call grid_data%AGDDumpForRestart(time, accept_reject_dt, imicro, XLSEC, XLSEC_old, force_write=.TRUE., write_only_current_value=.TRUE.)
accept_reject_dt = tdot
last_converged_iterations = max_nl_iterations

DO while(.not.(terminate_simulation_global))
  imicro = imicro +1
  if ( AmIMPIMaster()) then
     write(*,*) '***************************************************'
     write(*,*) 'STEP = ',imicro
   endif

   if (tdot.le.0._k_real) error stop "tdot is zero before  ITIMEMARCH"



if(new_bc_block) then
  call sim_all_macro_data%sim_time%resetDtToDt0()
else
  call sim_all_macro_data%sim_time%setDeltaTime(accept_reject_dt)
endif
   call sim_all_macro_data%sim_time%updateCurrentTime()
   call the_bc_object%ForceTime()

   call the_bc_object%computeCurrentBC( )

   if ( AmIMPIMaster()) write(*,*) "-------------------------------------------------------------"
   if ( AmIMPIMaster()) write(*,*) " TIME_MARCH DT", tdot, " time", time
   if ( AmIMPIMaster()) write(*,*) "-------------------------------------------------------------"



    call the_bc_object%printCurrentImposedBC()

   log_message = "STEP = " // trim(adjustl(int2string(imicro))) // " NEW DT FROM TIME MARCH is " // trim(adjustl(real2string(tdot)))
   call writeToLogFile(log_message)
   vel_grad_macro = i_vel_grad!

if(new_bc_block) then
  call computeInitialGuess()
endif
call phase_material_array%updateStiffness()
call phase_material_array%updateAverageStiffness(c066, s066, c0, s0)

if(imicro.eq.1) then

   ! initialize compute XLSEC_new and XMSEC_new
   XLSEC = c066*tdot
   XLSEC_3333 = c0*tdot
   XMSEC_3333 = s0/tdot
   ! initialize old x*sec
   XLSEC_old = XLSEC
   XMSEC_old_3333 = XMSEC_3333
   XLSEC_old_3333= XLSEC_3333
endif

   ! at the beginning of each step we reset delta values
   vel_grad_correction=0._k_real

   outer_loop_iteration=0
   erre_abs = 2*nl_strain_rate_rel_tol
   erre_rel = 2*nl_strain_rate_rel_tol
   errs_abs = 2*nl_stress_rel_tol
   errs_rel = 2*nl_stress_rel_tol
   ERRSBC=2*imposed_stressBC_abs_tol

   IINNERFAIL = 0
   converged = .false.
   !write(*,*) 'parallel test_3'
   DO WHILE (.not.(converged))
      ! LOOP OVER ITERATIONS IN A GIVEN STEP
      outer_loop_iteration=outer_loop_iteration+1

      !! AC
101   IF(IINNERFAIL.NE.0.OR.(.not.(ACCEPT_SOLUTION)).or.(outer_loop_iteration.EQ.max_nl_iterations+1)) then
          log_message = "STEP = " // trim(adjustl(int2string(imicro))) // " --LOOPING BACK: "
          if (IINNERFAIL.gt.0.) then
            log_message = TRIM(ADJUSTL(log_message)) // trim(adjustl(int2string(IINNERFAIL))) // " VOXELS DID NOT CONVERGED IN THE INNER LOOP. REDUCING TIME STEP FROM tdot = " // trim(adjustl(real2string(tdot)))
            if( AmIMPIMaster()) WRITE(*,*) trim(adjustl(log_message))
            call sim_all_macro_data%sim_time%reduceDeltaTime(2._k_real)
          endif
          if (.not.(ACCEPT_SOLUTION)) then
            log_message = TRIM(ADJUSTL(log_message)) // "AFTER STATE VARIABLES UPDATE. REDUCING TIME STEP FROM tdot = " // trim(adjustl(real2string(tdot)))
            if( AmIMPIMaster()) WRITE(*,*) trim(adjustl(log_message))
            call sim_all_macro_data%sim_time%setDeltaTime(accept_reject_dt)
            accept_solution = .true.
          endif
          if (outer_loop_iteration.EQ.max_nl_iterations+1) then
            log_message = TRIM(ADJUSTL(log_message)) // "MAXIMUM NUMBER OF ALLOWED NON LINEAR ITERATION EXCEEDED. REDUCING TIME STEP FROM tdot = " // trim(adjustl(real2string(tdot)))
            if( AmIMPIMaster()) WRITE(*,*) trim(adjustl(log_message))
            call sim_all_macro_data%sim_time%reduceDeltaTime(2._k_real)
          endif

          if( AmIMPIMaster()) WRITE(*,*) 'NEW TIME STEP IS: ', tdot
          call sim_all_macro_data%sim_time%updateCurrentTime()
          log_message = trim(adjustl(log_message)) // " TO tdot = " // trim(adjustl(real2string(tdot))) // " NEW time = " // trim(adjustl(real2string(time)))
          call writeToLogFile(log_message)
          ! reset everything stress to previous solution
          call all_mighty_grid_global%AMGResetAllStatefulVariables()

          ! reset X*sec
          XLSEC = XLSEC_old
          XMSEC = XMSEC_old
          XMSEC_3333 = XMSEC_old_3333
          XLSEC_3333 = XLSEC_old_3333
          rve_equivalent_total_strain_rate= rve_equivalent_total_strain_rate_old
          rve_stress_vm = rve_stress_vm_old
          rve_stress(:,:) = rve_stress_t(:,:)

           call the_bc_object%ForceTime()
           call the_bc_object%computeCurrentBC( )
           call the_bc_object%printCurrentImposedBC()

           IF (iMicro.eq.1.or.new_bc_block)  then
             call computeInitialGuess()
           else
             rve_equivalent_total_strain_rate= rve_equivalent_total_strain_rate_old
             rve_stress_vm = rve_stress_vm_old
             rve_stress = rve_stress_t
             call StressGuessNextOuterLoop()
           endif
          outer_loop_iteration = 1
          vel_grad_correction=0._k_real
    endif



    call MPIAverageGridMatrix(vel_grad, vel_grad_macro)
      do  iz = z_start_rank, z_end_rank
         do  iy = y_start_rank, y_end_rank 
            do  ix = x_start_rank, x_end_rank
              select case (IAL)
              case (1)
                call fft_grid_data_global%setFFTDataAtPointTensor2(ix,iy,iz, grid_stress(:,:,ix,iy,iz ))
              case default
                error stop "acceptable values for IAL is 1"
              end select
            end do
         end do
      end do

      ! fourier transform of stress field
      call fft_grid_data_global%executeDFTTensor2()

      do  iz = z_start_rank, z_end_rank
        do  iy = y_start_rank, y_end_rank 
           do  ix = x_start_rank, x_end_rank

              ix_global = all_mighty_grid_global%AMGXIndexdxLocal2Global(ix)
              iy_global = all_mighty_grid_global%AMGYIndexdxLocal2Global(iy)
              iz_global = all_mighty_grid_global%AMGZIndexdxLocal2Global(iz)

               if(ix_global.eq.(npts1/2+1).or.iy_global.eq.(npts2/2+1).or.iz_global.eq.(npts3/2+1)) then
                  g1=-XMSEC_3333
               else
                  if (fft_frequency_vector_norm(ix,iy,iz).ne.0._k_real) then
                    a = T4ijkl_T2jl(XLSEC_3333, fft_frequency_tensor_K(:,:,ix,iy,iz)) ! c0 is the viscous stiffness [MPa*s] so a is kind of a stress rate
                    call matrixInverseSymmetric(a)
                    g1 = -1._k_real*T2ik_T2jl(a, fft_frequency_tensor_K(:,:,ix,iy,iz))
                  else
                    g1 = 0._k_real
                  end if
               end if

               if (fft_frequency_vector_norm(ix,iy,iz).ne.0._k_real) then
                complex_tensor2 = fft_grid_data_global%getFFTDataAtPointComplexTensor2(ix,iy,iz)
               else
                 complex_tensor2 = cmplx(0.,0.,k_real)
               endif

               call fft_grid_data_global%setFFTDataAtPointTensor2(ix,iy,iz, T4ijkl_T2kl(g1, complex_tensor2))
            end do ! ix
         end do ! iy
      end do ! iz

      ! inverse fourier transform of complex strain field
      call fft_grid_data_global%executeInverseDFTTensor2()


      do  iz = z_start_rank, z_end_rank
        do  iy = y_start_rank, y_end_rank 
           do  ix = x_start_rank, x_end_rank
              select case (IAL)
              case (1)
                vel_grad(:,:,ix,iy,iz) = vel_grad(:,:,ix,iy,iz) + vel_grad_correction +&
                     fft_grid_data_global%getFFTDataAtPointRealTensor2(ix,iy,iz)*wgt
                total_strain_rate(:,:,ix,iy,iz) = getSymmetricPart(vel_grad(:,:,ix,iy,iz))
              case default
                error stop "acceptable values for IAL are 1 or 2"
              end select
            end do
         end do
      end do

      IINNERFAIL = 0
      call evpal ! here ew call the material
      call get_smacro

      errsbc_abs = ERRSBC
      errsbc_rel = errsbc_abs/rve_stress_vm
167 format (i5,10(E16.7))

      if ( AmIMPIMaster()) then
      if (outer_loop_iteration == 1 ) write(*,*) "Iter# "//&
                                 " ERRE_REL       "//&
                                 " ERRS_REL       "//&
                                 " ERRSBC         "//&
                                 " ERRE_ABS       "//&
                                 " ERRS_ABS       "//&
                                 " ERRSBC_ABS     "//&
                                 " ERRE_REL_MAX   "//&
                                 " ERRS_REL_MAX   "//&
                                 " ERRE_ABS_MAX   "//&
                                 " ERRS_ABS_MAX   "
      write(*,167) outer_loop_iteration,erre_rel,errs_rel,ERRSBC_rel, erre_abs, errs_abs, ERRSBC_abs, erre_rel_max, errs_rel_max, erre_abs_max, errs_abs_max
      endif
      ! CHECK OUTER LOOP CONVERGENCE
      converged = .true.

      !CHECK STRAIN RATE CONVERGENCE:

      ! compute average total strain rate
      call MPIAverageGridMatrix(total_strain_rate, avg_tot_strain_rate)


      avg_tot_strain_rate_norm = tensor2Norm(avg_tot_strain_rate)

      ! strain close or below machine tolerance
      if (avg_tot_strain_rate_norm.lt.min_meaningful_strain_rate_to_check_convergence) then
        ! strain rate below meaningful strain rate
        write(*,*) "using abs strain tolerance "
        write(*,*) "avg_tot_strain_rate_norm", avg_tot_strain_rate_norm
        converged = converged.and.converged
      else if (avg_tot_strain_rate_norm.lt.strain_rate_value_for_using_abs_tol) then
        ! using absolute tolerance
        write(*,*) "using abs strain tolerance "
        write(*,*) avg_tot_strain_rate_norm.lt.strain_rate_value_for_using_abs_tol
        write(*,*) "avg_tot_strain_rate_norm", avg_tot_strain_rate_norm
        converged = converged.and.( erre_abs.le.nl_strain_rate_abs_tol )
      else
        ! using relative tolerance
        converged = converged.and.( erre_rel.le.nl_strain_rate_rel_tol )
      endif


      ! check stress convergence
      converged = converged.and.( errs_rel.le.nl_stress_rel_tol )

      if ((.not.converged).and.(outer_loop_iteration.eq.max_nl_iterations).and.(errs_abs.lt.1e-6)) &
      converged = .TRUE.

      ! check BC convergence
      converged = converged.and.( errsbc_abs.le.imposed_stressBC_abs_tol )
      converged = converged.and.(outer_loop_iteration.ge.2)

      if (converged) then ! update state variables and write outputs


        ! what are we doing here?
        IF((update_texture).and.imicro.gt.1) THEN
          voxel_size(1) = voxel_size(1)*(1._k_real+vel_grad_macro(1,1)*tdot)
          voxel_size(2) = voxel_size(2)*(1._k_real+vel_grad_macro(2,2)*tdot) 
          voxel_size(3) = voxel_size(3)*(1._k_real+vel_grad_macro(3,3)*tdot)  
          call fft_grid_data_global%initFrequencyVectorAndTensor(voxel_size)
          voxel_volume = product(voxel_size)
        ENDIF

        ! update staggered state variables including porosity
        if (phase_hardening) then
         call phase_material_array%updateStateVariablesStaggered()
        endif

        ! update texture state variables
        if (update_texture) then

         call phase_material_array%updateRotationMatrix()
        !  call update_rotation_matrix !updating euler angles
        !  call update_stiffness(update_average_stiffness=.true.) ! rotate the stiffness
         call phase_material_array%updateStiffness()
         call phase_material_array%updateAverageStiffness(c066, s066, c0, s0)
         call phase_material_array%initTexture ! update shcmid and climb tensor
        else ! for all otehr cases we just update the stiffness
          call phase_material_array%updateStiffness()
          call phase_material_array%updateAverageStiffness(c066, s066, c0, s0)
        endif
        call UpdateRateFields()
        call acceptRejectSolution(accept_solution, accept_reject_dt)
        if(.not.accept_solution) goto 101

      endif
END DO !DO WHILE in OUTER LOOP

new_bc_block = .false.
bc_block_completed = .false.
! update xlsec
XLSEC = XLSEC_new
XMSEC = XMSEC_new
XMSEC_3333 = XMSEC_new_3333
XLSEC_3333 = XLSEC_new_3333
!store xlsec
XLSEC_old = XLSEC_new
XMSEC_old = XMSEC_new
XMSEC_old_3333 = XMSEC_new_3333
XLSEC_old_3333= XLSEC_new_3333

! check termination criteria
call checkSimulationTermination()

call the_bc_object%checkIfBlockCompleted(bc_block_completed)
call all_mighty_grid_global%AMGDumpForRestart(time, accept_reject_dt, imicro, XLSEC, XLSEC_old, &
      force_write=bc_block_completed.or.terminate_simulation_global, write_only_current_value=.FALSE.)
call all_mighty_grid_global%AMGDumpForRestart(time, accept_reject_dt, imicro, XLSEC, XLSEC_old, &
      force_write=bc_block_completed.or.terminate_simulation_global, write_only_current_value=.TRUE.)

! if we reacehd this point we are ready to go to the next step
! update state varaiables
call all_mighty_grid_global%AMGupdateAllStatefulVariables()
call sim_all_macro_data%sim_macro_field_averages%updateAverageQuantities()

! start output
call phase_material_array%writeInelastiStrainRateToFile()
! call writeDDaverages(dd_csv_writer, rho_m, rho_w, time, .FALSE.)
call phase_material_array%writeAverageQuantitiesCSV(write_headers=.FALSE.)

log_message = "STEP = " // trim(adjustl(int2string(imicro))) // &
              " dt = " // trim(adjustl(real2string(tdot))) //  &
              " time = " // trim(adjustl(real2string(time))) //  &
              " CONVERGED in " //  trim(adjustl(int2string(outer_loop_iteration))) // " iterations"
call writeToLogFile(log_message)

if (.not.terminate_simulation_global) &
  call the_bc_object%updateBC(terminate_simulation_global, log_message, new_bc_block)

if (terminate_simulation_global) then
  call writeToLogFile(log_message)
  call MPI_Barrier(MPI_COMM_WORLD, my_mpi_err)
  call fft_grid_data_global%destroyPlans()
  call MPI_Barrier(MPI_COMM_WORLD, my_mpi_err)
  call h5close_f(ierr)
  call MPI_Barrier(MPI_COMM_WORLD, my_mpi_err)
  call finalize_mpi()
  stop
endif

call sim_all_macro_data%sim_time%advanceTime()
call sim_all_macro_data%sim_time%setDeltaTime(accept_reject_dt)
call sim_all_macro_data%sim_time%updateCurrentTime()

call StressGuessNextOuterLoop()


rve_stress_t(:,:)=rve_stress(:,:)
rve_equivalent_total_strain_rate_older = rve_equivalent_total_strain_rate_old
rve_equivalent_total_strain_rate_old = rve_equivalent_total_strain_rate
rve_stress_vm_old = rve_stress_vm

last_converged_iterations = outer_loop_iteration
END DO ! ON THE SIMULATION TIME

if (output_err.and. AmIMPIMaster()) close(21)

log_message = "SIMULATION COMPLETED, BECAUSE TIME >= SIMUALTION TIME END"
call writeToLogFile(log_message)

call fft_grid_data_global%destroyPlans()
call h5close_f(ierr)
call finalize_mpi()

END PROGRAM EVPFFT
