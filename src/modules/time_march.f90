module time_march_mod
  use kinds
  use log_file_mod
  use global, only : total_strain, total_strain_old, &
                     total_strain_rate, total_strain_rate_old, &
                     grid_stress, grid_stress_old, &
                     grid_stress_rate, grid_stress_rate_old, &
                     tdot, tdot_max, &
                     phase_material_array, & 
                     all_mighty_grid_global
  use mpi_useful_routines_mod, only : MPIMaxIncrementGridTensor2, MPIComputeQualityTensor2
  use mpi_variables_mod, only : mpi_master_rank, mpi_rank, i_am_mpi_master
implicit none
real(k_real), pointer, dimension(:) :: stress_abs_rel_tol => null()
real(k_real), pointer, dimension(:) :: total_strain_abs_rel_tol => null()

real(k_real) :: max_allowed_stress_increment = 10._k_real
real(k_real), parameter :: time_step_safety_factor = 0.9_k_real

contains
  subroutine readGlobalTimemarchTolerances(optf_reader)
    use string_module
    use read_from_file_utils
    class(file_reader), intent(inout) :: optf_reader
    type(string_array) :: dummy_string_array
    call optf_reader%readLineAndCheckStringsAreEqual("--Time-March-Tolerances", dummy_string_array)
    call optf_reader%readVectorParameter("stress-absolute-and-relative-tol[MPa,unitless]", 2, stress_abs_rel_tol)
    call optf_reader%readVectorParameter("total-strain-absolute-and-relative-tol[mm/mm,unitless]", 2, total_strain_abs_rel_tol)

  end subroutine

  subroutine acceptRejectSolution(accept, new_dt)
    use number_to_string_mod
    use log_file_mod, only : writeToScreen, write_detailed_log_to_screen
    implicit none
    logical, intent(out) :: accept
    logical :: accept_var, accept_reject_global_vars
    real(k_real), intent(out) :: new_dt
    real(k_real) :: max_dt_var
    real(k_real) :: max_total_strain_increment
    real(k_real) :: max_stress_increment
    real(k_real) :: Q

    log_message = ""
    accept = .true.
    new_dt = tdot_max



    call phase_material_array%acceptRejectSolution(new_dt, accept)

    accept_reject_global_vars = .true.
    call writeToScreen("***************************************************************")
    call writeToScreen("*******GLOBAL STATE VARIABLES ACCEPT/REJECT SOLUTION***********")
    call writeToScreen("")


    ! check total strain incremenet
    call MPIComputeQualityTensor2(total_strain, total_strain_old, total_strain_rate, total_strain_rate_old, &
    all_mighty_grid_global%AMGgetNumVoxel()*9, tdot, total_strain_abs_rel_tol(2), total_strain_abs_rel_tol(1), Q, max_total_strain_increment)
    call computeMaxAllowedTimeStepQuality(Q, tdot, "total strain ",  accept_var, max_total_strain_increment, max_dt_var)
    accept_reject_global_vars = accept_reject_global_vars.and.accept_var
    new_dt = min(new_dt, max_dt_var)


    call MPIComputeQualityTensor2(grid_stress, grid_stress_old, grid_stress_rate, grid_stress_rate_old, &
    all_mighty_grid_global%AMGgetNumVoxel()*9, tdot, stress_abs_rel_tol(2), stress_abs_rel_tol(1), Q, max_stress_increment)
    call computeMaxAllowedTimeStepQuality(Q, tdot, "stress ",  accept_var, max_stress_increment, max_dt_var)
    accept_reject_global_vars = accept_reject_global_vars.and.accept_var
    new_dt = min(new_dt, max_dt_var)

    if (accept_reject_global_vars) call writeToScreen("solution might be ACCEPTED")
    if (.not.accept_reject_global_vars)  call writeToScreen("solution will be REJECTED")
    if (i_am_mpi_master.and.write_detailed_log_to_screen) write(*,*) "new allowable dt is ", new_dt
    call writeToScreen("")


    call writeToScreen("")
    call writeToScreen("*******GLOBAL STATE VARIABLES ACCEPT/REJECT SOLUTION END*******")
    call writeToScreen("***************************************************************")
    call writeToScreen("")

    accept = accept.and.accept_reject_global_vars
    ! here we give a pratival restriction: the new time step cannot be larger than twice the previous one
    new_dt = min(new_dt, tdot*2_k_real)


    call writeToScreen("")
    call writeToScreen("**************************************************************************************************************")
    call writeToScreen("**************************************************************************************************************")
    if(accept.and.(i_am_mpi_master.and.write_detailed_log_to_screen)) then
      if(i_am_mpi_master.and.write_detailed_log_to_screen) write(*,*) "solution has been ACCEPTED, an the new time step is ", new_dt
    elseif((.not.accept).and.(i_am_mpi_master.and.write_detailed_log_to_screen)) then
      if(i_am_mpi_master.and.write_detailed_log_to_screen) write(*,*) "solution has been REJECTED. We will loop back with a time step of ", new_dt
    endif
    call writeToScreen("**************************************************************************************************************")
    call writeToScreen("**************************************************************************************************************")
    call writeToScreen("")

    if(.not.accept) then
    log_message = "solution has been REJECTED because: " // trim(adjustl(log_message))
    call writeToLogFile(log_message)
  endif
  end subroutine

  subroutine computeMaxAllowedTimeStepLinear(current_increment_in, allowed_increment, current_tdot, var_name, accept_increment, max_tdot)
    use log_file_mod, only : writeToScreen, write_detailed_log_to_screen
    implicit none
    real(k_real), intent(in) :: current_increment_in, allowed_increment, current_tdot
    character(len=*), intent(in) :: var_name
    real(k_real), intent(out) :: max_tdot

    logical, intent(out) :: accept_increment

    real(k_real) :: current_increment

    current_increment = max(1e-10_k_real, current_increment_in)
    max_tdot = current_tdot*allowed_increment/current_increment * time_step_safety_factor

    accept_increment = current_increment.le.allowed_increment

    if (accept_increment) then
      if(i_am_mpi_master.and.write_detailed_log_to_screen) write(*,*) "the current increment for ", var_name , " is ", current_increment, " and is smaller then the max allowed increment ", allowed_increment
    else
      if(i_am_mpi_master.and.write_detailed_log_to_screen) write(*,*) "the current increment for ", var_name , " is ", current_increment, " and exceeds the maximum allowed increment ", allowed_increment
    endif
    if(i_am_mpi_master.and.write_detailed_log_to_screen) write(*,*) "the new max allowable time step for ", var_name , " is ", max_tdot
    if(i_am_mpi_master.and.write_detailed_log_to_screen) write(*,*) " "

  end subroutine

  subroutine computeMaxAllowedTimeStepQuality(Q, current_tdot, var_name, accept_increment, max_increment, max_tdot)
    use log_file_mod, only : log_message, writeToScreen, write_detailed_log_to_screen
    use number_to_string_mod
    implicit none
    real(k_real), intent(in) :: Q, current_tdot, max_increment
    character(len=*), intent(in) :: var_name
    real(k_real), intent(out) :: max_tdot

    logical, intent(out) :: accept_increment

    if (Q.gt.0._k_real) then
      max_tdot = current_tdot/sqrt(Q)
    else 
      max_tdot = 1e10_k_real
    endif

    accept_increment = current_tdot<1.1_k_real*max_tdot

    if (accept_increment) then
      max_tdot = max_tdot * time_step_safety_factor
    else 
      max_tdot = max_tdot/3._k_real
    endif
    
    if (accept_increment) then
      if(i_am_mpi_master.and.write_detailed_log_to_screen) write(*,*) "the current quality for ", var_name , " is ", Q, " and is smaller then 1 "
    else
      if(i_am_mpi_master.and.write_detailed_log_to_screen) write(*,*) "the current quality for ", var_name , " is ", Q, " and exceeds the maximum allowed increment of 1"
      if(i_am_mpi_master.and.write_detailed_log_to_screen) write(*,*) "the current increment for ", var_name , " is ", max_increment
    endif
    if(i_am_mpi_master.and.write_detailed_log_to_screen) write(*,*) "the new max allowable time step for ", var_name , " is ", max_tdot
    if(i_am_mpi_master.and.write_detailed_log_to_screen) write(*,*) " "

    if (.not.accept_increment) log_message = NEW_LINE('A')//"* the maximum "// var_name // " increment is: "// &
    trim(adjustl(real2string(max_increment))) // &
    trim(adjustl(log_message))

  end subroutine
end module
