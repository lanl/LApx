module termination_criterion_mod
  use kinds
implicit none
logical :: terminate_simulation_global
character(len=1000) :: termination_message_global

logical :: use_VMeq_strain_rate_increment_termination = .FALSE.
real(k_real) :: VMeq_strain_rate_increment_termination_percentage
real(k_real) :: min_VM_eq_strain_rate = 1e6_k_real
real(k_real) :: increase_from_minimum_creep_rate_check_after_time

logical :: use_min_allowed_eq_creep_rate_for_termination = .FALSE.
real(k_real) :: min_allowed_creep_rate
real(k_real) :: current_creep_rate
real(k_real) :: minimum_creep_rate_check_after_time

logical :: use_max_allowed_local_porosity = .FALSE.
real(k_real) :: maximum_allowed_local_porosity

logical :: use_maximum_accumualted_plastic_strain_criterion = .FALSE.
real(k_real) :: max_allowed_accumualted_plastic_strain

logical :: use_stress_percentage_drop_criterion = .FALSE.
real(k_real) :: allowed_stress_percentage_drop
integer :: stress_component_i, stress_component_j
real(k_real) :: max_recorded_stress = 0._k_real

contains
  subroutine readTerminationCriteriaOptions(optf_reader)
    use string_module
    use read_from_file_utils
    implicit none
    class(file_reader), intent(inout) :: optf_reader
    type(string_array) :: dummy_string_array

    call optf_reader%readLineAndCheckStringsAreEqual("--Creep-simulation-termination-options", dummy_string_array)

    call optf_reader%readParameter("-use-increase-from-minimum-creep-rate-criterion[TRUE/FALSE]", use_VMeq_strain_rate_increment_termination)
    call optf_reader%readParameter("percent-increase-from-minimum-creep-rate[%]", VMeq_strain_rate_increment_termination_percentage)
    call optf_reader%readParameter("check-after-time[s]", increase_from_minimum_creep_rate_check_after_time)

    if (use_VMeq_strain_rate_increment_termination) VMeq_strain_rate_increment_termination_percentage=VMeq_strain_rate_increment_termination_percentage/100._k_real

    call optf_reader%readParameter("-use-minimum-allowed-creep-rate-criterion[TRUE/FALSE]", use_min_allowed_eq_creep_rate_for_termination)
    call optf_reader%readParameter("minimum-allowed-creep-rate[1/s]", min_allowed_creep_rate)
    call optf_reader%readParameter("check-after-time[s]", minimum_creep_rate_check_after_time)

    call optf_reader%readParameter("-use-maximum-accumualted-plastic-strain-criterion[TRUE/FALSE]", use_maximum_accumualted_plastic_strain_criterion)
    call optf_reader%readParameter("maximum-allowed-accumualted-plastic-strain[%]", max_allowed_accumualted_plastic_strain)
    if (use_maximum_accumualted_plastic_strain_criterion) max_allowed_accumualted_plastic_strain=max_allowed_accumualted_plastic_strain/100._k_real

    call optf_reader%readParameter("-use-max-local-porosity-criterion[TRUE/FALSE]", use_max_allowed_local_porosity)
    call optf_reader%readParameter("max-local-porosity[unitless]", maximum_allowed_local_porosity)

    call optf_reader%readParameter("-use-stress-drop-percentage-criterion[TRUE/FALSE]", use_stress_percentage_drop_criterion)
    call optf_reader%readParameter("percentage-of-stress-drop", allowed_stress_percentage_drop)
    call optf_reader%readParameter("stress-component-i", stress_component_i)
    call optf_reader%readParameter("stress-component-j", stress_component_j)

  end subroutine

  subroutine checkSimulationTermination()
    use bc_objects_mod, only : i_scauchy_components
    use log_file_mod, only : writeToLogFile
    implicit none
    terminate_simulation_global = .FALSE.
    termination_message_global = ""

    ! creep termination criteria
    if(sum(i_scauchy_components).eq.9) then
      if (.not.(terminate_simulation_global).and.use_VMeq_strain_rate_increment_termination) &
        call checkMinimumStrainRateIncremenet(terminate_simulation_global, termination_message_global)

      if (.not.(terminate_simulation_global).and.use_min_allowed_eq_creep_rate_for_termination) &
        call checkMinimumAllowedStrainRate(terminate_simulation_global, termination_message_global)

    endif

    if (.not.(terminate_simulation_global).and.use_maximum_accumualted_plastic_strain_criterion) &
      call checkMaximumAllowedAccumulatedStrain(terminate_simulation_global, termination_message_global)

    ! if (.not.(terminate_simulation_global).and.use_max_allowed_local_porosity) &
    !   call checkMaximumLocalPorosity(terminate_simulation_global, termination_message_global)

    if (.not.(terminate_simulation_global).and.use_stress_percentage_drop_criterion) &
      call checkStressPercentageDrop(terminate_simulation_global, termination_message_global)

    if (terminate_simulation_global) &
      call writeToLogFile(termination_message_global)
    
  end subroutine

  subroutine checkMinimumStrainRateIncremenet(terminate_simualtion, termination_message)
    use global, only : sim_all_macro_data
    use tensor_math_mod, only : computeVMEquivalentStrain
    use number_to_string_mod
    implicit none
    real(k_real) , dimension(:,:), pointer :: total_strain_rate=> null()
    logical, intent(out) :: terminate_simualtion
    character(len=*), intent(out) :: termination_message
    real(k_real) :: eq_tot_strain_rate
    real(k_real), pointer :: time_ptr => null()

    nullify(time_ptr, total_strain_rate)
    call sim_all_macro_data%sim_time%getTimePointer(time_ptr)
    call sim_all_macro_data%sim_macro_field_averages%getAverageStrainRatePointer(total_strain_rate)

    terminate_simualtion = .FALSE.
    termination_message = ""

    if(time_ptr.le.increase_from_minimum_creep_rate_check_after_time) then
      min_VM_eq_strain_rate = 1e6_k_real

    else
      eq_tot_strain_rate = computeVMEquivalentStrain(total_strain_rate)
      write(*,*) "min_VM_eq_strain_rate ", min_VM_eq_strain_rate
      write(*,*) "eq_tot_strain_rate ", eq_tot_strain_rate
      if (eq_tot_strain_rate.lt.min_VM_eq_strain_rate) then
        min_VM_eq_strain_rate = eq_tot_strain_rate
      else
        write(*,*) "(eq_tot_strain_rate-min_VM_eq_strain_rate)/min_VM_eq_strain_rate ", (eq_tot_strain_rate-min_VM_eq_strain_rate)/min_VM_eq_strain_rate
        if ((eq_tot_strain_rate-min_VM_eq_strain_rate)/min_VM_eq_strain_rate.ge.VMeq_strain_rate_increment_termination_percentage) then
          terminate_simualtion = .TRUE.
          termination_message="SIMULATION COMPLETED, THE EQUIVALENT STRAIN RATE IS " // &
          trim(adjustl(real2string(VMeq_strain_rate_increment_termination_percentage))) // &
           " TIMES LARGER THAN THE MINIMUM CREEP RATE" // &
           trim(adjustl(real2string(min_VM_eq_strain_rate)))
        endif
      endif
    endif

  end subroutine

  subroutine checkMinimumAllowedStrainRate(terminate_simualtion, termination_message)
    use global, only : sim_all_macro_data, total_strain_rate
    use tensor_math_mod, only : computeVMEquivalentStrain
    use mpi_useful_routines_mod, only : MPIAverageGridMatrix
    use number_to_string_mod
    implicit none
    logical, intent(out) :: terminate_simualtion
    character(len=*), intent(out) :: termination_message
    real(k_real) :: eq_tot_strain_rate, avg_tot_strain_rate(3,3)

    real(k_real), pointer :: time_ptr => null()

    nullify(time_ptr)
    call sim_all_macro_data%sim_time%getTimePointer(time_ptr)
    call MPIAverageGridMatrix(total_strain_rate, avg_tot_strain_rate)
    terminate_simualtion = .FALSE.
    termination_message = ""
    eq_tot_strain_rate = computeVMEquivalentStrain(avg_tot_strain_rate)
    if(time_ptr.gt.minimum_creep_rate_check_after_time) then
        if (eq_tot_strain_rate.le.min_allowed_creep_rate) then
          terminate_simualtion = .TRUE.
          termination_message="SIMULATION COMPLETED, EQUIVALENT STRAIN RATE "// &
          trim(adjustl(real2string(eq_tot_strain_rate))) // &
           "SMALLER THAN THE MINIMUM ALLOWED STRAIN RATE "// &
           trim(adjustl(real2string(min_allowed_creep_rate)))
        endif
    endif

  end subroutine

  subroutine checkMaximumAllowedAccumulatedStrain(terminate_simualtion, termination_message)
    use global, only : ept
    use tensor_math_mod, only : computeVMEquivalentStrain
    use mpi_useful_routines_mod, only : MPIAverageGridMatrix
    use number_to_string_mod
    implicit none
    logical, intent(out) :: terminate_simualtion
    character(len=*), intent(out) :: termination_message
    real(k_real) :: eq_pl_strain, avg_plastic_strain(3,3)

    call MPIAverageGridMatrix(ept, avg_plastic_strain)

    terminate_simualtion = .FALSE.
    termination_message = ""
    eq_pl_strain = computeVMEquivalentStrain(avg_plastic_strain)
    if (eq_pl_strain.ge.max_allowed_accumualted_plastic_strain) then
      terminate_simualtion = .TRUE.
      termination_message="SIMULATION COMPLETED, THE EQUIVALENT PLASTIC STRAIN"// &
      trim(adjustl(real2string(eq_pl_strain))) // &
       "EXCEEDS THE ALLOWED MAXIMUM "// &
       trim(adjustl(real2string(max_allowed_accumualted_plastic_strain)))
    endif

  end subroutine

  ! subroutine checkMaximumLocalPorosity(terminate_simualtion, termination_message)
  !   use global, only : effective_porosity
  !   use number_to_string_mod
  !   use mpi_useful_routines_mod, only : MPIMaxGridScalar
  !   implicit none
  !   logical, intent(out) :: terminate_simualtion
  !   character(len=*), intent(out) :: termination_message
  !   real(k_real) :: max_local_porosity

  !   terminate_simualtion = .FALSE.
  !   termination_message = ""
  !   call MPIMaxGridScalar(effective_porosity, max_local_porosity)
  !   if (max_local_porosity.ge.maximum_allowed_local_porosity) then
  !     terminate_simualtion = .TRUE.
  !     termination_message="SIMULATION COMPLETED, MAXIMUM LOCAL POROSITY "// &
  !     trim(adjustl(real2string(max_local_porosity))) // &
  !      "EXCEEDS THE ALLOWED MAXIMUM "// &
  !      trim(adjustl(real2string(maximum_allowed_local_porosity)))
  !   endif

  ! end subroutine

  subroutine checkStressPercentageDrop(terminate_simualtion, termination_message)
    use global, only : grid_stress
    use tensor_math_mod, only : computeVMEquivalentStrain
    use mpi_useful_routines_mod, only : MPIAverageGridMatrix
    use number_to_string_mod
    implicit none
    logical, intent(out) :: terminate_simualtion
    character(len=*), intent(out) :: termination_message
    real(k_real) :: current_stress_tensor(3,3), current_stress, current_stress_drop_pct

    terminate_simualtion = .FALSE.
    termination_message = ""

    call MPIAverageGridMatrix(grid_stress, current_stress_tensor)
    current_stress = current_stress_tensor(stress_component_i, stress_component_j)

    if (current_stress.ge.max_recorded_stress) then
      max_recorded_stress = current_stress
    else 
      if (max_recorded_stress.ne.0._k_real) then
        current_stress_drop_pct = (1._k_real-current_stress/max_recorded_stress)*100._k_real
        if (current_stress_drop_pct.gt.allowed_stress_percentage_drop) then
          terminate_simualtion = .TRUE.
          termination_message="SIMULATION COMPLETED, STRESS DROP "// &
          trim(adjustl(real2string(current_stress_drop_pct))) // &
          "GRATER THAN THE ALLOWED STRESS DROP "// &
          trim(adjustl(real2string(allowed_stress_percentage_drop)))
        endif
      endif
    endif

  end subroutine
end module
