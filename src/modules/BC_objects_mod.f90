module bc_objects_mod
#include "macro_debug.fpp"
  use kinds
  use read_from_file_utils
  use string_module
  use print_utils_mod, only : printToScreen
  use polymorphic_dtype_array_mod, only : dtype_array_ptr
  use simumaltion_macro_data_mod, only : sim_all_macro_data_obj
  implicit none

  ! boudnary condition variables
  real(k_real), dimension(:,:), pointer :: i_scauchy => null(), &
                                           i_scauchy_rate => null(), &
                                           i_vel_grad => null()
  real(k_real), pointer :: i_tempetature => null(), i_tempetature_rate => null()

  integer, dimension(:,:), pointer :: i_scauchy_components => null(), &
                                      i_vel_grad_comp => null()


  type general_stress_displacement_temperature_rate_bc
    class(sim_all_macro_data_obj), pointer :: sim_all_macro_data => null()

    logical :: is_initialized=.FALSE., &
               use_time_end_criterion=.FALSE.

    real(k_real) :: block_duration, block_time_begin, block_time_end
    real(k_real) :: temperature_rate, temperature_old

    integer, dimension(6) :: imposed_stress_rate_components_voigt
    real(k_real), dimension(6) :: imposed_stress_rate_values_voigt, stress_old_voigt

    integer, dimension(3,3) :: imposed_stress_rate_components
    real(k_real), dimension(3,3) :: imposed_stress_rate_values, stress_old

    integer, dimension(3,3) :: imposed_velocity_gradient_components
    integer, dimension(6) :: imposed_velocity_gradient_components_voigt
    real(k_real), dimension(3,3) :: imposed_velocity_gradient_values, displacement_gradient_old
    type(string_type)  :: block_ends_criterion

    logical :: stress_controlled =.FALSE.
    logical :: displacement_controlled =.FALSE.
    logical :: mixed =.FALSE.

    real(k_real), pointer :: time => null(), &
                             dt => null()
  contains
    procedure :: linkBaseObjects => linkBaseObjectsBCBlockGeneral
    procedure :: readBCBlock => readBCBlockGeneral
    procedure :: checkBC
    procedure :: initializeBlock !must be called at the point in time we first need this block bc
    procedure :: checkIfBlockCompleted
    procedure :: computeCurrentBC => computeCurrentBCRateGeneral
    procedure :: updateBC => updateBCRateGeneral
    procedure :: forceTime => forceTimeBCRateGeneral
  end type

  type, extends(dtype_array_ptr) :: boundary_condition_array_type
    class(sim_all_macro_data_obj), pointer :: sim_all_macro_data => null()
    integer :: idx_last_used_bc = 0

    real(k_real), pointer, dimension(:,:) :: avg_stress_ptr => null(), avg_u_grad_ptr => null()
    real(k_real), pointer, dimension(:,:) :: current_imposed_stress_ptr => null(), &
                                             current_imposed_stress_rate_ptr => null(), &
                                             current_imposed_velgrad_ptr => null(), &
                                             current_imposed_dispgrad_ptr => null()
    integer, pointer, dimension(:,:) :: current_imposed_stress_componets_ptr => null(), &
                                        current_imposed_velgrad_components_ptr => null()

    real(k_real), pointer :: current_imposed_temperature_ptr => null(), &
                             current_imposed_temperature_rate_ptr => null()

    real(k_real) :: time_simulation_begin, temperature_simulation_begin, &
                    stress_simulation_begin(3,3), disp_grad_simulation_begin(3,3)
    class(general_stress_displacement_temperature_rate_bc), pointer :: current_bc_ptr=> null()
    integer :: num_bc_blocks
    logical :: im_allocated = .FALSE.

  contains
    procedure :: linkBaseObjects => linkBaseObjectsBCObject
    procedure :: getTimeSimulationBegins
    procedure :: initgridPointers => initgridPointersBCObject
    procedure :: addElement => addBoundaryConditionBlock
    procedure :: getElementPtr => getBoundaryConditionBlockPtr
    procedure :: readBCFile
    procedure :: computeCurrentBC => computeCurrentBCBCObject
    procedure :: updateBC => updateBCBCObject
    procedure :: checkIfBlockCompleted => checkIfBlockCompletedBCBCObject
    procedure :: initBC => initBCObject
    procedure :: forceTime => forceTimeBCObject
    procedure :: printCurrentImposedBC
    procedure :: getTemperaturePointer
  end type

contains

  subroutine getTimeSimulationBegins(this, time_simulation_begin)
    class(boundary_condition_array_type), intent(in) :: this
    real(k_real), intent(out) :: time_simulation_begin
    
    time_simulation_begin = this%time_simulation_begin
  end subroutine

  subroutine linkBaseObjectsBCBlockGeneral(this, sim_all_macro_data)
    class(general_stress_displacement_temperature_rate_bc), intent(inout) :: this
    class(sim_all_macro_data_obj), intent(in), target :: sim_all_macro_data

    this%sim_all_macro_data => sim_all_macro_data
    call sim_all_macro_data%sim_time%getTimePointer(this%time)
    call sim_all_macro_data%sim_time%getDeltaTimePointer(this%dt)
  end subroutine

  subroutine forceTimeBCRateGeneral(this, force_time, time_to_be_forced)
    class(general_stress_displacement_temperature_rate_bc), intent(in) :: this
    logical, intent(out) :: force_time
    real(k_real), intent(out) :: time_to_be_forced

    force_time =.FALSE.
    if ((this%time .gt. this%block_time_end).and. &
        ((this%time-this%dt).lt.this%block_time_end)) then
      force_time =.TRUE.
      time_to_be_forced = this%block_time_end
    end if
  end subroutine

  subroutine readBCBlockGeneral(this, f_reader, block_id)
    use mpi_useful_routines_mod, only : AmIMPIMaster
    use log_file_mod, only : write_detailed_log_to_screen
    use number_to_string_mod
    implicit none
    class(file_reader), intent(inout) :: f_reader
    class(general_stress_displacement_temperature_rate_bc), intent(inout) :: this
    integer, intent(in) :: block_id
    type(string_array) :: dummy_string_array

    __DECL_UNUSED_INTEGER__

    __SUPPRESS_UNUSED_INTEGER__(block_id)
    call f_reader%readParameter("-Block-ends-criterion", this%block_ends_criterion)

    select case(this%block_ends_criterion%getString())
    case ("time")
      this%use_time_end_criterion = .TRUE.
      call f_reader%readParameter("-Block-duration[s]", this%block_duration)

    case default
      if (AmIMPIMaster().and.write_detailed_log_to_screen) then
        write(*,*) "unrecognized termination criterion ", this%block_ends_criterion%getString()
        write(*,*) "available criterion are: "
        write(*,*) "time "
      endif
      error stop "abort"
    end select
    call f_reader%skipEmptyLine()

    ! read Tempearture BC
    call f_reader%readLineAndCheckStringsAreEqual("-Temperature-BC", dummy_string_array)
    call f_reader%readParameter("Imposed-Temperature-Rate[K/s]", this%temperature_rate)
    call f_reader%skipEmptyLine()

    ! read Stress BC
    call f_reader%readLineAndCheckStringsAreEqual("-Stress-BC", dummy_string_array)
    call f_reader%readLineAndCheckStringsAreEqual("Imposed-Stress-Rate-Components", dummy_string_array)
    call f_reader%readStress6FromFile(this%imposed_stress_rate_components_voigt)
    call f_reader%readLineAndCheckStringsAreEqual("Imposed-Stress-Rate-Values[MPa/s]", dummy_string_array)
    call f_reader%readStress6FromFile(this%imposed_stress_rate_values_voigt)
    call f_reader%skipEmptyLine()

    ! read Displacement BC
    call f_reader%readLineAndCheckStringsAreEqual("-Displacement-BC", dummy_string_array)
    call f_reader%readLineAndCheckStringsAreEqual("Imposed-Velocity-Gradient-Components", dummy_string_array)
    call f_reader%readDisplacement33FromFile(this%imposed_velocity_gradient_components)
    call f_reader%readLineAndCheckStringsAreEqual("Imposed-Velocity-Gradient-Values[1/s]", dummy_string_array)
    call f_reader%readDisplacement33FromFile(this%imposed_velocity_gradient_values)

    call this%checkBC
  end subroutine

  subroutine checkBC(this)
    use mpi_useful_routines_mod, only : AmIMPIMaster
    use voigt_indicial_conversion_mod, only : VectorVoigtToTensor2
    implicit none
    class(general_stress_displacement_temperature_rate_bc), intent(inout) :: this

    ! first transform the stress flag6 into a rank two tensor and do the same for the imposed rate
    this%imposed_stress_rate_components = VectorVoigtToTensor2(this%imposed_stress_rate_components_voigt)
    this%imposed_stress_rate_values = VectorVoigtToTensor2(this%imposed_stress_rate_values_voigt)

    !check we are not overconstraining the system
    if (any((this%imposed_stress_rate_components + this%imposed_velocity_gradient_components).gt.1).and.AmIMPIMaster()) then
      write(*,*) "we cannot impose the same component of the strain rate and of the stress rate!"
      call printToScreen(this%imposed_stress_rate_components, "Imposed stress rate components")
      call printToScreen(this%imposed_velocity_gradient_components, "Imposed velocity gradient components")
      call printToScreen(this%imposed_stress_rate_components+this%imposed_velocity_gradient_components, "sum of boundary conditions")
      error stop "Bounday conditions are illegal, the problem is overconstrained"
    endif

    !check we are not underconstraining the problem
    if (any((this%imposed_stress_rate_components + this%imposed_velocity_gradient_components).lt.1).and.AmIMPIMaster()) then
      write(*,*) "one or more stress or strain rate components are not constrained"
      call printToScreen(this%imposed_stress_rate_components, "Imposed stress rate components")
      call printToScreen(this%imposed_velocity_gradient_components, "Imposed velocity gradient components")
      call printToScreen(this%imposed_stress_rate_components+this%imposed_velocity_gradient_components, "sum of boundary conditions")
      error stop "Bounday conditions are illegal, the problem is underconstrained"
    endif

    ! check we are imposing 2 diagonal components of the velocity gradient
    if ((this%imposed_velocity_gradient_components(1,1)+this%imposed_velocity_gradient_components(2,2)+this%imposed_velocity_gradient_components(3,3).eq.2).and.AmIMPIMaster()) then
      write(*,*) "error on imposing diagonal components of the velocity gradient"
      call printToScreen(this%imposed_velocity_gradient_components, "Imposed velocity gradient components")
      error stop "You cannot impose 2 diagonal components of the velcoit gradient"
    endif

    !check the imposed velocity gradient is symmetric
    if (any(this%imposed_velocity_gradient_values*this%imposed_velocity_gradient_components.ne.transpose(this%imposed_velocity_gradient_values*this%imposed_velocity_gradient_components)).and.AmIMPIMaster()) then
      write(*,*) " the imposed velocity gradient is non symmetric"
      call printToScreen(this%imposed_velocity_gradient_values, "Imposed velocity gradient components")
      call printToScreen(this%imposed_velocity_gradient_values-transpose(this%imposed_velocity_gradient_values), "udot - udot^T")
      error stop
    endif

    ! if we are here we passed all the BC checks
    this%stress_controlled =.FALSE.
    this%displacement_controlled =.FALSE.
    if (sum(this%imposed_stress_rate_components).eq.9) this%stress_controlled =.TRUE.
    if (sum(this%imposed_velocity_gradient_components).eq.9) this%displacement_controlled =.TRUE.
    if ((.not.(this%stress_controlled)).and.(.not.(this%displacement_controlled)) ) this%mixed =.TRUE.

  end subroutine

  subroutine initializeBlock(this, stress_begin, displacement_gradient_begin, temperature_begin, time_begin_in)
    use voigt_indicial_conversion_mod, only : Tensor2ToVectorVoigt
    implicit none
    class(general_stress_displacement_temperature_rate_bc), intent(inout) :: this
    real(k_real), intent(in) :: stress_begin(3,3), displacement_gradient_begin(3,3), temperature_begin
    real(k_real), intent(in), optional :: time_begin_in
    if (present(time_begin_in)) then
      if (time_begin_in.ne.this%time) error stop "initializeBCBlock time_begin != this%time. Abort"
      this%block_time_begin = time_begin_in
    else
      this%block_time_begin = this%time
    endif

    if (this%use_time_end_criterion) this%block_time_end = this%block_time_begin + this%block_duration
    this%stress_old = stress_begin
    this%stress_old_voigt = Tensor2ToVectorVoigt(this%stress_old)
    this%displacement_gradient_old = displacement_gradient_begin
    this%temperature_old = temperature_begin
  end subroutine


  subroutine checkIfBlockCompleted(this, block_completed)
    use voigt_indicial_conversion_mod, only : Tensor2ToVectorVoigt
    implicit none
    class(general_stress_displacement_temperature_rate_bc), intent(inout) :: this
    logical, intent(out) :: block_completed
    real(k_real), parameter :: bc_time_tol = 1e-10
    block_completed = .FALSE.
    if (this%use_time_end_criterion) then
     if (abs(this%time - this%block_time_end).le.bc_time_tol) then
       block_completed=.TRUE.
     endif
    endif
  end subroutine

  subroutine computeCurrentBCRateGeneral(this, stress, stress_rate, vel_grad, u_grad, istress_rate_components, ivel_grad_components, temperature, temperature_rate)
    implicit none
    class(general_stress_displacement_temperature_rate_bc), intent(in) :: this
    real(k_real), intent(out), optional, dimension(3,3) :: stress, &
                                                           stress_rate, &
                                                           vel_grad, &
                                                           u_grad
    real(k_real), intent(out), optional :: temperature, temperature_rate
    integer, intent(out), optional, dimension(3,3) :: istress_rate_components, ivel_grad_components

    real(k_real), dimension(3,3)  :: stress_temp, &
                                     stress_rate_temp, &
                                     vel_grad_temp, &
                                     u_grad_temp
    real(k_real) :: temperature_temp, temperature_rate_temp

    stress_temp = this%stress_old + &
             this%imposed_stress_rate_values *this%dt * this%imposed_stress_rate_components
    stress_rate_temp = this%imposed_stress_rate_values
    vel_grad_temp = this%imposed_velocity_gradient_values*this%imposed_velocity_gradient_components
    u_grad_temp = this%displacement_gradient_old +  &
             this%imposed_velocity_gradient_values*this%dt*this%imposed_velocity_gradient_components
    temperature_temp = this%temperature_old + this%temperature_rate* this%dt
    temperature_rate_temp = this%temperature_rate
    if (present(stress)) stress = stress_temp
    if (present(stress_rate)) stress_rate = stress_rate_temp
    if (present(vel_grad)) vel_grad = vel_grad_temp
    if (present(u_grad)) u_grad = u_grad_temp
    if (present(istress_rate_components)) istress_rate_components = this%imposed_stress_rate_components
    if (present(ivel_grad_components)) ivel_grad_components = this%imposed_velocity_gradient_components
    if (present(temperature)) temperature = temperature_temp
    if (present(temperature_rate)) temperature_rate = temperature_rate_temp

  end subroutine

  subroutine updateBCRateGeneral(this, block_completed)
    implicit none
    class(general_stress_displacement_temperature_rate_bc), intent(inout) :: this
    logical, intent(out) :: block_completed
    real(k_real) :: dummy_tensor2(3,3)

    block_completed = .FALSE.

    call this%computeCurrentBC(this%stress_old, dummy_tensor2, this%displacement_gradient_old, temperature=this%temperature_old)
    call this%checkIfBlockCompleted(block_completed)
  end subroutine

  !------------------------------------------------------------------------------!
  !                  BOUNDARY CONDITION ARRAY SUBROUTINE                         !
  !------------------------------------------------------------------------------!

  subroutine readBCFile(this, bc_filename)
    use number_to_string_mod
    use mpi_useful_routines_mod, only : AmIMPIMaster
    use log_file_mod, only : write_detailed_log_to_screen
    use voigt_indicial_conversion_mod, only : VectorVoigtToTensor2
    implicit none
    class(boundary_condition_array_type), intent(inout) :: this
    class(general_stress_displacement_temperature_rate_bc), pointer :: new_element_general => null()
    type(file_reader) :: f_reader
    character(len=*), intent(in) :: bc_filename
    type(string_array) :: dummy_string_array
    integer :: i, block_id
    type(string_type)  :: block_type
    logical :: new_simulation, finish_incomplete_simulation, start_from_preload
    integer :: start_from_sum
    real(k_real) :: initial_stress(6)

    call f_reader%openReadTextFile(bc_filename)
    call f_reader%readLineAndCheckStringsAreEqual("--Boundary-Condition-File", dummy_string_array)
    call f_reader%readParameter("-new-simulation", new_simulation)
    call f_reader%readParameter("-finish-incomplete-simulation", finish_incomplete_simulation)
    call f_reader%readParameter("-start-from-preload", start_from_preload)
    call f_reader%skipEmptyLine()

    ! sanity checks
    start_from_sum = 0
    if (new_simulation) start_from_sum = start_from_sum +1
    if (finish_incomplete_simulation) start_from_sum = start_from_sum +1
    if (start_from_preload) start_from_sum = start_from_sum +1

    if (start_from_sum>1) error stop "Boundary codnition file. More than one starting option is true. Abort"
    if (start_from_sum.eq.0) error stop "Boundary codnition file. You must select one starting option. Abort"
    
    call f_reader%readLineAndCheckStringsAreEqual("--Initial-State", dummy_string_array)
    call f_reader%readParameter("-Time-Simulation-Begins[s]", this%time_simulation_begin)
    call  f_reader%readParameter("-Initial-Temperature[K]", this%temperature_simulation_begin)
    call f_reader%readLineAndCheckStringsAreEqual("-Initial-Stress[MPa]", dummy_string_array)

    call f_reader%readStress6FromFile(initial_stress)
    this%stress_simulation_begin = VectorVoigtToTensor2(initial_stress)
    call f_reader%readLineAndCheckStringsAreEqual("-Initial-Displacement-Gradient[unitless]", dummy_string_array)
    call f_reader%readDisplacement33FromFile(this%disp_grad_simulation_begin)
    call f_reader%skipEmptyLine()

    ! sanity checks
    if (this%temperature_simulation_begin.lt.0._k_real) error stop   "Boundary codnition file. You cannot input a negative temeprature!"
    if (new_simulation) then
      if (any(this%stress_simulation_begin.ne.0._k_real)) error stop &
        "Boundary codnition file. You cannot impose any non zero initial stress component for a new simulation"
      if (any(this%disp_grad_simulation_begin.ne.0._k_real)) error stop &
        "Boundary codnition file. You cannot impose any non zero initial displacement gradient component for a new simulation"
    endif

    call f_reader%readParameter("--Number-of-BC-Blocks", this%num_bc_blocks)
    call f_reader%skipEmptyLine()

    do i = 1, this%num_bc_blocks
      call f_reader%readParameter("--Block", block_id)
      if (block_id.ne.i) then
        write(*,*) "readBCFile the boundary condition block id in the file is "//&
        trim(adjustl(int2string(block_id)))// &
        ", instead i was expecting block number "//trim(adjustl(int2string(i)))
        error stop "error in the boudnary condition file"
      endif
      call f_reader%readParameter("-Block-Type", block_type)
      select case(block_type%getString())
      case ("rate-general")
        nullify(new_element_general)
        allocate(new_element_general)
        call new_element_general%readBCBlock(f_reader, block_id)
        call this%addElement(new_element_general)
        nullify(new_element_general)
      case default
        if (AmIMPIMaster().and.write_detailed_log_to_screen) write(*,*) "unrecognized boundary condition block type ", block_type%getString()
        if (AmIMPIMaster().and.write_detailed_log_to_screen) write(*,*) "available block types are: "
        if (AmIMPIMaster().and.write_detailed_log_to_screen) write(*,*) "rate-general "
        error stop "abort"
      end select
      if (i.ne.this%num_bc_blocks) call f_reader%skipEmptyLine()
    enddo
    call f_reader%closeTextFile()
  end subroutine

  subroutine addBoundaryConditionBlock(this, new_element)
    class(boundary_condition_array_type), intent(inout) :: this
    class(general_stress_displacement_temperature_rate_bc), pointer, intent(inout) :: new_element

    call this%extend()
    this%all_pt(this%n_elements)%pt => new_element
    nullify(new_element)
  end subroutine

  subroutine getBoundaryConditionBlockPtr(this, idx, ptr)
    class(boundary_condition_array_type), intent(in) :: this
    integer, intent(in) :: idx
    class(general_stress_displacement_temperature_rate_bc), pointer, intent(out) :: ptr

    call this%checkElementExist(idx)

    associate(elem => this%all_pt(idx)%pt)
    select type(elem)
      class is (general_stress_displacement_temperature_rate_bc)
        ptr => elem
      class default
        error stop "getBoundaryConditionBlockPtr wrong class"
    end select
    end associate
  end subroutine

  subroutine linkBaseObjectsBCObject(this, sim_all_macro_data)
    class(boundary_condition_array_type), intent(inout) :: this
    class(sim_all_macro_data_obj), intent(inout), target :: sim_all_macro_data
    class(general_stress_displacement_temperature_rate_bc), pointer :: bc_ptr => null()
    integer :: i
    this%sim_all_macro_data => sim_all_macro_data

    ! allocate pointers
    allocate(this%current_imposed_stress_ptr(3,3), &
             this%current_imposed_stress_rate_ptr(3,3), &
             this%current_imposed_velgrad_ptr(3,3), &
             this%current_imposed_dispgrad_ptr(3,3), &
             this%current_imposed_stress_componets_ptr(3,3), &
             this%current_imposed_velgrad_components_ptr(3,3), &
             this%current_imposed_temperature_ptr, &
             this%current_imposed_temperature_rate_ptr )

    i_scauchy => this%current_imposed_stress_ptr
    i_scauchy_rate => this%current_imposed_stress_rate_ptr
    i_vel_grad => this%current_imposed_velgrad_ptr
    i_scauchy_components => this%current_imposed_stress_componets_ptr
    i_vel_grad_comp => this%current_imposed_velgrad_components_ptr
    i_tempetature => this%current_imposed_temperature_ptr
    i_tempetature_rate => this%current_imposed_temperature_rate_ptr
    do i = 1, this%num_bc_blocks
      nullify(bc_ptr)
      call this%getElementPtr(i, bc_ptr)
      call bc_ptr%linkBaseObjects(sim_all_macro_data)
    enddo

    this%im_allocated=.true.

    call sim_all_macro_data%sim_macro_field_averages%setTemperaturePointer(this%current_imposed_temperature_ptr)
  end subroutine

  subroutine initGridPointersBCObject(this)
    class(boundary_condition_array_type), intent(inout) :: this

    call this%sim_all_macro_data%sim_macro_field_averages%getAverageStressPointer(this%avg_stress_ptr)
    call this%sim_all_macro_data%sim_macro_field_averages%getAverageDisplacementGradientPointer(this%avg_u_grad_ptr)


  end subroutine



  subroutine computeCurrentBCBCObject(this, stress, stress_rate, vel_grad, u_grad, istress_rate_components, ivel_grad_components, temperature, temperature_rate)
    implicit none
    class(boundary_condition_array_type), intent(inout) :: this
    real(k_real), intent(out), optional, dimension(3,3) :: stress, &
                                                           stress_rate, &
                                                           vel_grad, &
                                                           u_grad
    real(k_real), intent(out), optional :: temperature, temperature_rate
    integer, intent(out), optional :: istress_rate_components(3,3), ivel_grad_components(3,3)
    real(k_real),  dimension(3,3) :: stress_temp, &
                                     stress_rate_temp, &
                                     vel_grad_temp, &
                                     u_grad_temp
    real(k_real) :: temperature_temp, temperature_rate_temp

    integer, dimension(3,3) :: istress_rate_components_temp, ivel_grad_components_temp

    if (.not.associated(this%current_bc_ptr)) error stop "computeCurrentBCBCObject this%current_bc_ptr not associated"
    call this%current_bc_ptr%computeCurrentBC(stress_temp, stress_rate_temp, vel_grad_temp, u_grad_temp, istress_rate_components_temp, ivel_grad_components_temp, &
                                              temperature_temp, temperature_rate_temp)

    this%current_imposed_velgrad_ptr = vel_grad_temp
    this%current_imposed_stress_ptr = stress_temp
    this%current_imposed_stress_rate_ptr = stress_rate_temp
    this%current_imposed_velgrad_components_ptr = ivel_grad_components_temp
    this%current_imposed_stress_componets_ptr = istress_rate_components_temp
    this%current_imposed_temperature_ptr = temperature_temp
    this%current_imposed_temperature_rate_ptr = temperature_rate_temp

    if (present(stress)) stress = stress_temp
    if (present(stress_rate)) stress_rate = stress_temp
    if (present(vel_grad)) vel_grad = vel_grad_temp
    if (present(u_grad)) u_grad = u_grad_temp
    if (present(istress_rate_components)) istress_rate_components = istress_rate_components_temp
    if (present(ivel_grad_components)) ivel_grad_components = ivel_grad_components_temp
    if (present(temperature)) temperature = temperature_temp
    if (present(temperature_rate)) temperature_rate = temperature_rate_temp
  end subroutine

  subroutine updateBCBCObject(this, terminate_simulation, terminate_message, new_bc_block)
    use log_file_mod
    use number_to_string_mod
    implicit none
    class(boundary_condition_array_type), intent(inout) :: this
    logical, intent(inout) :: terminate_simulation, new_bc_block
    character(len=*), intent(inout) :: terminate_message
    logical :: block_completed
    real(k_real), dimension(3,3) :: stress_prev, u_grad_prev, &
                                    new_bloc_stress_begin, new_block_u_grad_begin
    integer, dimension(3,3) :: i_stress_rate_comp_prev, i_vel_grad_comp_prev
    real(k_real) :: temperature_prev

    integer :: i,j
    new_bc_block = .false.
    call this%current_bc_ptr%checkIfBlockCompleted(block_completed)
    if (block_completed) then
      if (this%idx_last_used_bc.eq.this%num_bc_blocks) then
        if (.not.(terminate_simulation)) then
          terminate_simulation = .TRUE.
          terminate_message="SIMULATION COMPLETED, ALL BOUNDARY CONDITIONS BLOCKS HAVE BEEN EXECUTED"
        else
          terminate_message=terminate_message//NEW_LINE('A')//"SIMULATION COMPLETED, ALL BOUNDARY CONDITIONS BLOCKS HAVE BEEN EXECUTED"
        endif
      else

        log_message = "BOUNDARY CONDITION BLOCK " // trim(adjustl(int2string(this%idx_last_used_bc))) &
                // "COMPLETED "
        call writeToLogFile(log_message)

        !get last BCs info from previous block
        call this%current_bc_ptr%computeCurrentBC( stress = stress_prev, &
                                                   u_grad= u_grad_prev, &
                                                   istress_rate_components = i_stress_rate_comp_prev, &
                                                   ivel_grad_components = i_vel_grad_comp_prev, &
                                                   temperature = temperature_prev)
        ! the new stress begin must be equal to the imposed bc for the imposed components,
        ! and equal to the the average values for the non imposed components
        new_bloc_stress_begin = stress_prev
        new_block_u_grad_begin = u_grad_prev
        do j=1,3
          do i=1,3
            if (i_stress_rate_comp_prev(i,j).eq.0) new_bloc_stress_begin(i,j) = this%avg_stress_ptr(i,j)
            if (i_vel_grad_comp_prev(i,j).eq.0) new_block_u_grad_begin(i,j) = this%avg_u_grad_ptr(i,j)
          enddo
        enddo

        ! actually move to the next block
        this%idx_last_used_bc = this%idx_last_used_bc + 1
        call this%getElementPtr(this%idx_last_used_bc, this%current_bc_ptr)
        call this%current_bc_ptr%initializeBlock(new_bloc_stress_begin, new_block_u_grad_begin, temperature_prev)
        log_message = "NOW MOVING TO BOUNDARY CONDITION BLOCK " // trim(adjustl(int2string(this%idx_last_used_bc)))
        call writeToLogFile(log_message)
        new_bc_block = .true.
      endif
    else
      call this%current_bc_ptr%updateBC(block_completed)
    endif
  end subroutine

  subroutine checkIfBlockCompletedBCBCObject(this, block_completed)
    class(boundary_condition_array_type), intent(in) :: this
    logical, intent(out) :: block_completed

    block_completed = .FALSE.
    call this%current_bc_ptr%checkIfBlockCompleted(block_completed) 
  end subroutine

  subroutine forceTimeBCObject(this)
    use number_to_string_mod
    use log_file_mod
    implicit none
    class(boundary_condition_array_type), intent(in) :: this
    logical :: force_time
    real(k_real) :: time_to_be_forced

    call this%current_bc_ptr%forceTime(force_time, time_to_be_forced)
    if (force_time) then
      call this%sim_all_macro_data%sim_time%ForceTime(time_to_be_forced)
      log_message = "FORCING END OF BOUDANRY CONDITION BLOCK " // trim(adjustl(int2string(this%idx_last_used_bc))) &
              // ": END BLOCK TIME IS " // trim(adjustl(real2string(time_to_be_forced)))
      call writeToLogFile(log_message)
    endif
  end subroutine

  subroutine initBCObject(this)
    implicit none
    class(boundary_condition_array_type), intent(inout) :: this

    this%idx_last_used_bc = 1
    call this%getElementPtr(this%idx_last_used_bc, this%current_bc_ptr)
    call this%current_bc_ptr%initializeBlock(this%stress_simulation_begin,&
                                             this%disp_grad_simulation_begin, &
                                             this%temperature_simulation_begin, &
                                             this%time_simulation_begin)
    call this%computeCurrentBC()
  end subroutine

  subroutine printCurrentImposedBC(this)
    use log_file_mod
    use number_to_string_mod
    use log_file_mod
    use print_utils_mod, only : writeToString
    implicit none
    class(boundary_condition_array_type), intent(in) :: this


    log_message = "-------------------------------------------------------------"
    log_message = trim(adjustl(log_message)) //new_line("A") // " IMPOSED BOUDNARY CONDITIONS "
    log_message = trim(adjustl(log_message)) //new_line("A") //  "-------------------------------------------------------------"
    call writeToLogFile(log_message)

    call writeToString(this%current_imposed_temperature_rate_ptr, "Imposed Temperature Rate", log_message)
    call writeToLogFile(log_message)

    call writeToString(this%current_imposed_temperature_ptr, "Imposed Temperature", log_message)
    call writeToLogFile(log_message)

    call writeToString(this%current_imposed_stress_rate_ptr, "Imposed Stress Rate", log_message)
    call writeToLogFile(log_message)

    call writeToString(this%current_imposed_stress_componets_ptr, "Imposed Stress Rate Components", log_message)
    call writeToLogFile(log_message)

    call writeToString(this%current_imposed_stress_ptr, "Imposed Stress", log_message)
    call writeToLogFile(log_message)

    call writeToString(this%current_imposed_velgrad_ptr, "Imposed Velocity Gradient", log_message)
    call writeToLogFile(log_message)

    call writeToString(this%current_imposed_velgrad_components_ptr, "Imposed Velocity Gradient Components", log_message)
    call writeToLogFile(log_message)

    log_message = "-------------------------------------------------------------"
    log_message = trim(adjustl(log_message)) //new_line("A") // " END IMPOSED BOUDNARY CONDITIONS "
    log_message = trim(adjustl(log_message)) //new_line("A") //  "-------------------------------------------------------------"
    call writeToLogFile(log_message)

  end subroutine

  subroutine getTemperaturePointer(this, temperature_ptr)
    class(boundary_condition_array_type), intent(in) :: this
    real(k_real), intent(out), pointer :: temperature_ptr

    if (.not.(this%im_allocated)) error stop "getTemperaturePointer not intialized. Abort"
    if (associated(temperature_ptr)) error stop "getTemperaturePointer temperature_ptr already associated. Abort!"
    if (.not.associated(this%current_imposed_temperature_ptr)) error stop "getTemperaturePointer ptr not associated. Abort!"

    temperature_ptr => this%current_imposed_temperature_ptr
  end subroutine

end module bc_objects_mod
