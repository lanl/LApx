module simumaltion_macro_data_mod
use all_grid_data_mod, only : all_grid_data
use grid_data_types, only : griddata_tensor2
use kinds
use csv_writer_mod, only : csv_writer
implicit none

type sim_time_obj
  real(k_real), private, pointer :: time=> null(), &
                                    time_old=> null(), &
                                    dt=> null(), &
                                    dt0=> null(), &
                                    dt_old => null(), &
                                    dt_reference => null(), &
                                    dtmin=> null(), &
                                    dtmax=> null()

  logical :: initialized =.false.
  logical :: im_allocated =.false.
  logical :: forced_time =.false.
contains
  procedure :: alloc => allocTimeObject
  procedure :: init => initTimeObject
  procedure :: getTimePointer
  procedure :: getTimeOldPointer
  procedure :: getDeltaTimePointer
  procedure :: getDeltaTimeOldPointer
  procedure :: getDeltaTimeReferencePointer
  procedure :: getDeltaTMinPointer
  procedure :: getDeltaTMaxPointer
  procedure :: resetDtToDt0
  procedure :: reduceDeltaTime
  procedure :: setDeltaTime
  procedure :: setTimeFromReload
  procedure :: updateCurrentTime
  procedure :: increaseDeltaTime
  procedure :: advanceTime
  procedure :: ForceTime
end type

type sim_average_quantities_obj
  class(all_grid_data), pointer :: grid_data => null()
  class(griddata_tensor2), pointer :: vel_grad => null(), &
                                      dis_grad_tot => null(), &
                                      plastic_strain => null(), &
                                      plastic_strain_rate => null(), &
                                      stress => null()

  real(k_real), dimension(:,:), pointer :: vel_grad_avg => null(), &
                                           dis_grad_tot_avg => null(), &
                                           plastic_strain_avg => null(), &
                                           plastic_strain_rate_avg => null(), &
                                           stress_avg => null()

  real(k_real), dimension(:,:), pointer :: strain_avg => null(), &
                                           strain_rate_avg => null(), &
                                           elastic_strain_avg => null(), &
                                           elastic_strain_rate_avg => null()

  real(k_real), pointer :: eq_strain => null()
  real(k_real), pointer :: eq_strain_rate => null()
  real(k_real), pointer :: eq_elastic_strain => null()
  real(k_real), pointer :: eq_elastic_strain_rate => null()
  real(k_real), pointer :: eq_plastic_strain => null()
  real(k_real), pointer :: eq_plastic_strain_rate => null()
  real(k_real), pointer :: eq_stress => null()

  real(k_real), pointer :: temperature => null(), time => null()

  logical :: initialized =.false.
  logical :: im_allocated =.false.
  logical :: write_headers =.true.

  type(csv_writer) :: macro_avg_writer
contains
  procedure :: alloc => allocAverageQuantities
  procedure :: init => initAverageQuantities
  procedure :: initGridPointers => initGridPointersAverageQuantities
  procedure :: updateAverageQuantities
  procedure :: getAverageStressPointer
  procedure :: getAverageDisplacementGradientPointer
  procedure :: getAverageStrainRatePointer
  procedure :: getAverageEqStrainRatePointer
  procedure :: getAverageStrainPointer
  procedure :: getAveragePlasticStrainPointer
  procedure :: writeToFile => writeToFileAveragaValues
  procedure :: setTemperaturePointer
end type


type sim_all_macro_data_obj
  type(sim_time_obj), pointer :: sim_time => null()
  ! type(sim_temperature_obj), pointer :: sim_temperature => null()
  type(sim_average_quantities_obj), pointer :: sim_macro_field_averages => null()

contains
  procedure :: init => initAllSimMacrobejct
end type

contains

  !----------------------------------------------------------------------------!
  !------------------------TIME OBJECT SUBROUTINES-----------------------------!
  !----------------------------------------------------------------------------!

  subroutine allocTimeObject(this)
    class(sim_time_obj), intent(inout) :: this

    allocate( this%time, &
              this%time_old, &
              this%dt, &
              this%dt0, &
              this%dt_old, &
              this%dt_reference, &
              this%dtmin, &
              this%dtmax)
    this%im_allocated = .true.
    this%dt=0._k_real
  end subroutine


  subroutine initTimeObject(this, time0, dt, dtmin, dtmax)
    class(sim_time_obj), intent(inout) :: this
    real(k_real), intent(in) :: time0, dt, dtmin, dtmax

    if (.not.(this%im_allocated)) error stop "initTimeObject not allocated. Abort! "
    if (this%initialized) error stop "initTimeObject already intialized"

    if (dt<=0._k_real) error stop "initTimeObject, dt <0. Abort! "
    if (time0<0._k_real) error stop "initTimeObject, time0 <0. Abort! "
    if (dtmin<=0._k_real) error stop "initTimeObject, dtmin <=0. Abort! "
    if (dtmax<=0._k_real) error stop "initTimeObject, dtmax <=0. Abort! "
    if (dtmax<dtmin) error stop "initTimeObject, dtmax<dtmin. Abort! "

    this%time_old = time0
    this%dt0 = dt
    this%dt_old = dt
    this%dt_reference = dt
    this%time = time0
    this%dtmin = dtmin
    this%dtmax = dtmax
    this%initialized = .true.

  end subroutine

  subroutine getTimePointer(this, time_ptr)
    class(sim_time_obj), intent(in) :: this
    real(k_real), intent(out), pointer :: time_ptr

    if (.not.(this%im_allocated)) error stop "TimeObject not intialized"
    if (associated(time_ptr)) error stop "getTimePointer time_ptr already associated. Abort!"
    if (.not.associated(this%time)) error stop "getTimePointer ptr not associated. Abort!"
    time_ptr => this%time
  end subroutine

  subroutine resetDtToDt0(this)
    implicit none
    class(sim_time_obj), intent(inout) :: this
    this%dt = this%dt0
  end subroutine

  subroutine getTimeOldPointer(this, time_old_ptr)
    class(sim_time_obj), intent(in) :: this
    real(k_real), intent(out), pointer :: time_old_ptr

    if (.not.(this%im_allocated)) error stop "TimeObject not intialized"
    if (associated(time_old_ptr)) error stop "getTimeOldPointer time_old_ptr already associated. Abort!"
    if (.not.associated(this%time_old)) error stop "getTimeOldPointer ptr not associated. Abort!"
    time_old_ptr => this%time_old
  end subroutine

  subroutine ForceTime(this, time_to_force)
    class(sim_time_obj), intent(inout) :: this
    real(k_real), intent(in) :: time_to_force

    this%dt = time_to_force - this%time_old
    this%time = time_to_force
    this%forced_time = .true.
  end subroutine

  subroutine getDeltaTimePointer(this, dt_ptr)
    class(sim_time_obj), intent(in) :: this
    real(k_real), intent(out), pointer :: dt_ptr

    if (.not.(this%im_allocated)) error stop "TimeObject not intialized"
    if (associated(dt_ptr)) error stop "getDeltaTimePointer dt_ptr already associated. Abort!"
    if (.not.associated(this%dt)) error stop "getDeltaTimePointer ptr not associated. Abort!"

    dt_ptr => this%dt
  end subroutine

  subroutine getDeltaTimeOldPointer(this, dt_ptr)
    class(sim_time_obj), intent(in) :: this
    real(k_real), intent(out), pointer :: dt_ptr

    if (.not.(this%im_allocated)) error stop "TimeObject not intialized"
    if (associated(dt_ptr)) error stop "getDeltaTimeOldPointer dt_ptr already associated. Abort!"
    if (.not.associated(this%dt_old)) error stop "getDeltaTimeOldPointer ptr not associated. Abort!"

    dt_ptr => this%dt_old
  end subroutine

  subroutine getDeltaTMinPointer(this, dt_ptr)
    class(sim_time_obj), intent(in) :: this
    real(k_real), intent(out), pointer :: dt_ptr

    if (.not.(this%im_allocated)) error stop "TimeObject not intialized"
    if (associated(dt_ptr)) error stop "getDeltaTMinPointer dt_ptr already associated. Abort!"
    if (.not.associated(this%dtmin)) error stop "getDeltaTMinPointer ptr not associated. Abort!"

    dt_ptr => this%dtmin
  end subroutine

  subroutine getDeltaTMaxPointer(this, dt_ptr)
    class(sim_time_obj), intent(in) :: this
    real(k_real), intent(out), pointer :: dt_ptr

    if (.not.(this%im_allocated)) error stop "TimeObject not intialized"
    if (associated(dt_ptr)) error stop "getDeltaTMaxPointer dt_ptr already associated. Abort!"
    if (.not.associated(this%dtmax)) error stop "getDeltaTMaxPointer ptr not associated. Abort!"

    dt_ptr => this%dtmax
  end subroutine

  subroutine getDeltaTimeReferencePointer(this, dt_ptr)
    class(sim_time_obj), intent(in) :: this
    real(k_real), intent(out), pointer :: dt_ptr

    if (.not.(this%im_allocated)) error stop "TimeObject not intialized"
    if (associated(dt_ptr)) error stop "getDeltaTimeReferencePointer dt_ptr already associated. Abort!"
    if (.not.associated(this%dt_old)) error stop "getDeltaTimeReferencePointer ptr not associated. Abort!"

    dt_ptr => this%dt_reference
  end subroutine

  subroutine setDeltaTime(this, new_dt)
    use log_file_mod
    use number_to_string_mod
    implicit none
    class(sim_time_obj), intent(inout) :: this
    real(k_real), intent(in) :: new_dt

    this%dt = min(new_dt, this%dtmax)
    if (this%dt.lt.this%dtmin) then
      log_message = "SIMULATION ABORTED: Request dt (" // trim(adjustl(real2string(this%dt))) //  &
      ") lower than dtmin" // trim(adjustl(real2string(this%dtmin))) // " ."
      call writeToLogFile(log_message)
      error stop "request dt lower thatn dtmin. Abort!"
    endif
  end subroutine

  subroutine setTimeFromReload(this, new_time, new_dt)
    implicit none
    class(sim_time_obj), intent(inout) :: this
    real(k_real), intent(in) :: new_time, new_dt
    this%time = new_time
    this%dt = new_dt
    this%time_old = new_time - new_dt
  end subroutine

  subroutine reduceDeltaTime(this, reduce_factor)
    class(sim_time_obj), intent(inout) :: this
    real(k_real), intent(in) :: reduce_factor
    this%dt = this%dt/reduce_factor
    if (this%dt.lt.this%dtmin) error stop "request dt lower thatn dtmin. Abort!"
    this%time = this%time_old + this%dt
  end subroutine

  subroutine increaseDeltaTime(this, increase_factor)
    class(sim_time_obj), intent(inout) :: this
    real(k_real), intent(in) :: increase_factor
    this%dt = min(this%dt*increase_factor, this%dtmax)
  end subroutine

  subroutine updateCurrentTime(this)
    class(sim_time_obj), intent(inout) :: this
    this%time = this%time_old + this%dt
  end subroutine

  subroutine advanceTime(this)
    class(sim_time_obj), intent(inout) :: this
    if (this%forced_time) then
       this%dt = this%dt_old
       this%forced_time = .false.
    endif

    this%time_old = this%time
    this%time = this%time_old + this%dt

    this%dt_old = this%dt


    if (this%dt<=0._k_real) then
      write(*,*) "advanceTime generated a negative dt"
      write(*,*) "new time ", this%time
      write(*,*) "previous time ", this%time_old
      stop
    endif
  end subroutine

  !----------------------------------------------------------------------------!
  !--------------------TEMPERATURE OBJECT SUBROUTINES--------------------------!
  !----------------------------------------------------------------------------!

  ! subroutine allocTemperatureObject(this)
  !   class(sim_temperature_obj), intent(inout) :: this
  !   allocate(this%temperature)
  !   this%im_allocated =.true.
  ! end subroutine
  !
  ! subroutine initTemperatureObject(this, temperature0)
  !   class(sim_temperature_obj), intent(inout) :: this
  !   real(k_real), intent(in) :: temperature0
  !
  !   if (.not.(this%im_allocated)) error stop "initTemperatureObject not allocated. Abort! "
  !   if (temperature0<0._k_real) error stop "initTemperatureObject, temperature0 <0. Abort! "
  !
  !   this%temperature = temperature0
  !   this%initialized = .true.
  ! end subroutine
  !
  ! subroutine getTemperaturePointer(this, temperature_ptr)
  !   class(sim_temperature_obj), intent(in) :: this
  !   real(k_real), intent(out), pointer :: temperature_ptr
  !
  !   if (.not.(this%im_allocated)) error stop "getTemperaturePointer not intialized. Abort"
  !   if (associated(temperature_ptr)) error stop "getTemperaturePointer temperature_ptr already associated. Abort!"
  !   if (.not.associated(this%temperature)) error stop "getTemperaturePointer ptr not associated. Abort!"
  !   temperature_ptr => this%temperature
  ! end subroutine



  !----------------------------------------------------------------------------!
  !------------------AVG QUANTITIES OBJECT SUBROUTINES-------------------------!
  !----------------------------------------------------------------------------!

  subroutine allocAverageQuantities(this)
    class(sim_average_quantities_obj), intent(inout) :: this

    allocate(this%strain_rate_avg(3,3), &
             this%strain_avg(3,3), &
             this%elastic_strain_avg(3,3), &
             this%elastic_strain_rate_avg(3,3), &
             this%eq_strain, &
             this%eq_strain_rate, &
             this%eq_elastic_strain, &
             this%eq_elastic_strain_rate, &
             this%eq_plastic_strain, &
             this%eq_plastic_strain_rate, &
             this%eq_stress)

    this%im_allocated =.true.
  end subroutine

  subroutine initAverageQuantities(this, grid_data)
    class(sim_average_quantities_obj), intent(inout) :: this
    class(all_grid_data), intent(in), target :: grid_data

    this%grid_data => grid_data
    this%initialized = .true.

    call this%macro_avg_writer%createNewCSV("simulation_macro_averages.csv")
  end subroutine

  subroutine initGridPointersAverageQuantities(this)
    class(sim_average_quantities_obj), intent(inout) :: this
    if (.not.(this%initialized)) error stop  "initGridPointersAverageQuantities not initialized"

    call this%grid_data%getTensor2PointerToVariableByName("velocity_gradient", this%vel_grad)
    call this%grid_data%getAvgTensor2DataPointerByName("velocity_gradient", this%vel_grad_avg )
    call this%grid_data%getTensor2PointerToVariableByName("total_displacement_gradient", this%dis_grad_tot)
    call this%grid_data%getAvgTensor2DataPointerByName("total_displacement_gradient", this%dis_grad_tot_avg )
    call this%grid_data%getTensor2PointerToVariableByName("plastic_strain", this%plastic_strain )
    call this%grid_data%getAvgTensor2DataPointerByName("plastic_strain", this%plastic_strain_avg )
    call this%grid_data%getTensor2PointerToVariableByName("plastic_strain_rate", this%plastic_strain_rate )
    call this%grid_data%getAvgTensor2DataPointerByName("plastic_strain_rate", this%plastic_strain_rate_avg )
    call this%grid_data%getTensor2PointerToVariableByName("cauchy_stress", this%stress )
    call this%grid_data%getAvgTensor2DataPointerByName("cauchy_stress", this%stress_avg )

  end subroutine

  subroutine updateAverageQuantities(this)
    use tensor_math_mod, only : tensor2Norm
    use change_tensor_basis, only : chg_basis_tensor2_to_vector6
    use print_utils_mod, only : printToScreen
    use log_file_mod, only : writeToScreen
    implicit none
    class(sim_average_quantities_obj), intent(inout) :: this
    real(k_real) :: v6(6)
    call this%vel_grad%computeAverage()
    this%strain_rate_avg = 0.5*( this%vel_grad_avg+transpose(this%vel_grad_avg))
    this%eq_strain_rate = tensor2Norm(this%strain_rate_avg )

    call this%dis_grad_tot%computeAverage()
    this%strain_avg = 0.5*( this%dis_grad_tot_avg+transpose(this%dis_grad_tot_avg))
    this%eq_strain = tensor2Norm(this%strain_avg)

    call this%plastic_strain%computeAverage()
    this%eq_plastic_strain = tensor2Norm(this%plastic_strain_avg)

    call this%plastic_strain_rate%computeAverage()
    this%eq_plastic_strain_rate = tensor2Norm(this%plastic_strain_rate_avg)

    this%elastic_strain_avg = this%strain_avg - this%plastic_strain_avg
    this%eq_elastic_strain = tensor2Norm(this%elastic_strain_avg)

    this%elastic_strain_rate_avg = this%strain_rate_avg - this%plastic_strain_rate_avg
    this%eq_elastic_strain_rate = tensor2Norm(this%elastic_strain_rate_avg)

    call this%stress%computeAverage()
    this%eq_stress = tensor2Norm(this%stress_avg)


    call writeToScreen("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    call writeToScreen(" PRINTING AVERAGE VALUES")

    call printToScreen(this%eq_strain_rate,  "avg eq strain rate: ")
    call printToScreen(this%strain_avg, "Total Strain")
    call printToScreen(this%strain_rate_avg, "Total Strain Rate")
    call printToScreen(this%elastic_strain_avg, "Elastic Strain")
    call printToScreen(this%strain_rate_avg-this%plastic_strain_rate_avg, "Elastic Strain Rate")
    call printToScreen(this%plastic_strain_avg, "Plastic Strain")
    call printToScreen(this%plastic_strain_rate_avg, "Plastic Strain Rate")
    call chg_basis_tensor2_to_vector6(this%plastic_strain_rate_avg, v6)
    call printToScreen(v6, "Plastic Strain Rate B-basis")
    call printToScreen(this%stress_avg, "Stress")

    call writeToScreen(" END OF AVERAGE VALUES")
    call writeToScreen("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    ! call this%writeToFile()

  end subroutine

  subroutine setTemperaturePointer(this, temperature_ptr)
    class(sim_average_quantities_obj), intent(inout) :: this
    real(k_real), pointer, intent(in) :: temperature_ptr

    if (.not.associated(temperature_ptr)) error stop  "temperature_ptr not associated "
    if (associated(this%temperature)) error stop  "temperature pointer already associated "
    this%temperature => temperature_ptr

  end subroutine

  subroutine getAverageStressPointer(this, avg_stress_ptr)
    class(sim_average_quantities_obj), intent(in) :: this
    real(k_real), pointer, dimension(:,:), intent(out) :: avg_stress_ptr
    if (.not.(this%initialized)) error stop  "getAverageStressPointer not initialized"
    if (associated(avg_stress_ptr)) error stop  "avg_stress_ptr already associated "

    avg_stress_ptr => this%stress_avg

  end subroutine

  subroutine getAverageDisplacementGradientPointer(this, avg_disp_grad_ptr)
    class(sim_average_quantities_obj), intent(in) :: this
    real(k_real), pointer, dimension(:,:), intent(out) :: avg_disp_grad_ptr
    if (.not.(this%initialized)) error stop  "getAverageDisplacementGradientPointer not initialized"
    if (associated(avg_disp_grad_ptr)) error stop  "avg_disp_grad_ptr already associated "

    avg_disp_grad_ptr => this%dis_grad_tot_avg

  end subroutine

  subroutine getAverageStrainRatePointer(this, avg_strain_rate_ptr)
    class(sim_average_quantities_obj), intent(in) :: this
    real(k_real), pointer, dimension(:,:), intent(out) :: avg_strain_rate_ptr
    if (.not.(this%initialized)) error stop  "getAverageStrainRatePointer not initialized"
    if (associated(avg_strain_rate_ptr)) error stop  "avg_strain_rate_ptr already associated "

    avg_strain_rate_ptr => this%strain_rate_avg

  end subroutine

  subroutine getAverageStrainPointer(this, avg_strain_ptr)
    class(sim_average_quantities_obj), intent(in) :: this
    real(k_real), pointer, dimension(:,:), intent(out) :: avg_strain_ptr
    if (.not.(this%initialized)) error stop  "getAverageStrainPointer not initialized"
    if (associated(avg_strain_ptr)) error stop  "avg_strain_ptr already associated "

    avg_strain_ptr => this%strain_avg

  end subroutine

  subroutine getAveragePlasticStrainPointer(this, avg_pl_strain_ptr)
    class(sim_average_quantities_obj), intent(in) :: this
    real(k_real), pointer, dimension(:,:), intent(out) :: avg_pl_strain_ptr
    if (.not.(this%initialized)) error stop  "getAveragePlasticStrainPointer not initialized"
    if (associated(avg_pl_strain_ptr)) error stop  "avg_pl_strain_ptr already associated "

    avg_pl_strain_ptr => this%plastic_strain_avg

  end subroutine

  subroutine getAverageEqStrainRatePointer(this, avg_eq_strain_rate_ptr)
    class(sim_average_quantities_obj), intent(in) :: this
    real(k_real), pointer, intent(out) :: avg_eq_strain_rate_ptr
    if (.not.(this%initialized)) error stop  "getAverageStrainRatePointer not initialized"
    if (associated(avg_eq_strain_rate_ptr)) error stop  "avg_strain_rate_ptr already associated "

    avg_eq_strain_rate_ptr => this%eq_strain_rate

  end subroutine


subroutine initAllSimMacrObejct(this)
  class(sim_all_macro_data_obj), intent(inout) :: this
  allocate(this%sim_time, this%sim_macro_field_averages)
  call this%sim_time%alloc()
  ! call this%sim_temperature%alloc()
  call this%sim_macro_field_averages%alloc()

  call this%sim_time%getTimePointer(this%sim_macro_field_averages%time)
  ! call this%sim_temperature%getTemperaturePointer(this%sim_macro_field_averages%temperature)
end subroutine

subroutine writeToFileAveragaValues(this)
  use write_to_file_utils_mod
  implicit none
  class(sim_average_quantities_obj), intent(inout) :: this
  integer :: num_write, idx

  num_write = 1
  if (this%write_headers) num_write = num_write +1

  do idx=1,num_write
    call this%macro_avg_writer%openToAppendCSV()
    call this%macro_avg_writer%AppendScalar(this%time, "time", this%write_headers)
    call this%macro_avg_writer%AppendScalar(this%temperature, "temperature", this%write_headers)
    call this%macro_avg_writer%AppendTensor(this%strain_avg, "strain", "total_strain", this%write_headers)
    call this%macro_avg_writer%AppendTensor(this%strain_rate_avg, "strain", "total_strain_rate", this%write_headers)
    call this%macro_avg_writer%AppendTensor(this%elastic_strain_avg, "strain", "elastic_strain", this%write_headers)
    call this%macro_avg_writer%AppendTensor(this%elastic_strain_rate_avg, "strain", "elastic_strain_rate", this%write_headers)
    call this%macro_avg_writer%AppendTensor(this%plastic_strain_avg, "strain", "plastic_strain", this%write_headers)
    call this%macro_avg_writer%AppendTensor(this%plastic_strain_rate_avg, "strain", "plastic_strain_rate", this%write_headers)
    call this%macro_avg_writer%AppendTensor(this%stress_avg, "stress", "cauchy_stress", this%write_headers)
    call this%macro_avg_writer%closeCSVFile()
    this%write_headers = .false.
  enddo

end subroutine

end module
