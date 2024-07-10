#include "macro_debug.fpp"

submodule(climb_mod) exponential_climb_smod
implicit none
contains
!
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


  call this%initParametersClimbBase(phase_id, common_material_parameter_ptr,&
                                    use_damage, n_gauss, n_std_dev_gauss_integration, &
                                    elasticity_obj, crystal_paremeters_ptr)
end subroutine

module subroutine readMaterialParametersFromFileExponentialClimb(this, matf_reader)
  use read_from_file_utils, only : file_reader
  implicit none
  class(exponential_climb), intent(inout) :: this
  class(file_reader), intent(inout) :: matf_reader
  real(k_real), dimension(:), pointer :: mode_dep_parameter => null()

  call matf_reader%readVectorParameter( "beta-dot-0[unitless]", this%n_slip_modes, mode_dep_parameter)
  call this%convertModeToSlipSystemParameter(mode_dep_parameter, this%bdot0)
  deallocate(mode_dep_parameter); nullify(mode_dep_parameter)
  call matf_reader%readVectorParameter( "n-exponent[unitless]", this%n_slip_modes, mode_dep_parameter)
  call this%convertModeToSlipSystemParameter(mode_dep_parameter, this%n_exp)
  deallocate(mode_dep_parameter); nullify(mode_dep_parameter)
  call matf_reader%readParameter("n-exponent-damage[unitless]", this%n_exp_damage)
  call matf_reader%readVectorParameter( "initial-critical-resolved-shear-stress[MPa]", this%n_slip_modes, mode_dep_parameter)
  call this%convertModeToSlipSystemParameter(mode_dep_parameter, this%crss0)
  deallocate(mode_dep_parameter); nullify(mode_dep_parameter)
  call matf_reader%readParameter( "rss-std-deviation[MPa]", this%rss_std_dev)
end subroutine

module subroutine  addFieldVariablesToGridExponentialClimb(this)
  use grid_data_var_type_mod
  implicit none
  class(exponential_climb), intent(inout) :: this

  call this%cp_climb_base%addFieldVariablesToGrid()
  associate (all_grid_data_vars => this%grid_data)

  call all_grid_data_vars%addVar("crss_climb", ss_scalar)
  end associate
end subroutine

module subroutine initGridPointersExponentialClimb(this)
  implicit none
  class(exponential_climb), intent(inout) :: this

  call this%cp_climb_base%initGridPointers()

  call this%grid_data%getSSScalarDataPointerByName("crss_climb", &
    this%crss_grid)

end subroutine

module subroutine setPointDataExponentialClimb(this, ix,iy,iz)
  implicit none
  class(exponential_climb), intent(inout) :: this
  integer, intent(in) :: ix,iy,iz

  !calling the parent class to do its job
  call this%cp_climb_base%setPointData(ix,iy,iz)
  this%crss_ptr => this%crss_grid(:,ix,iy,iz)
end subroutine
!
module subroutine computeGammaDotAndDGammaDotDRssSSExponentialClimb(this, rss, gdot, dgdot_drss)
  implicit none
  class(exponential_climb), intent(inout) :: this
  real(k_real), intent(in) :: rss
  real(k_real), intent(out) :: gdot
  real(k_real), intent(out) :: dgdot_drss

  associate (bdot0 => this%bdot0, &
             crss => this%crss_ptr, &
             n => this%n_exp, &
             ss_idx => this%ss_idx )


    gdot = bdot0(ss_idx)*(abs(rss/crss(ss_idx)))**n(ss_idx) * &
                            sign(1._k_real,rss)

    ! we are assuming dcrss/drss = 0
    dgdot_drss = n(ss_idx) * bdot0(ss_idx)*&
                        (abs(rss/crss(ss_idx)))**(n(ss_idx)-1) / crss(ss_idx)

  end associate
end subroutine

module subroutine initStateVariablesAtMaterialPointExponentialClimb(this)
  use test_utils_mod
  implicit none
  class(exponential_climb), intent(inout) :: this

  this%crss_ptr = this%crss0
  if (this%n_gauss > 1) then
  this%std_dev_ptr = this%rss_std_dev
  endif
end subroutine

module subroutine updateStateVariablesAtMaterialPointInnerLoopExponentialClimb(this)
  implicit none
  class(exponential_climb), intent(inout) :: this
  __DECL_CLASS_UNUSED_THIS__

  __SUPPRESS_CLASS_UNUSED_THIS__
  ! nothing to do here
end subroutine

module subroutine updateStateVariablesAtMaterialPointStaggeredExponentialClimb(this)
  implicit none
  class(exponential_climb), intent(inout) :: this
  __DECL_CLASS_UNUSED_THIS__

  __SUPPRESS_CLASS_UNUSED_THIS__
  ! nothing to do here
end subroutine

module subroutine acceptRejectSolutionExponentialClimb(this, dt_max, accept_solution_flag)
  use mpi_useful_routines_mod, only : MPIMaxRelativeIncrementGridSSScalar, MPIMaxIncrementGridSSScalar
  use mpi_variables_mod, only : i_am_mpi_master
  use time_march_mod, only : computeMaxAllowedTimeStepLinear
  implicit none
  class(exponential_climb), intent(inout) :: this
  real(k_real), intent(out) :: dt_max
  logical, intent(out) :: accept_solution_flag

  if (i_am_mpi_master) write(*,*) "*********************************"
  if (i_am_mpi_master) write(*,*) "acceptRejectSolutionExponentialClimb"

  accept_solution_flag = .true.
  dt_max = this%dt_max

  if (i_am_mpi_master) write(*,*) "nothing to check"

  if (accept_solution_flag.and.i_am_mpi_master) write(*,*) "solution might be ACCEPTED"
  if (.not.accept_solution_flag.and.i_am_mpi_master) write(*,*) "solution will be REJECTED"
  if (i_am_mpi_master) write(*,*) "new allowable dt is ", dt_max
  if (i_am_mpi_master) write(*,*) "*********************************"
  if (i_am_mpi_master) write(*,*) ""

end subroutine

end submodule
