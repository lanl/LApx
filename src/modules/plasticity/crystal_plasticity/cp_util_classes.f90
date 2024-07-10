module cp_util_classes_mod
use kinds
use linear_interpolation_mod, only : piecewise_linear_interpolation_array_type
implicit none

type, extends(piecewise_linear_interpolation_array_type) :: slipmode_temperaure_piecewise_linear_parameter
  integer, dimension(:), pointer :: slip_idx_to_mode, i_start_mode, i_end_mode
  integer :: n_slip_modes =-1
  real(k_real) :: old_temperature = -1._k_real
  real(k_real), pointer, dimension(:) :: old_values => null()
contains
  procedure :: readFromFile => readFromFileModeTemperaturePiecewiseLinearParameter
  procedure :: intepolateAllSlipSystems => intepolateAllSSTemperaturePiecewiseLinearParameter
end type

contains
  subroutine readFromFileModeTemperaturePiecewiseLinearParameter(this, matf_reader, n_slip_modes, n_ss_per_mode, parameter_name, value_unit)
    use read_from_file_utils, only : file_reader
    implicit none
    class(slipmode_temperaure_piecewise_linear_parameter), intent(inout) :: this
    class(file_reader), intent(inout) :: matf_reader
    integer, intent(in) ::n_slip_modes, n_ss_per_mode(:)
    character(len=*), intent(in) :: parameter_name, value_unit
    character(len=9), parameter :: mode_base_string="slip-mode"
    integer :: mode_idx, previous_i_end, i_start, i_end

    this%n_slip_modes = n_slip_modes
    allocate(this%slip_idx_to_mode(sum(n_ss_per_mode)), &
             this%i_start_mode(n_slip_modes), &
             this%i_end_mode(n_slip_modes), &
             this%old_values(sum(n_ss_per_mode)))

    i_start = 1
    i_end = n_ss_per_mode(1)
    do mode_idx =1,this%n_slip_modes
      this%slip_idx_to_mode(i_start:i_end) = mode_idx
      this%i_start_mode(mode_idx) = i_start
      this%i_end_mode(mode_idx) = i_end
      if (mode_idx < n_slip_modes) then
        previous_i_end = i_end
        i_start = i_end+1
        i_end = previous_i_end + n_ss_per_mode(mode_idx+1)
      endif
    end do

    call matf_reader%readMultiLinearInterpParameter(n_slip_modes, parameter_name, value_unit, mode_base_string, this)

  end subroutine

  subroutine intepolateAllSSTemperaturePiecewiseLinearParameter(this, temperature, vals)
    use linear_interpolation_mod, only : piecewise_linear_interpolation
    implicit none
    class(slipmode_temperaure_piecewise_linear_parameter), intent(inout) :: this
    real(k_real), intent(in) :: temperature
    real(k_real), dimension(:), intent(out) :: vals
    class(piecewise_linear_interpolation), pointer :: interpolator =>null()
    real(k_real) :: val
    integer :: i

    if (temperature.eq.this%old_temperature) then
      vals = this%old_values
    else
      do i=1,this%n_slip_modes
        nullify(interpolator)
        call this%getElementPtr(i, interpolator)
        call interpolator%interpolate(temperature, val)
        vals(this%i_start_mode(i):this%i_end_mode(i)) = val
      enddo
      this%old_temperature = temperature
      this%old_values = vals
    endif

  end subroutine

end module cp_util_classes_mod
