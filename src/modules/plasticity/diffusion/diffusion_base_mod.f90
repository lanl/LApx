module diffusion_base_mod
use inelastic_strain_mod, only : inelastic_strain_base
use all_grid_data_mod, only : all_grid_data
use common_material_parameter_mod, only : common_material_parameter
use kinds
implicit none

type, extends(inelastic_strain_base) :: diffusion_base
contains
  procedure :: initParametersDiffusionBase
  procedure :: addFieldVariablesToGrid => addFieldVariablesToGridDiffusionBase
  procedure :: initGridPointers => initGridPointersDiffusionBase
  procedure :: setPointData => setPointDataDiffusionBase
  procedure :: computeLeblondStress => computeLeblondStressDiffusionBase
end type
contains

subroutine computeLeblondStressDiffusionBase(this, stress6)
    class(diffusion_base), intent(inout) :: this
    integer :: i, j
    real(k_real), intent(in) :: stress6(6)
    associate (sigma => this%stress6, &
               leblond_stress => this%leblond_stress, &
               dleblond_stress_dstress => this%dleblond_stress_dstress, &
               f => this%porosity_ptr)

    ! this is the modified stress
    do i=1,5
      sigma(i) = 1._k_real/(1._k_real-f)*(&
      (1._k_real+2._k_real/3._k_real*f)*1.5_k_real * stress6(i) )
    enddo

    sigma(6) = f/(1._k_real-f)* 0.75_k_real * stress6(6)

    ! this is the derivative of the modifed stress w.r.t. the original stress
    do j=1,6; do i=1,6
      if (i==j) then
        if (i.le.5) then
          dleblond_stress_dstress(i,i) = 1._k_real/(1._k_real-f)*&
                                        ((1._k_real+2._k_real/3._k_real*f)*1.5_k_real )
        else
          dleblond_stress_dstress(6,6) = f/(1._k_real-f)* &
                                        0.75_k_real
        endif
      else
        dleblond_stress_dstress(i,j) = 0._k_real
      endif
    enddo; enddo

    end associate
end subroutine

subroutine initParametersDiffusionBase(this, phase_id, common_material_parameter_ptr, use_damage)
  implicit none
  class(diffusion_base), intent(inout) :: this
  integer, intent(in) :: phase_id
  class(common_material_parameter), pointer, intent(in) :: common_material_parameter_ptr
  logical, intent(in) :: use_damage
  integer, parameter :: n_guass =1
  real(k_real), parameter :: n_std_dev_gauss_integration = 0._k_real

  call this%inelastic_strain_rate_name%setString("diffusion_strain_rate")
  call this%initParametersInelasticStrainBase(phase_id, common_material_parameter_ptr, use_damage, n_guass, n_std_dev_gauss_integration)
end subroutine

subroutine  addFieldVariablesToGridDiffusionBase(this)
  use grid_data_var_type_mod
  implicit none
  class(diffusion_base), intent(inout) :: this
  call this%inelastic_strain_base%addFieldVariablesToGrid()

  associate (all_grid_data_vars => this%grid_data)

end associate
end subroutine

subroutine initGridPointersDiffusionBase(this)
  implicit none
  class(diffusion_base), intent(inout) :: this

  call this%inelastic_strain_base%initGridPointers()

end subroutine

subroutine setPointDataDiffusionBase(this, ix, iy, iz)
  implicit none
  class(diffusion_base), intent(inout) :: this
  integer, intent(in) ::  ix, iy, iz
  call this%inelastic_strain_base%setPointData(ix, iy, iz)
end subroutine

end module diffusion_base_mod
