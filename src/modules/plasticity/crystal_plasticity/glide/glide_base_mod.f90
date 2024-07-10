module glide_base_mod
use kinds
use cp_base_mod, only : cp_base
use all_grid_data_mod, only : all_grid_data
use common_material_parameter_mod, only : common_material_parameter

#include "macro_debug.fpp"

! all glide deformation mode should inheirit from cp_glide_base
type, extends(cp_base) :: cp_glide_base
  real(k_real), pointer, dimension(:,:) :: h_alpha_beta => null() ! the self and latent hardening matrix
contains
  procedure :: initParametersGlideBase
  procedure :: initGridPointers => initGridPointersGlideBase
  procedure :: setPointData=>setPointDataGlideBase
  procedure :: addFieldVariablesToGrid => addFieldVariablesToGridGlideBase
  procedure :: initSchmidClimbTensorInCrystalAxes => initSchmidClimbTensorInCrystalAxesGlideBase
  procedure :: computeGammaDotAndDGammaDotDRssSS=> computeGammaDotAndDGammaDotDRssSSGlideBase
  procedure :: getRotationRateContribution => getRotationRateContributionGlideBase

end type

contains

  subroutine initParametersGlideBase(this, phase_id, common_material_parameter_ptr, use_damage, n_gauss, n_std_dev_gauss_integration, &
                                                elasticity_obj, crystal_paremeters_ptr)
    use stiffness_base_mod, only : stiffness_base   
    use cp_base_mod, only : crystal_paremeters_type                                  
    implicit none
    class(cp_glide_base), intent(inout) :: this
    integer, intent(in) :: phase_id
    class(common_material_parameter), pointer, intent(in) :: common_material_parameter_ptr
    logical, intent(in) :: use_damage
    integer, intent(in) :: n_gauss
    real(k_real), intent(in) :: n_std_dev_gauss_integration
    class(stiffness_base), pointer, intent(in) :: elasticity_obj
    class(crystal_paremeters_type), pointer, intent(in) :: crystal_paremeters_ptr

    this%n_schmid_components = 5
    call this%initParametersCPBase(phase_id, common_material_parameter_ptr, use_damage, n_gauss, n_std_dev_gauss_integration, &
     elasticity_obj, crystal_paremeters_ptr)
    call this%inelastic_strain_rate_name%setString("glide_strain_rate")

    allocate(this%h_alpha_beta(this%n_ss, this%n_ss))
    this%h_alpha_beta = crystal_paremeters_ptr%hardening_matrix

  end subroutine

  subroutine  addFieldVariablesToGridGlideBase(this)
    use grid_data_var_type_mod
    implicit none
    class(cp_glide_base), intent(inout) :: this

    call this%cp_base%addFieldVariablesToGrid()
    associate (all_grid_data_vars => this%grid_data)

    call all_grid_data_vars%addVar("schmid_tensor", ss_vector6)
    call all_grid_data_vars%addVar("gamma_dot", ss_scalar)

  end associate
  end subroutine

  subroutine initGridPointersGlideBase(this)
    class(cp_glide_base), intent(inout) :: this
    call this%cp_base%initGridPointers()

    call this%grid_data%getSSVector6DataPointerByName("schmid_tensor", this%schmid_grid)
    call this%grid_data%getSSScalarDataPointerByName("gamma_dot", this%gammadot_grid)

  end subroutine

  subroutine setPointDataGlideBase(this, ix,iy,iz)
    implicit none
    class(cp_glide_base), intent(inout) :: this
    integer, intent(in) :: ix,iy,iz

    call this%cp_base%setPointData(ix,iy,iz)
    this%schmid_ptr => this%schmid_grid(:,:,ix,iy,iz)
    this%gammadot_ptr => this%gammadot_grid(:,ix,iy,iz)

  end subroutine


  subroutine computeGammaDotAndDGammaDotDRssSSGlideBase(this, rss, gdot, dgdot_drss)
    class(cp_glide_base), intent(inout) :: this
    real(k_real), intent(in) :: rss
    real(k_real), intent(out) :: gdot, dgdot_drss
    __DECL_CLASS_UNUSED_THIS__
    __DECL_UNUSED_REAL__

    __SUPPRESS_CLASS_UNUSED_THIS__
    __SUPPRESS_UNUSED_REAL__(rss)
    __SUPPRESS_UNUSED_REAL_OUT__(gdot)
    __SUPPRESS_UNUSED_REAL_OUT__(dgdot_drss)
    error stop "If you end up here you means you did not override computeGammaDotAndDGammaDotDRssSSGlideBase in your class"

  end subroutine

  subroutine initSchmidClimbTensorInCrystalAxesGlideBase(this)
    use embedded_slip_systems_mod, only : computeSchmidTensorFromNormalAndDirection
    implicit none
    class(cp_glide_base), intent(inout) :: this
    associate(ss_idx=>this%ss_idx, &
              n_ss=> this%n_ss, &
              direction_ca => this%slip_system_direction_ca, &
              normal_ca => this%slip_system_normal_ca, &
              schmid_tensor_ca => this%full_schmid_climb_tensor_ca )

    do ss_idx=1,n_ss
      call computeSchmidTensorFromNormalAndDirection(normal_ca(:,ss_idx), direction_ca(:,ss_idx), schmid_tensor_ca(:,:,ss_idx))
    enddo

    end associate

  end subroutine

  subroutine getRotationRateContributionGlideBase(this, R_dot, ix, iy, iz)
    use tensor_math_mod, only : mat66InnerProdct
    implicit none
    class(cp_glide_base), intent(inout) :: this
    real(k_real), dimension(3,3), intent(out) :: R_dot
    integer, intent(in) :: ix, iy, iz
    real(k_real), dimension(3,3) :: skew_schmid_sample

    call this%setPointData(ix, iy, iz)

    associate(ss_idx => this%ss_idx, &
              R => this%R_crystal2sample_ptr, &
              skew_schmid_tensor_ca => this%skew_schmid_climb_tensor_ca, &
              gdot => this%gammadot_ptr)

    R_dot = 0.
    do ss_idx = 1, this%n_ss
      skew_schmid_sample = matmul(R,matmul(skew_schmid_tensor_ca(:,:,ss_idx), transpose(R)))
      R_dot =  R_dot + gdot(ss_idx) * skew_schmid_sample
    enddo

    end associate

  end subroutine

end module glide_base_mod
