module climb_base_mod
use kinds
use cp_base_mod, only : cp_base
use common_material_parameter_mod, only : common_material_parameter
use all_grid_data_mod, only : all_grid_data
use print_utils_mod
#include "macro_debug.fpp"


implicit none

type, extends(cp_base) :: cp_climb_base
integer :: cd_idx = 0 !-> climb_direction_idx
integer :: n_cd = 0  !-> number_climb_direction
real(k_real), pointer :: climb_direction_vector_ca(:,:) => null(), &
                         climb_direction_tensor_ca(:,:,:) => null()
real(k_real), dimension(:,:,:,:,:), pointer :: climb_direction_tensor_grid => null()
real(k_real), dimension(:,:), pointer :: climb_direction_tensor_ptr => null()

real(k_real), pointer, dimension(:) :: cds(:) => null() ! climb direction stress
real(k_real), pointer, dimension(:,:) :: dcds_dstress => null() ! sclimb direction stress dstress
real(k_real), pointer, dimension(:,:,:) :: dclimbdirection_dstress => null() ! sclimb direction stress dstress

contains
  procedure :: initParametersClimbBase
  procedure :: addFieldVariablesToGrid => addFieldVariablesToGridClimbBase
  procedure :: initGridPointers => initGridPointersClimbBase
  procedure :: setPointData => setPointDataClimbBase
  procedure :: initSchmidClimbTensorInCrystalAxes => initSchmidClimbTensorInCrystalAxesClimbBase
  procedure :: computeGammaDotAndDGammaDotDRssSS => computeGammaDotAndDGammaDotDRssSSClimbBase
  procedure :: LeblondComputeResolvedStressAndDirection => LeblondComputeResolvedStressAndDirectionClimb
  procedure :: computeEpsilonDotAndDepsilonDotDStress => computeEpsilonDotAndDepsilonDotDStressClimb
  procedure :: updateSchmidClimbTensorTensorAtMaterialPoint => updateSchmidClimbTensorTensorAtMaterialPointClimbBase
  procedure :: computeClimbDirectionStress 
end type


contains
  subroutine initParametersClimbBase(this, phase_id, common_material_parameter_ptr, &
                                     use_damage, n_gauss, n_std_dev_gauss_integration, &
                                     elasticity_obj, crystal_paremeters_ptr)

    use embedded_slip_systems_mod, only : crystal_type_enum, HCP, computeClimbDirectionTensor
    use HCP_SS_mod, only : getHCP_climb_directions
    use stiffness_base_mod, only : stiffness_base
    use cp_base_mod, only : crystal_paremeters_type                        
    implicit none

    class(cp_climb_base), intent(inout) :: this
    integer, intent(in) :: phase_id
    class(common_material_parameter), pointer, intent(in) :: common_material_parameter_ptr
    logical, intent(in) :: use_damage
    integer, intent(in) :: n_gauss
    real(k_real), intent(in) :: n_std_dev_gauss_integration
    integer :: n_cd_temp
    class(stiffness_base), pointer, intent(in) :: elasticity_obj
    class(crystal_paremeters_type), pointer, intent(in) :: crystal_paremeters_ptr

    this%n_cd = 0
    this%n_schmid_components = 6
    call this%initParametersCPBase(phase_id, common_material_parameter_ptr, &
               use_damage, n_gauss, n_std_dev_gauss_integration, &
       elasticity_obj, crystal_paremeters_ptr)

    call this%inelastic_strain_rate_name%setString("climb_strain_rate")
    select case(crystal_paremeters_ptr%getCrystalType())
    case (HCP)
      associate (cd_idx => this%cd_idx, &
                n_cd => this%n_cd)
        call getHCP_climb_directions(this%climb_direction_vector_ca, crystal_paremeters_ptr%getCAratio(), n_cd_temp )
        
        n_cd = n_cd_temp
        allocate( this%climb_direction_tensor_ca(3,3,n_cd)  )
        allocate(this%cds(n_cd))
        allocate(this%dcds_dstress(6,n_cd))
        allocate(this%dclimbdirection_dstress(6,6, n_cd))
        write(*,*) "n_cd ", n_cd
        
      end associate

    case default
    end select

    if (this%n_cd>0) then
    endif

    call this%initSchmidClimbTensorInCrystalAxes()

  end subroutine

  subroutine computeClimbDirectionStress(this, stress6, with_damage_in)
    use print_utils_mod, only : printToScreen
    use tensor_math_mod, only : doubleContraction
    implicit none
    class(cp_climb_base), intent(inout) :: this
    real(k_real), target, intent(in) ::stress6(6)
    logical, intent(in), optional :: with_damage_in
    logical :: with_damage

    with_damage = .TRUE.
    if (present(with_damage_in)) then
      with_damage = with_damage_in
    else
       with_damage = this%use_damage
    endif

    with_damage = .FALSE.
    associate (cds => this%cds, &
               cdt => this%climb_direction_tensor_ptr, &
               cd_idx => this%cd_idx )

    do cd_idx =1,this%n_cd
      cds(cd_idx) = doubleContraction( stress6, cdt(:,cd_idx) )

    enddo

    ! ALERT compute leblond stress modifies the resolved shear stress,
    ! computes lebolnd stress, and also takes care of all the required derivatives

    if (with_damage.and.(this%porosity_ptr.ne.0._k_real)) then
      call this%computeLeblondStress(stress6)
    else
      this%dcds_dstress(:,:) = this%climb_direction_tensor_ptr(:,:)
      this%dclimbdirection_dstress(:,:,:) = 0._k_real
    end if

    end associate

  end subroutine

  subroutine  addFieldVariablesToGridClimbBase(this)
    use grid_data_var_type_mod
    implicit none
    class(cp_climb_base), intent(inout) :: this

    call this%cp_base%addFieldVariablesToGrid()
    associate (all_grid_data_vars => this%grid_data)

    call all_grid_data_vars%addVar("climb_tensor", ss_vector6)
    call all_grid_data_vars%addVar("beta_dot", ss_scalar)
    if (this%n_cd > 0) &
    call all_grid_data_vars%addVar("climb_direction_tensor", generic_matrix, additional_var_dimensions=(/6,this%n_cd/))
    end associate

  end subroutine

  subroutine initGridPointersClimbBase(this)
    class(cp_climb_base), intent(inout) :: this

    call this%cp_base%initGridPointers()
    call this%grid_data%getSSVector6DataPointerByName("climb_tensor", this%schmid_grid)
    call this%grid_data%getSSScalarDataPointerByName("beta_dot", this%gammadot_grid)
    if (this%n_cd > 0) &
    call this%grid_data%getGenericMatrixDataPointerByName("climb_direction_tensor", this%climb_direction_tensor_grid)
    
  end subroutine

  subroutine setPointDataClimbBase(this, ix,iy,iz)
    implicit none
    class(cp_climb_base), intent(inout) :: this
    integer, intent(in) :: ix,iy,iz

    call this%cp_base%setPointData(ix,iy,iz)
    this%schmid_ptr => this%schmid_grid(:,:,ix,iy,iz)
    this%gammadot_ptr => this%gammadot_grid(:,ix,iy,iz)
    if (this%n_cd > 0) &
    this%climb_direction_tensor_ptr => this%climb_direction_tensor_grid(:,:,ix,iy,iz)
  end subroutine

  subroutine computeEpsilonDotAndDepsilonDotDStressClimb(this, stress6, epsilon_dot, depsilon_dot_dstress)
    use tensor_math_mod, only : vectorNOuterProduct
    class(cp_climb_base), intent(inout) :: this
    

    real(k_real), target, intent(in) ::stress6(6)
    real(k_real), target, intent(out) :: epsilon_dot(6), depsilon_dot_dstress(6,6)
    real(k_real), dimension(this%n_ss) :: rss_without_damage, gdot_without_damage, dgammadot_without_damage_drss
    real(k_real), dimension(6,this%n_ss) :: drss_without_damage_dstress
    integer :: j
    !!! compute strain

    ! WE PRECOMPUTE EVERYTHING WE NEED BEFORE COMPUTING THE STRAIN AND ITS DERIVATIVE

    ! the following are needed to consistently remove teh hydrostaic part of climb
    ! the idea is for epsilon dot 6 we remove the climb component of the undamaged material
    if (this%use_damage.and.(this%porosity_ptr.ne.0._k_real)) then
      call this%computeResolvedStressAndReslovedStressDStress(stress6, with_damage_in=.FALSE.)
      rss_without_damage = this%rss
      drss_without_damage_dstress = this%drss_dstress
    endif
    
    call this%computeResolvedStressAndReslovedStressDStress(stress6)
    ! compute resolved stress might change the location this%drss_dstress is pointing to
    ! associate must be after it
    associate (n_sch => this%n_schmid_components, &
               schmid => this%schmid_ptr, &
               gdot => this%gammadot_ptr, &
               ss_idx => this%ss_idx, &
               drss_dstress => this%drss_dstress, &
               dgdot_drss => this%dgammadot_drss)

    !TODO WE NEED TO FIND A BETTER WAY FOR DOING THE SWITCH OF PROCEDURE
    if (this%use_damage.and.(this%porosity_ptr.ne.0._k_real)) then
     do ss_idx=1,this%n_ss
       call this%computeGammaDotGaussAndDGammaDotDRssCPBase(rss_without_damage(ss_idx), gdot_without_damage(ss_idx), dgammadot_without_damage_drss(ss_idx))
       call this%computeGammaDotGaussAndDGammaDotDRssCPBase(this%rss(ss_idx), gdot(ss_idx), this%dgammadot_drss(ss_idx))
      !  if (this%ix*this%iy*this%iz.eq.1.and.i_am_mpi_master) then
      !   write(*,*) "ss idx ", ss_idx
      !   write(*,*) "gdot_without_damage",  gdot(ss_idx)
      !   write(*,*) "dgammadot_drss", gdot_without_damage(ss_idx)
      !   endif
     enddo
    else 
      do ss_idx=1,this%n_ss
        call this%computeGammaDotGaussAndDGammaDotDRssCPBase(this%rss(ss_idx), gdot(ss_idx), this%dgammadot_drss(ss_idx))
        ! if (this%ix*this%iy*this%iz.eq.1.and.i_am_mpi_master) then
        !   write(*,*) "ss idx ", ss_idx
        !   write(*,*) "gdot_without_damage",  gdot(ss_idx)
        !   endif
      enddo
    endif

    epsilon_dot = 0
    depsilon_dot_dstress = 0

    if (this%use_damage.and.(this%porosity_ptr.ne.0._k_real)) then ! if we use damage
      do ss_idx =1,this%n_ss
        ! compute epsilon dot
        epsilon_dot = epsilon_dot + &
          gdot(ss_idx) * drss_dstress(:,ss_idx) ! drss dstress is now dlambda_dsigma
        epsilon_dot(6) = epsilon_dot(6) - gdot_without_damage(ss_idx) * drss_without_damage_dstress(6, ss_idx) ! remove hydrostatic part


        ! compute d epsilon dot d stress ! in this case we have an additional term
        depsilon_dot_dstress = depsilon_dot_dstress + &
          dgdot_drss(ss_idx) * vectorNOuterProduct(drss_dstress(:,ss_idx), drss_dstress(:,ss_idx), 6 ) + &
          gdot(ss_idx)*this%ddirection_dstress(:,:,ss_idx)

        ! remove hydrostatic part
        do j=1,6
          depsilon_dot_dstress(6,j) = depsilon_dot_dstress(6,j) - dgammadot_without_damage_drss(ss_idx) * drss_without_damage_dstress(j,ss_idx) *drss_without_damage_dstress(6,ss_idx)
        enddo

      enddo
    else
        do ss_idx =1,this%n_ss
          epsilon_dot(1:n_sch) = epsilon_dot(1:n_sch) + gdot(ss_idx) * schmid(1:n_sch,ss_idx)
          depsilon_dot_dstress(1:n_sch,1:n_sch) = depsilon_dot_dstress(1:n_sch,1:n_sch) + &
            vectorNOuterProduct(schmid(1:n_sch,ss_idx), drss_dstress(1:n_sch,ss_idx), n_sch) * dgdot_drss(ss_idx)
        enddo
        epsilon_dot(6) = 0._k_real ! remove hydrostatic part
        depsilon_dot_dstress(6,:) = 0._k_real ! remove hydrostatic part 
    endif

    this%inelastic_strain_rate_ptr(:) = epsilon_dot
    
    end associate

  end subroutine

  subroutine computeGammaDotAndDGammaDotDRssSSClimbBase(this, rss, gdot, dgdot_drss)
    class(cp_climb_base), intent(inout) :: this
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

  subroutine initSchmidClimbTensorInCrystalAxesClimbBase(this)
    use embedded_slip_systems_mod, only : computeClimbTensorFromNormalDirectionAndPsi, &
                                          computeClimbDirectionTensor
    implicit none
    class(cp_climb_base), intent(inout) :: this
    associate(ss_idx=>this%ss_idx, &
              n_ss=> this%n_ss, &
              cd_idx => this%cd_idx, &
              n_cd => this%n_cd, &
              direction_ca => this%slip_system_direction_ca, &
              normal_ca => this%slip_system_normal_ca, &
              climb_tensor_ca => this%full_schmid_climb_tensor_ca, &
              psi_angle => atan(this%ratio_edge_screw), &
              climb_direction_vector_ca => this%climb_direction_vector_ca, &
              climb_direction_tensor_ca => this%climb_direction_tensor_ca )

    do ss_idx=1,n_ss
      call computeClimbTensorFromNormalDirectionAndPsi(normal_ca(:,ss_idx), direction_ca(:,ss_idx), psi_angle, climb_tensor_ca(:,:,ss_idx))
    enddo

    if (n_cd > 0) then
      do cd_idx=1,n_cd
        call computeClimbDirectionTensor(climb_direction_vector_ca(:,cd_idx), climb_direction_tensor_ca(:,:,cd_idx)) 
      enddo
    end if

    end associate

  end subroutine

  subroutine LeblondComputeResolvedStressAndDirectionClimb(this, stress6, schmid6, lambda, n, porosity, rss, flow_dir, dflowdir_dsigma )
    class(cp_climb_base), intent(in) :: this
    real(k_real), intent(in) :: stress6(6), schmid6(6), lambda, n, porosity
    real(k_real), intent(out) :: rss, flow_dir(6), dflowdir_dsigma(6,6)
    real(k_real) :: rss_dev, dLamabda_dsi(6), dLamabda_dsidsj(6,6) ! rss_full, schmid6_dev(6), 

    call this%computeLambdaStressDerivatives(stress6, schmid6, lambda, n, porosity, &
                                              dLamabda_dsi, dLamabda_dsidsj, rss_dev_out = rss_dev )

    rss = lambda + stress6(6) * schmid6(6)
    flow_dir = dLamabda_dsi
    flow_dir(6) = flow_dir(6)+ schmid6(6)*schmid6(6)
    dflowdir_dsigma = dLamabda_dsidsj

    ! rss_full = sum(stress6 * schmid6)
    ! rss =  rss_full * lambda/rss_dev ! this is the value that will be used int the constituive law
    ! schmid6_dev = schmid6; schmid6_dev(6) = 0._k_real

    ! ! the flow direction is the derivative of the rss w.r.t. sigma
    ! flow_dir = schmid6*lambda/rss_dev & !first term
    !            + rss_full/rss_dev * dLamabda_dsi & !second term
    !            - rss_full / rss_dev**2 * lambda * schmid6_dev !third term

    ! do j=1,6
    !   do i=1,6
    ! dflowdir_dsigma(i,j) = (& ! derivatives of the first term
    !                         schmid6(i)*dLamabda_dsi(j)/rss_dev &
    !                         - schmid6(i)*lambda/rss_dev**2*schmid6_dev(j) ) + &
    !                        ( & ! derivatives of the second term
    !                        schmid6(j)/rss_dev * dLamabda_dsi(i) &
    !                        + rss_full/rss_dev * dLamabda_dsidsj(i,j) &
    !                        - rss_full/rss_dev**2*schmid6_dev(j) * dLamabda_dsi(i) ) &
    !                        + ( & ! derivatives of the third term
    !                        2._k_real*rss_full/rss_dev**3*schmid6_dev(j)*lambda*schmid6_dev(i)  &
    !                        - schmid6(j)/rss_dev**2*lambda*schmid6_dev(i) &
    !                        - rss_full/rss_dev**2 * dLamabda_dsi(j) * schmid6_dev(i) )

    ! enddo; enddo
  end subroutine

  subroutine updateSchmidClimbTensorTensorAtMaterialPointClimbBase(this)
    use change_tensor_basis, only : chg_basis_vector5_to_tensor2, chg_basis_vector6_to_tensor2, &
                                    chg_basis_tensor2_to_vector5, chg_basis_tensor2_to_vector6
    implicit none
    class(cp_climb_base), intent(inout) :: this
    real(k_real), dimension(3,3) :: tensor2

    call this%cp_base%updateSchmidClimbTensorTensorAtMaterialPoint() ! do usual stuff

    ! we do this only if we have climb directions 
    if (this%n_cd > 0) then
      associate(cd_idx=>this%cd_idx, &
                R => this%R_crystal2sample_ptr, &
                climb_direction_tensor_ca => this%climb_direction_tensor_ca)
      
      do cd_idx = 1, this%n_cd
        tensor2 = matmul(R,matmul(climb_direction_tensor_ca(:,:,cd_idx), transpose(R)))
        call chg_basis_tensor2_to_vector6(tensor2, this%climb_direction_tensor_ptr(:,cd_idx))
      enddo

      end associate
    endif

  end subroutine

end module climb_base_mod
