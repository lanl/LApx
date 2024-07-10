module cp_base_mod
  use kinds
  use tensor_math_mod, only : doubleContraction, vectorNOuterProduct, mat66InnerProdct
  use inelastic_strain_mod, only : inelastic_strain_base
  use all_grid_data_mod, only : all_grid_data
  use print_utils_mod, only : printToScreen
  use common_material_parameter_mod, only : common_material_parameter
  use test_utils_mod
  use stiffness_base_mod, only : stiffness_base  
  use embedded_slip_systems_mod, only : crystal_type_enum, NONE
  use string_module, only : string_type, string_array
  implicit none

#include "macro_debug.fpp"

type :: crystal_paremeters_type
  type(string_array) :: slip_system_types_str
  real(k_real) :: ratio_edge_screw
  integer, dimension(:), pointer :: n_ss_per_mode => null()
  integer :: n_ss_total, n_slip_modes
  real(k_real), pointer, dimension(:) :: burgVectorL => null(), &
                                         self_hardening => null(), &
                                         latent_hardening => null()
  real(k_real) :: latent_hardening_other_modes
  real(k_real), pointer, dimension(:) :: burgers_vector_length_per_mode, &
                                         burgers_vector_length_per_ss

  real(k_real) :: n_std_dev_gauss_integration, ca_ratio
  real(k_real), pointer, dimension(:,:) :: hardening_matrix  => null()
  integer :: n_gauss
  real(k_real), pointer, dimension(:,:) :: schmid_ca => null(), climb_ca => null()
  real(k_real), pointer, dimension(:,:) :: ss_direction_ca  => null(), ss_normal_ca => null()
  integer, pointer, dimension(:) :: slip_system_2_slip_mode => null()

  ! dislocation loop paramters
  integer :: number_disloc_loop_type
  real(k_real), pointer, dimension(:,:) :: dislocation_loops_ca => null()
  real(k_real), pointer, dimension(:,:) :: loop_crystallography_factor => null()

  integer(kind(crystal_type_enum)) :: crystal_type=NONE
  logical :: use_damage = .false.
  contains

  procedure :: readCrystalParameters
  procedure :: getCrystalType
  generic, public :: setCrystalType => setCrystalTypeString, setCrystalTypeEnum
  procedure, private :: setCrystalTypeString
  procedure, private :: setCrystalTypeEnum
  procedure :: setCAratio
  procedure :: getCAratio
  procedure :: addSSTypeString

end type

type, extends(inelastic_strain_base) :: cp_base
  real(k_real), pointer, dimension(:,:,:,:,:) :: schmid_grid => null() ! pointer to the grid shcmid/climb vector (2nd index is the slip system)
  real(k_real), pointer, dimension(:,:,:,:) :: gammadot_grid => null() ! pointer to flow magnitude (might be gamm_dot or beta_dot)
  real(k_real), pointer, dimension(:,:,:) :: std_dev_grid => null() ! grid pointer to the variance
  real(k_real), pointer, dimension(:,:,:,:,:) :: gaussian_probability_grid => null() ! a pointer to the current gaussian probability
  real(k_real), pointer, dimension(:,:,:,:,:) ::  R_crystal2sample_grid => null()
  real(k_real), pointer, dimension(:,:,:,:,:) ::  stiffness_B_basis_grid => null()
  real(k_real), pointer, dimension(:,:,:,:,:) ::  gauss_legendre_weight_grid => null()
  real(k_real), pointer, dimension(:,:,:,:,:) ::  total_integration_weight_grid => null()

  real(k_real), pointer, dimension(:,:) :: schmid_ptr => null() ! a pointer to the shcmid/ climb vector for all slip systems
  real(k_real), pointer, dimension(:) :: gammadot_ptr => null() ! a pointer to gammadot (might be gammadot or betadot)
  real(k_real), pointer :: std_dev_ptr => null() ! pointer to the local variance
  real(k_real), pointer, dimension(:,:) :: gaussian_probability_ptr => null() ! a pointer to the current gaussian probability
  real(k_real), pointer, dimension(:,:) ::  R_crystal2sample_ptr => null()
  real(k_real), pointer, dimension(:,:) ::  stiffness_B_basis_ptr => null()
  real(k_real), pointer, dimension(:,:) ::  gauss_legendre_weight_ptr => null()
  real(k_real), pointer, dimension(:,:) ::  total_integration_weight_ptr => null()

  real(k_real), pointer, dimension(:,:) :: slip_system_direction_ca=>null(), &
                                           slip_system_normal_ca=>null()
  real(k_real), pointer, dimension(:,:,:) :: full_schmid_climb_tensor_ca, &
                                             sym_schmid_climb_tensor_ca =>null(), &
                                             skew_schmid_climb_tensor_ca =>null()
  real(k_real) :: ratio_edge_screw


  type(stiffness_base), pointer :: elasticity_obj => null() ! elasticity object
  real(k_real), pointer, dimension(:,:,:,:) :: stiffness_crystal_axes_T4_ptr => null()
  real(k_real), pointer :: shear_mod_phase_ptr => null(), &
                           poisson_phase_ptr => null()
  real(k_real), pointer, dimension(:) :: undamged_ss_shear_modulus => null()
  real(k_real), pointer, dimension(:) :: ss_shear_modulus_ptr => null()


  integer :: n_ss = 0 ! the total number of slip systems
  integer, pointer, dimension(:) :: n_ss_per_mode => null(), & ! the number of slip system for each mode
                                    idx_mod_start => null(), &
                                    idx_mod_end => null()
  integer :: n_slip_modes = 0  ! the total number of slip modes
  integer :: n_schmid_components = -1 ! the number of components of the schmid vector
  integer :: ss_idx = -1 ! the slip system index we are working on

 

  ! workspace variables where we do calculation.
  ! later we might decide to link these pointers to griddata variables, if we need the field values
  real(k_real), pointer, dimension(:) :: rss => null() ! storage space for the computed resolved stress
  real(k_real), pointer, dimension(:) :: rss_gaussian => null() ! storage space for the computed resolved stress
  real(k_real), pointer, dimension(:) :: dgammadot_drss => null() ! storage space for the computed dslip_drss
  real(k_real), pointer, dimension(:) :: gammadot_gaussian => null()
  real(k_real), pointer, dimension(:) :: dgammadot_gaussian_drss_gaussian => null()
  real(k_real), pointer, dimension(:,:) :: drss_dstress => null() ! storage space for the computed drss_dstress

  ! damage
  real(k_real), pointer, dimension(:,:,:) :: ddirection_dstress => null()
  real(k_real) :: qH, & !-> hydrostatic stress scaling factor
                  qVM  !-> VM stress scaling factor
contains

  procedure :: initParametersCPBase     !-> initialize pointers and allocate space
  procedure :: addFieldVariablesToGrid => addFieldVariablesToGridCPBase
  procedure :: initGridPointers => initGridPointersCPBase     !-> initialize pointers and allocate space

  procedure :: computeResolvedStressAndReslovedStressDStress!-> a generic procedure to compute the drss_dstress

  procedure :: computeGammaDotAndDGammaDotDRssSS => computeGammaDotAndDGammaDotDRssSSCPBase  !-> compute gamma_dot and dgamma_dot_drss for a single slip system

  procedure :: computeGammaDotGaussAndDGammaDotDRssCPBase

  ! compute computeEpsilonDotCPBase and computeDEpsilonDotDstressCPBase must become
  ! external procedure and must splitted in two, one without damage and one when dmaged is used
  procedure :: computeEpsilonDotAndDepsilonDotDStress => computeEpsilonDotAndDepsilonDotDStressCPBase

  procedure :: setPointData => setPointDataCPBase !-> set all pointers to the current grid point
  procedure :: getStrainRateaAndStressJacobian => getStrainRateaAndStressJacobianCPBase !-> the procedure

  ! LPS procedures
  procedure :: computeLeblondStress => computeLeblondStressCPBase
  procedure :: computeQ
  procedure :: computeh
  procedure :: computeM
  procedure :: computeLeblondStressResidual
  procedure :: NRLeblondStressCPBase
  procedure :: LeblondComputeResolvedStressAndDirection =>  LeblondComputeResolvedStressAndDirectionCPBase
  procedure :: computeLambdaStressDerivatives

  ! support procedure:
  procedure :: convertModeToSlipSystemParameter
  procedure :: getModeBySSIndex

  procedure :: initTexture => initTextureCPBase
  procedure :: updateSchmidClimbTensorTensorAtMaterialPoint => updateSchmidClimbTensorTensorAtMaterialPointCPBase

  ! stiffness related calculation
  procedure :: getShearModulusPoissonRatioAndBulkModulusFromStiffness
  procedure :: updateUndamagedSSShearModulus

  ! test routine
  procedure :: computeDGammaDotDRssSSFD
  procedure :: initSchmidClimbTensorInCrystalAxes => initSchmidClimbTensorInCrystalAxesCPBase

  procedure :: updatePhaseParameters  => updatePhaseParametersCPBase

  procedure :: AveragePerSlipModeSSScalar
  procedure :: writeAvgToCSVByModeAndTotal
end type

interface
  ! damage related routines
  module subroutine computeLeblondStressCPBase(this, stress6)
    implicit none
    class(cp_base), intent(inout) :: this
    real(k_real), intent(in) :: stress6(6)
  end subroutine

  module subroutine computeQ(this, stress6, rss, lambda, schmid, f, q4, Q, dQdsi, dQdLambda, dQ2dsidsj, dQ2dsidLambda, dQ2dLambda2)
    implicit none
    class(cp_base), intent(in) :: this
    real(k_real), intent(in) :: stress6(6), rss, lambda, schmid(5), f, q4
    real(k_real), intent(out) :: Q ! the calue fo Q
    real(k_real), optional, intent(out) :: dQdLambda, & ! used for the jacobian
                                            dQdsi(6), & ! needed to compute the flow direction
                                            dQ2dsidsj(6,6), dQ2dsidLambda(6), dQ2dLambda2 !second derivatives, needed to compute dflow_dstress
  end subroutine

  module subroutine computeM(this, lambda, s6, qM, M, dMdsi, dMdLambda, dM2dLambda2, dM2dsidLambda)
    implicit none
    class(cp_base), intent(in) :: this
    real(k_real), intent(in) :: lambda, s6, qM
    real(k_real), intent(out) :: M
    real(k_real), optional, intent(out) :: dMdLambda, & ! used for the jacobian
                                           dMdsi(6),  & ! needed to compute the flow direction
                                           dM2dLambda2, dM2dsidLambda(6) !second derivatives, needed to compute dflow_dstress
  end subroutine

  module subroutine computeH(this, n, M, h, dhdM, dh2dM2 )
    implicit none
    class(cp_base), intent(in) :: this
    real(k_real), intent(in) :: n, M
    real(k_real), intent(out) :: h
    real(k_real), intent(out), optional :: dhdM, dh2dM2
  end subroutine

  module subroutine computeLeblondStressResidual(this, rss, lambda, schmid, stress6, n, f, R, dR_dLambda, dR_dsigma, dR2_dsidsj, dR2_dsidLambda, dR2_dLambda2)
    implicit none
    class(cp_base), intent(in) :: this
    real(k_real), intent(in)  :: rss, lambda, schmid(5), stress6(6), n, f
    real(k_real), intent(out)  :: R, dR_dLambda
    real(k_real), optional, intent(out)  :: dR_dsigma(6), dR2_dsidsj(6,6), dR2_dsidLambda(6), dR2_dLambda2
  end subroutine

  module subroutine NRLeblondStressCPBase(this, stress6, schmid5, n, porosity, lambda)
    implicit none
    class(cp_base), intent(in) :: this
    real(k_real), intent(in) :: stress6(6), schmid5(5), n, porosity
    real(k_real), intent(out) :: lambda
  end subroutine

  module subroutine LeblondComputeResolvedStressAndDirectionCPBase(this, stress6, schmid6, lambda, n, porosity, rss, flow_dir, dflowdir_dsigma )
    implicit none
    class(cp_base), intent(in) :: this
    real(k_real), intent(in) :: stress6(6), schmid6(6), lambda, n, porosity
    real(k_real), intent(out) :: rss, flow_dir(6), dflowdir_dsigma(6,6)
  end subroutine

  module subroutine computeLambdaStressDerivatives(this, stress6, schmid6, lambda, n, porosity, dLambda_dsi, dLambda2_dsidsj, rss_dev_out )
    implicit none
    class(cp_base), intent(in) :: this
    real(k_real), intent(in) :: stress6(6), schmid6(6), lambda, n, porosity
    real(k_real), intent(out) :: dLambda_dsi(6), dLambda2_dsidsj(6,6)
    real(k_real), optional, intent(out) :: rss_dev_out
    end subroutine
end interface

contains

  subroutine initParametersCPBase(this, phase_id, common_material_parameter_ptr, use_damage, &
                                  n_gauss, n_std_dev_gauss_integration, &
                                  elasticity_obj, crystal_paremeters_ptr)
    use tensor_math_mod, only : getSymmetricPart, getSkewPart
    use stiffness_base_mod, only : stiffness_base
    implicit none
    class(cp_base), intent(inout) :: this
    integer, intent(in) :: phase_id
    class(common_material_parameter), pointer, intent(in) :: common_material_parameter_ptr
    logical, intent(in) :: use_damage
    integer, intent(in) :: n_gauss
    real(k_real), intent(in) :: n_std_dev_gauss_integration
    class(stiffness_base), pointer, intent(in) :: elasticity_obj
    class(crystal_paremeters_type), pointer, intent(in) :: crystal_paremeters_ptr
    integer :: n_slip_modes, idx_mode

    n_slip_modes = sum(shape(crystal_paremeters_ptr%n_ss_per_mode))
    if (n_slip_modes<1) error stop "initCPBase n_slip_modes <1. Abort!"
    this%n_slip_modes = n_slip_modes

    if (any(crystal_paremeters_ptr%n_ss_per_mode<1)) error stop "initCPBase some slip mode less then 1 slip system Abort!"
    allocate(this%n_ss_per_mode(crystal_paremeters_ptr%n_slip_modes), &
             this%idx_mod_start(crystal_paremeters_ptr%n_slip_modes), &
             this%idx_mod_end(crystal_paremeters_ptr%n_slip_modes))
    this%n_ss_per_mode = crystal_paremeters_ptr%n_ss_per_mode
    this%n_ss = sum(crystal_paremeters_ptr%n_ss_per_mode)

    this%idx_mod_start(1) = 1
    this%idx_mod_end(1) = crystal_paremeters_ptr%n_ss_per_mode(1)
    if (this%n_slip_modes>1) then 
      do idx_mode=2,this%n_slip_modes
        this%idx_mod_start(idx_mode) = this%idx_mod_end(idx_mode-1) +1
        this%idx_mod_end(idx_mode) = this%idx_mod_end(idx_mode-1) + crystal_paremeters_ptr%n_ss_per_mode(idx_mode)
      enddo
    endif

    if (.not.(associated(elasticity_obj))) error stop "WKKCTGlide elasticity_obj not associated"
    this%elasticity_obj => elasticity_obj
    this%stiffness_crystal_axes_T4_ptr => this%elasticity_obj%stiffness_crystal_axes
    this%shear_mod_phase_ptr => this%elasticity_obj%G_avg
    this%poisson_phase_ptr => this%elasticity_obj%nu_avg
    allocate(this%undamged_ss_shear_modulus(this%n_ss), &
             this%ss_shear_modulus_ptr(this%n_ss))
    this%undamged_ss_shear_modulus = 0._k_real

    allocate(this%slip_system_direction_ca(3, this%n_ss), &
             this%slip_system_normal_ca(3, this%n_ss), &
             this%sym_schmid_climb_tensor_ca(3,3, this%n_ss), &
             this%skew_schmid_climb_tensor_ca(3,3, this%n_ss), &
             this%full_schmid_climb_tensor_ca(3,3, this%n_ss) )

    this%slip_system_direction_ca = crystal_paremeters_ptr%ss_direction_ca
    this%slip_system_normal_ca = crystal_paremeters_ptr%ss_normal_ca
    this%ratio_edge_screw = crystal_paremeters_ptr%ratio_edge_screw


    allocate (this%rss(this%n_ss), &
              this%dgammadot_drss(this%n_ss))

    call this%initParametersInelasticStrainBase(phase_id, common_material_parameter_ptr, &
                          use_damage, n_gauss, n_std_dev_gauss_integration)

    ! init schmid and climb tensor in crystal axes
    call this%initSchmidClimbTensorInCrystalAxes
    associate(ss_idx=>this%ss_idx, &
              n_ss=> this%n_ss, &
              schmid_climb_tensor_ca => this%full_schmid_climb_tensor_ca, &
              sym_schmid_climb_ca=> this%sym_schmid_climb_tensor_ca, &
              skew_schmid_climb_ca => this%skew_schmid_climb_tensor_ca)
    do ss_idx=1,n_ss
      sym_schmid_climb_ca(:,:,ss_idx) = getSymmetricPart(schmid_climb_tensor_ca(:,:,ss_idx))
      skew_schmid_climb_ca(:,:,ss_idx) = getSkewPart(schmid_climb_tensor_ca(:,:,ss_idx))
    enddo

    end associate

    allocate(this%rss_gaussian(this%n_gauss),&
             this%gammadot_gaussian(this%n_gauss), &
             this%dgammadot_gaussian_drss_gaussian(this%n_gauss) )

    ! allocate damage variables
    allocate( this%drss_dstress(6,this%n_ss), & ! need to allocate drss_dstress
              this%ddirection_dstress(6,6, this%n_ss))
    ! and zero them out
    this%drss_dstress = 0._k_real
    this%ddirection_dstress = 0._k_real

  end subroutine

  subroutine  addFieldVariablesToGridCPBase(this)
    use grid_data_var_type_mod
    implicit none
    class(cp_base), intent(inout) :: this

    call this%inelastic_strain_base%addFieldVariablesToGrid()
    associate (all_grid_data_vars => this%grid_data)
    ! call all_grid_data_vars%addVar("gaussian_probability", ss_generic_vector, &
    !                                  additional_var_dimensions=(/this%n_gauss/))
    ! call all_grid_data_vars%addVar("gauss_legendre_weight", ss_generic_vector, &
    !                                  additional_var_dimensions=(/this%n_gauss/))
    ! call all_grid_data_vars%addVar("total_integration_weight", ss_generic_vector, &
    !                                  additional_var_dimensions=(/this%n_gauss/))
    call all_grid_data_vars%addVar("R_crystal2sample", tensor2)
    call all_grid_data_vars%addVar("stiffness", matrix66)
    call all_grid_data_vars%addVar("standard_deviation", scalar, stateful_level=2)
    end associate

  end subroutine

  subroutine initGridPointersCPBase(this)
    implicit none
    class(cp_base), intent(inout) :: this

    call this%inelastic_strain_base%initGridPointers()

    ! call this%grid_data%getSSGenericVectorDataPointerByName("gaussian_probability", this%gaussian_probability_grid)
    ! call this%grid_data%getSSGenericVectorDataPointerByName("gauss_legendre_weight", this%gauss_legendre_weight_grid)
    ! call this%grid_data%getSSGenericVectorDataPointerByName("total_integration_weight", this%total_integration_weight_grid)
    call this%grid_data%getTensor2DataPointerByName("R_crystal2sample", this%R_crystal2sample_grid)
    call this%grid_data%getMatrix66DataPointerByName("stiffness", this%stiffness_B_basis_grid)
    call this%grid_data%getScalarDataPointerByName("standard_deviation", this%std_dev_grid)
    allocate(this%total_integration_weight_ptr(this%n_gauss,this%n_ss))


  end subroutine

  subroutine setPointDataCPBase(this, ix, iy, iz)
    use change_tensor_basis, only : chg_basis_tensor2_to_vector6
    implicit none
    class(cp_base), intent(inout) :: this
    integer, intent(in) ::  ix, iy, iz

    call this%inelastic_strain_base%setPointData(ix, iy, iz)

    this%schmid_ptr => this%schmid_grid(:,:,ix, iy, iz)
    this%gammadot_ptr => this%gammadot_grid(:,ix,iy,iz)
    call chg_basis_tensor2_to_vector6(this%stress_ptr, this%stress6)
    this%std_dev_ptr => this%std_dev_grid(ix,iy,iz)
    ! this%gaussian_probability_ptr => this%gaussian_probability_grid(:,:,ix,iy,iz)
    ! this%gauss_legendre_weight_ptr => this%gauss_legendre_weight_grid(:,:,ix,iy,iz)
    ! this%total_integration_weight_ptr => this%total_integration_weight_grid(:,:,ix,iy,iz)
    this%R_crystal2sample_ptr => this%R_crystal2sample_grid(:,:,ix,iy,iz)
    this%stiffness_B_basis_ptr => this%stiffness_B_basis_grid(:,:,ix,iy,iz)

    if (this%use_damage) then
      this%ss_shear_modulus_ptr = this%undamged_ss_shear_modulus*(1._k_real - this%porosity_ptr )
    else 
      this%ss_shear_modulus_ptr = this%undamged_ss_shear_modulus
    endif

    if (this%use_damage) this%qH = this%n_ss**(1._k_real/(this%n_exp_damage))
    if (this%use_damage) this%qVM = this%n_ss**(1._k_real/(this%n_exp_damage+1._k_real))
  end subroutine

  subroutine computeResolvedStressAndReslovedStressDStress(this, stress6, with_damage_in)
    use print_utils_mod, only : printToScreen
    implicit none
    class(cp_base), intent(inout) :: this
    real(k_real), target, intent(in) ::stress6(6)
    logical, intent(in), optional :: with_damage_in
    logical :: with_damage

    with_damage = .TRUE.
    if (present(with_damage_in)) then
      with_damage = with_damage_in
    else
       with_damage = this%use_damage
    endif

    associate (n_sch => this%n_schmid_components, &
              rss => this%rss, &
              schmid => this%schmid_ptr, &
              ss_idx => this%ss_idx )

    do ss_idx =1,this%n_ss
      rss(ss_idx) = doubleContraction( stress6(1:n_sch), schmid(1:n_sch,ss_idx) )
    enddo

    ! ALERT compute leblond stress modifies the resolved shear stress,
    ! computes lebolnd stress, and also takes care of all the required derivatives

    if (with_damage.and.(this%porosity_ptr.ne.0._k_real)) then
      call this%computeLeblondStress(stress6)
    else
      this%drss_dstress(1:n_sch,:) = this%schmid_ptr(1:n_sch,:)
      this%ddirection_dstress(:,:,:) = 0._k_real
    end if

    end associate

  end subroutine

  subroutine computeEpsilonDotAndDepsilonDotDStressCPBase(this, stress6, epsilon_dot, depsilon_dot_dstress)
    implicit none
    class(cp_base), intent(inout) :: this
    real(k_real), target, intent(in) ::stress6(6)
    real(k_real), target, intent(out) :: epsilon_dot(6), depsilon_dot_dstress(6,6)
    !!! compute strain

    ! if (this%ix.eq.1.and.this%iy.eq.1.and.this%iz.eq.1.and.this%inelastic_strain_rate_name%getString().eq."dg_climb") then
      ! write(*,*) this%inelastic_strain_rate_name%getString()
    ! ! call printToScreen(this%gammadot_ptr, "gammadot_ptr")
    ! endif

    ! WE PRECOMPUTE EVERYTHING WE NEED BEFORE COMPUTING THE STRAIN AND ITS DERIVATIVE
    call this%computeResolvedStressAndReslovedStressDStress(stress6)

    ! compute resolved stress might change the location this%drss_dstress is pointing to
    ! associate must be after it
    associate (n_sch => this%n_schmid_components, &
               schmid => this%schmid_ptr, &
               gdot => this%gammadot_ptr, &
               ss_idx => this%ss_idx, &
               drss_dstress => this%drss_dstress, &
               dgdot_drss => this%dgammadot_drss)

    do ss_idx=1,this%n_ss
      call this%computeGammaDotGaussAndDGammaDotDRssCPBase(this%rss(ss_idx), this%gammadot_ptr(ss_idx), this%dgammadot_drss(ss_idx))
    enddo

    epsilon_dot = 0
    depsilon_dot_dstress = 0

    if (this%use_damage.and.(this%porosity_ptr.ne.0._k_real)) then ! if we use damage
      do ss_idx =1,this%n_ss
        ! compute epsilon dot
        epsilon_dot = epsilon_dot + &
          gdot(ss_idx) * drss_dstress(:,ss_idx) ! drss dstress is now dlambda_dsigma

        ! compute d epsilon dot d stress ! in this case we have an additional term
        depsilon_dot_dstress = depsilon_dot_dstress + &
          dgdot_drss(ss_idx) * vectorNOuterProduct(drss_dstress(:,ss_idx), drss_dstress(:,ss_idx), 6 ) + &
          gdot(ss_idx)*this%ddirection_dstress(:,:,ss_idx)
      enddo
    else
        do ss_idx =1,this%n_ss
          epsilon_dot(1:n_sch) = epsilon_dot(1:n_sch) + gdot(ss_idx) * schmid(1:n_sch,ss_idx)
          depsilon_dot_dstress(1:n_sch,1:n_sch) = depsilon_dot_dstress(1:n_sch,1:n_sch) + &
            vectorNOuterProduct(schmid(1:n_sch,ss_idx), drss_dstress(1:n_sch,ss_idx), n_sch) * dgdot_drss(ss_idx)
        enddo
    endif

    end associate

  end subroutine

  subroutine computeGammaDotAndDGammaDotDRssSSCPBase(this, rss, gdot, dgdot_drss)
    implicit none
    class(cp_base), intent(inout) :: this
    real(k_real), intent(in) :: rss
    real(k_real), intent(out) :: gdot, dgdot_drss
    __DECL_CLASS_UNUSED_THIS__
    __DECL_UNUSED_REAL__

    __SUPPRESS_CLASS_UNUSED_THIS__
    __SUPPRESS_UNUSED_REAL__(rss)
    __SUPPRESS_UNUSED_REAL_OUT__(gdot)
    __SUPPRESS_UNUSED_REAL_OUT__(dgdot_drss)
    error stop "If you end up here you means you did not override computeGammaDotBase in your class"

  end subroutine

  subroutine computeDGammaDotDRssSSFD(this, rss, gdot, dgdot_drss)
    implicit none
    class(cp_base), intent(inout) :: this
    real(k_real), intent(inout) :: rss, gdot, dgdot_drss
    real(k_real) :: gdot_perturbed, dummy, dgdot_drss_FD
    real(k_real), parameter :: FD_tol =1e-6

    call this%computeGammaDotAndDGammaDotDRssSS( rss+FD_tol, gdot_perturbed, dummy)
    dgdot_drss_FD = (gdot_perturbed-gdot)/FD_tol

    write(*,*) "SS: ", this%ss_idx,  " igauss: ", this%igauss
    write(*,*) "jacobian difference: ", abs(dgdot_drss - dgdot_drss_FD)
    write(*,*) "dgdot_drss analytical: ", dgdot_drss
    write(*,*) "dgdot_drss FD: ", dgdot_drss_FD

    call this%computeGammaDotAndDGammaDotDRssSS( rss, gdot, dgdot_drss)

  end subroutine

  subroutine computeGammaDotGaussAndDGammaDotDRssCPBase(this, rss, gdot, dgdot_drss)
    use probability_mod, only : gaussianProbabilityFromMeanAndStdDev
    use mpi_variables_mod, only : mpi_master_rank
    use math_constants, only : PI
    implicit none
    class(cp_base), intent(inout) :: this
    real(k_real), intent(in) :: rss
    real(k_real), intent(out) :: gdot, dgdot_drss

    associate( gaussian=>this%gaussian, &
                n_gauss => this%n_gauss, &
                rss_gaussian => this%rss_gaussian, &
                gdot_gaussian => this%gammadot_gaussian,&
                dgodtgauss_drssgauss => this%dgammadot_gaussian_drss_gaussian, &
                igauss => this%igauss, &
                std_deviation => this%std_dev_ptr, &
                ! pdf_local => this%gaussian_probability_ptr(:,this%ss_idx), &
                ! gauss_legendre_weights => this%gauss_legendre_weight_ptr(:,this%ss_idx), &
                total_integration_weights => this%total_integration_weight_ptr(:,this%ss_idx), &
                N_std_dev => this%n_std_dev_gauss_integration)

    if (n_gauss>1) then

      call this%getGaussianValueAndIntegrationWeights(rss, std_deviation, rss_gaussian, total_integration_weights)
      ! ! we need to compute the rss for the gaussian
      ! rss_min = rss-N_std_dev*std_deviation
      ! rss_max = rss+N_std_dev*std_deviation

      ! call gaussian%getXiRealSpace(rss_min,  rss_max, rss_gaussian)

      do igauss=1,this%n_gauss
        ! call gaussianProbabilityFromMeanAndStdDev(rss, std_deviation, rss_gaussian(igauss), pdf_local(igauss))
        call this%computeGammaDotAndDGammaDotDRssSS( rss_gaussian(igauss), gdot_gaussian(igauss), dgodtgauss_drssgauss(igauss))
      enddo


      ! total_integration_weights = N_std_dev*std_deviation * pdf_local * this%integration_weights
      gdot = sum(gdot_gaussian*total_integration_weights)
      dgdot_drss = sum(dgodtgauss_drssgauss*total_integration_weights)

    else
      igauss = 1
      call this%computeGammaDotAndDGammaDotDRssSS( rss, gdot_gaussian(1), dgodtgauss_drssgauss(1))
      gdot = gdot_gaussian(1)
      dgdot_drss = dgodtgauss_drssgauss(1)
      ! pdf_local = 1._k_real
      !   ! save integration weights
      ! gauss_legendre_weights = 1._k_real
      !let's save the weighted probability. This will be used to update state variables
      total_integration_weights = 1._k_real
    endif




    end associate
  end subroutine


  subroutine getStrainRateaAndStressJacobianCPBase(this, stress6, epsilon_dot, depsilon_dot_dstress, ix, iy, iz)
    implicit none
    class(cp_base), intent(inout) :: this
    real(k_real), target, intent(in) ::stress6(6)
    real(k_real), dimension(6), intent(out) :: epsilon_dot
    real(k_real), dimension(6,6), intent(out) :: depsilon_dot_dstress
    integer, intent(in) ::  ix, iy, iz

    ! first we set all grid pointers to the local material point
    call this%setPointData(ix, iy, iz)

    call this%computeEpsilonDotAndDepsilonDotDStress(stress6, epsilon_dot, depsilon_dot_dstress)

    this%inelastic_strain_rate_ptr(:) = epsilon_dot
  end subroutine

  subroutine convertModeToSlipSystemParameter(this, param_values_per_mode, param)
    implicit none
    class(cp_base), intent(inout) :: this
    real(k_real), dimension(this%n_slip_modes), intent(in) :: param_values_per_mode
    real(k_real), pointer, dimension(:), intent(inout) :: param
    integer :: i, istart, iend, mode_idx
    if (associated(param)) error stop "convertModeToSlipSystemParameter: already associated parameter"
    allocate(param(this%n_ss))

    mode_idx = 1
    istart = 1
    iend = this%n_ss_per_mode(1)
    do i=1,this%n_slip_modes
      param(istart:iend) = param_values_per_mode(mode_idx)
      mode_idx = mode_idx+1
      if (mode_idx.le.this%n_slip_modes) then
        istart = istart + this%n_ss_per_mode(mode_idx-1)
        iend = iend + this%n_ss_per_mode(mode_idx)
      endif
    enddo

  end subroutine

  function getModeBySSIndex(this, ss_idx) result(i_mode)
    implicit none
    class(cp_base), intent(inout) :: this
    integer, intent(in) :: ss_idx
    integer :: i_mode
    integer :: i, mode_idx_end(this%n_slip_modes)

    mode_idx_end(1) = this%n_ss_per_mode(1)
    if (this%n_slip_modes > 1) then
    do i=2,this%n_slip_modes
      mode_idx_end(i) = mode_idx_end(i-1) + this%n_ss_per_mode(i)
    enddo
    endif

    do i=1,this%n_slip_modes
      if (ss_idx<=mode_idx_end(i)) then
        i_mode = i
        return
      endif
    enddo

    error stop "getModeBySSIndex couldn't find the right mode index"
  end function

  subroutine initTextureCPBase(this)
    implicit none
    class(cp_base), intent(inout) :: this
    integer, dimension(3) :: nx_ny_nz

    nx_ny_nz = shape(this%phase_fraction_grid)

    associate (ix=>this%ix, nx=>nx_ny_nz(1),&
               iy=>this%iy, ny=>nx_ny_nz(2),&
               iz=>this%iz, nz=>nx_ny_nz(3))

      do iz=1,nz
        do iy=1,ny
          do ix=1,nx
            call this%setPointData(ix,iy,iz)
            if (this%phase_fraction_grid(ix,iy,iz).gt.0._k_real) then
            call this%updateSchmidClimbTensorTensorAtMaterialPoint()
            endif
          enddo
        enddo
      enddo

    end associate

  end subroutine

  subroutine initSchmidClimbTensorInCrystalAxesCPBase(this)
    implicit none
    class(cp_base), intent(inout) :: this
    __DECL_CLASS_UNUSED_THIS__
    error stop " if you end up here it means you didn't override comptueSchmidClimbTensorInCrystalAxes in Glide or Climb Base classes"
    __SUPPRESS_CLASS_UNUSED_THIS__
  end subroutine

  subroutine updateSchmidClimbTensorTensorAtMaterialPointCPBase(this)
    use change_tensor_basis, only : chg_basis_vector5_to_tensor2, chg_basis_vector6_to_tensor2, &
                                    chg_basis_tensor2_to_vector5, chg_basis_tensor2_to_vector6
    implicit none
    class(cp_base), intent(inout) :: this
    real(k_real), dimension(3,3) :: tensor2
    associate(ss_idx=>this%ss_idx, &
              R => this%R_crystal2sample_ptr, &
              sym_schmid_climb_ca=> this%sym_schmid_climb_tensor_ca )

    do ss_idx = 1, this%n_ss

      tensor2 = matmul(R,matmul(sym_schmid_climb_ca(:,:,ss_idx), transpose(R)))

      select case(this%n_schmid_components)
      case(5)
        call chg_basis_tensor2_to_vector5(tensor2, this%schmid_ptr(1:5,ss_idx))
      case(6)
        call chg_basis_tensor2_to_vector6(tensor2, this%schmid_ptr(:,ss_idx))
      case  default
        error stop "the number of components not recognized"
      end select
    enddo
    end associate

  end subroutine

  subroutine updateUndamagedSSShearModulus(this)
    use tensor_math_mod, only : T1ij_T1kl_T2ijkl
    implicit none
    class(cp_base), intent(inout) :: this

    associate(ss_idx=>this%ss_idx, &
              G_ss => this%undamged_ss_shear_modulus, &
              schmid_climb_tensor_ca => this%sym_schmid_climb_tensor_ca, &
              sym_schmid_climb_ca=> this%sym_schmid_climb_tensor_ca, &
              C3333 => this%stiffness_crystal_axes_T4_ptr )

      do ss_idx = 1, this%n_ss
        G_ss(ss_idx) = T1ij_T1kl_T2ijkl(sym_schmid_climb_ca(:,:,ss_idx), C3333)
      enddo

    end associate

  endsubroutine

  subroutine getShearModulusPoissonRatioAndBulkModulusFromStiffness(this, shear_modulus, poisson_ratio, bulk_modulus)
    use voigt_indicial_conversion_mod, only : Tensor4ToMatrixVoigt
    use change_tensor_basis, only : chg_basis_matrix66_to_tensor4
    use tensor_math_mod, only : rotateTensor4
    implicit none
    class(cp_base), intent(in) :: this
    real(k_real), intent(out) :: shear_modulus, poisson_ratio, bulk_modulus
    real(k_real), dimension(3,3,3,3) :: stiffnessT4, stiffnessT4_crystal_axes
    real(k_real), dimension(6,6) :: stiffness66_crystal_axes

    call chg_basis_matrix66_to_tensor4(this%stiffness_B_basis_ptr, stiffnessT4)
    call rotateTensor4(stiffnessT4, transpose(this%R_crystal2sample_ptr),  stiffnessT4_crystal_axes)
    stiffness66_crystal_axes = Tensor4ToMatrixVoigt(stiffnessT4_crystal_axes)
    associate(C11 => stiffness66_crystal_axes(1,1), &
              C12 => stiffness66_crystal_axes(1,2), &
              C44 => stiffness66_crystal_axes(4,4))

    shear_modulus = C44
    poisson_ratio = C12/(C11+C12)
    bulk_modulus = (C11+2._k_real*C12)/3._k_real
    end associate
  end subroutine


  subroutine updatePhaseParametersCPBase(this)
    implicit none
    class(cp_base), intent(inout) :: this
    call this%updateUndamagedSSShearModulus()

! #ifdef __DEBUG__
!     write(*,*) "updating phase parameters CPBase"
!     write(*,*) "undamaged ss shear modulus ", this%undamged_ss_shear_modulus
! #endif

  end subroutine

  subroutine AveragePerSlipModeSSScalar(this, my_grid_data, per_mode_average, total)
    use mpi_useful_routines_mod, only : MPIAverageVoxelWeightGridVector
    implicit none
    class(cp_base), intent(in) :: this
    real(k_real), dimension(:,:,:,:), intent(in) :: my_grid_data
    real(k_real), intent(out) :: per_mode_average(this%n_slip_modes)
    real(k_real) :: avg_by_ss(this%n_ss)
    real(k_real), optional :: total
    integer :: idx
    call MPIAverageVoxelWeightGridVector(my_grid_data, this%phase_fraction_grid, avg_by_ss)
    do idx=1,this%n_slip_modes
      per_mode_average(idx) = sum(avg_by_ss(this%idx_mod_start(idx):this%idx_mod_end(idx))) & 
                             /int2real(this%n_ss_per_mode(idx))
    enddo

    if (present(total)) total = sum(avg_by_ss)
  end subroutine


  subroutine writeAvgToCSVByModeAndTotal(this, csv_writer_obj,  var_name, ss_grid_data, write_header, write_total, scaling_factor)
    use mpi_useful_routines_mod, only : MPIAverageVoxelWeightGridVector
    use csv_writer_mod, only : csv_writer  
    
    implicit none
    class(cp_base), intent(in) :: this
    class(csv_writer), intent(inout) :: csv_writer_obj
    character(len=*), intent(in) :: var_name
    logical, intent(in) :: write_header, write_total
    real(k_real), dimension(:,:,:,:), intent(in) :: ss_grid_data
    real(k_real), intent(in), optional :: scaling_factor
    real(k_real) :: avg_by_mode(this%n_slip_modes), total, scaling_factor_
    
    scaling_factor_ = 1._k_real
    if (present(scaling_factor)) scaling_factor_ = scaling_factor

    call this%AveragePerSlipModeSSScalar(ss_grid_data, avg_by_mode, total)
    call csv_writer_obj%AppendVector(this%n_slip_modes, scaling_factor_*avg_by_mode, trim(adjustl(var_name))//"_avg_mode", write_header)
    if (write_total) &
    call csv_writer_obj%AppendScalar(scaling_factor_*total, trim(adjustl(var_name))//"_total", write_header)

  end subroutine
  subroutine readCrystalParameters(this, matf_reader )

    use common_material_parameter_mod
    use read_from_file_utils
    use print_utils_mod, only : printToScreen
    use tensor_math_mod, only : vector3OuterProduct, doubleContraction,doubleContractionBetweenTwoMatrix
    use embedded_slip_systems_mod, only : getNormalizedSlipSystemDirectionandNormal
    use embedded_slip_systems_mod, only : getNormalizedLoopsNormal
    use log_file_mod, only : write_detailed_log_to_screen
    use embedded_slip_systems_mod, only : crystal_type_enum, stringToEnumCrystalType, &
              FCC, BCC, HCP
    use mpi_useful_routines_mod, only : AmIMPIMaster
    implicit none
    class(crystal_paremeters_type) :: this
    class(file_Reader), intent(inout) :: matf_reader
    real(k_real), allocatable, dimension(:,:) ::  ss_direction_ca_temp, ss_normal_ca_temp
    

    ! temporary variables only used by this subroutine
    type(string_type) :: crystal_type
    type(string_array) :: slip_modes,dislocation_loops
    real(k_real), pointer, dimension(:,:) :: temporary_ss_direction_ca=>null(), temporary_ss_normal_ca=>null()
    real(k_real), pointer, dimension(:,:) :: temporary_loop_normal_ca=>null()
    type(string_array) :: dummy_string_array
    type(string_type) :: hardening_matrix_filename
    integer :: i, n_ss , is, ie
    integer :: n_loop, j, ss_idx
    integer(kind(crystal_type_enum)) :: crystal_type_int
    integer :: mode_idx_i, mode_idx_j
    
    real(k_real), dimension(3,3) :: ss_normal_ca_tensor,loop_normal_ca_tensor
    logical :: read_hardening_matrix_from_file
    type(file_Reader) :: hardening_matrix_reader

    read_hardening_matrix_from_file = .False.
    if (associated(this%ss_direction_ca)) error stop "readCrystalParameters ss_direction_ca already associated"
    if (associated(this%ss_normal_ca)) error stop "readCrystalParameters ss_normal_ca already associated"

    call matf_reader%readLineAndCheckStringsAreEqual("--Crystal-Parameters", dummy_string_array)
    call matf_reader%readParameter("crystal-type", crystal_type)
    call stringToEnumCrystalType(crystal_type%getString(),  crystal_type_int)
    call this%setCrystalType(crystal_type_int)

    select case (crystal_type_int)
    case (FCC)
    ! nothing to do
    case (BCC)
    ! nothing to do
    case (HCP)
    call matf_reader%readParameter("c/a-ratio", this%ca_ratio)
    case default
    if (AmIMPIMaster().and.write_detailed_log_to_screen) write(*,*) "the selected crystal type ", crystal_type%getString(), "is not implemented"
    if (AmIMPIMaster().and.write_detailed_log_to_screen) write(*,*) "available options are: FCC!"
    stop
    end select

    call matf_reader%readParameter("ratio-edge-screw", this%ratio_edge_screw)
    call matf_reader%readParameter("n-slip-modes", this%n_slip_modes)

    !sanity check before other improvements
    if(this%n_slip_modes<1) error stop "you need at least 1 slip modes to use crystal plasticity. Abort"
      allocate(this%n_ss_per_mode(this%n_slip_modes))
    call matf_reader%readVectorParameter("slip-modes", this%n_slip_modes, slip_modes)

    do i =1,this%n_slip_modes
      if(.not.(slip_modes%strings(i)%startsWith(crystal_type%getString()))) then
        if (AmIMPIMaster().and.write_detailed_log_to_screen) write(*,*) " you can't select the slip mode ", slip_modes%getStringByIndex(i), " for a ", crystal_type%getString(), " crystal type"
        error stop "Abort!"
      endif
    call this%addSSTypeString(slip_modes%getStringByIndex(i))
    enddo



    associate (n_ss_total => this%n_ss_total)
    n_ss_total = 0
    do i = 1,this%n_slip_modes
      select case (crystal_type%getString())
        case ("HCP")
        call getNormalizedSlipSystemDirectionAndNormal(slip_modes%getStringByIndex(i), n_ss, temporary_ss_direction_ca, temporary_ss_normal_ca, this%ca_ratio)
      case default
        call getNormalizedSlipSystemDirectionAndNormal(slip_modes%getStringByIndex(i), n_ss, temporary_ss_direction_ca, temporary_ss_normal_ca)
      end select
      this%n_ss_per_mode(i) = n_ss
      if (i==1) then
        allocate(this%ss_direction_ca(3,n_ss), this%ss_normal_ca(3,n_ss))
      else
        allocate(ss_direction_ca_temp(3,n_ss_total), ss_normal_ca_temp(3,n_ss_total))
        ss_direction_ca_temp = this%ss_direction_ca
        ss_normal_ca_temp = this%ss_normal_ca
        deallocate(this%ss_direction_ca, this%ss_normal_ca)
        allocate(this%ss_direction_ca(3,n_ss_total+n_ss), this%ss_normal_ca(3,n_ss_total+n_ss))
        this%ss_direction_ca(:,1:n_ss_total) = ss_direction_ca_temp
        this%ss_normal_ca(:,1:n_ss_total) = ss_normal_ca_temp
        deallocate(ss_direction_ca_temp, ss_normal_ca_temp)
      endif
      this%ss_direction_ca(:,n_ss_total+1:n_ss_total+n_ss) = temporary_ss_direction_ca
      this%ss_normal_ca(:,n_ss_total+1:n_ss_total+n_ss) = temporary_ss_normal_ca

      deallocate(temporary_ss_direction_ca, temporary_ss_normal_ca)
      nullify(temporary_ss_direction_ca, temporary_ss_normal_ca)
      n_ss_total = n_ss_total + n_ss
    enddo

    allocate(this%slip_system_2_slip_mode(n_ss_total))
    is = 1
    ie = 0
    do i = 1,this%n_slip_modes
      ie = is + this%n_ss_per_mode(i) - 1
      this%slip_system_2_slip_mode(is:ie) = i
      is = is + this%n_ss_per_mode(i)
    enddo
    
    allocate(this%hardening_matrix(this%n_ss_total, this%n_ss_total))
    call matf_reader%readParameter("read-hardening-matrix-from-file[TRUE/FALSE]", read_hardening_matrix_from_file)
    if (read_hardening_matrix_from_file) then
      call matf_reader%readParameter("hardening-matrix-filename", hardening_matrix_filename)
      call hardening_matrix_reader%openReadTextFile(hardening_matrix_filename%getString())
      call hardening_matrix_reader%readMatrix(this%n_ss_total, this%n_ss_total, this%hardening_matrix)
      ! call printToScreen( this%hardening_matrix, "hardening_matrix")
      call hardening_matrix_reader%closeTextFile()
    else
      call matf_reader%readVectorParameter("self-hardening-coeff", this%n_slip_modes, this%self_hardening)
      call matf_reader%readVectorParameter("latent-hardening-coeff", this%n_slip_modes, this%latent_hardening)
      call matf_reader%readParameter("latent-hardening-other-modes", this%latent_hardening_other_modes)

      ! assemble hardening matrix
      ! compute self and altent hardening matrix
      do j=1,this%n_ss_total
        mode_idx_j = this%slip_system_2_slip_mode(j)
        do i=j,this%n_ss_total
          mode_idx_i = this%slip_system_2_slip_mode(i)
          if (mode_idx_i == mode_idx_j) then
            if (i==j) then
              this%hardening_matrix(i,j) =  this%self_hardening(mode_idx_i)
            else
              this%hardening_matrix(i,j) =  this%latent_hardening(mode_idx_i)
            endif
          else
            this%hardening_matrix(i,j) =  this%latent_hardening_other_modes
          endif
          this%hardening_matrix(j,i) =  this%hardening_matrix(i,j)
        enddo
      enddo
    endif 

    call matf_reader%readVectorParameter("burger-vector-length", this%n_slip_modes, this%burgVectorL)

    call matf_reader%readParameter("n-dislocation-loop-type", this%number_disloc_loop_type)
    write(*,*) this%number_disloc_loop_type, "n_dislocation_loop_type"

    if(this%number_disloc_loop_type.eq.0) then
      call matf_reader%readLineAndCheckStringsAreEqual("dislocation-loops", dummy_string_array)
    else
      if (this%number_disloc_loop_type.lt.0) then
        write(*,*) " the number of dislocation loop types must be >=0, instead I have ", this%number_disloc_loop_type
        error stop "Abort"
      endif
    
      call matf_reader%readVectorParameter("dislocation-loops", this%number_disloc_loop_type, dislocation_loops)
      allocate(this%loop_crystallography_factor(n_ss_total,this%number_disloc_loop_type))
      this%loop_crystallography_factor = 0._k_real
      do i = 1,this%number_disloc_loop_type
        call getNormalizedLoopsNormal(dislocation_loops%getStringByIndex(i), n_loop, temporary_loop_normal_ca)
        do ss_idx = 1, n_ss_total
          ss_normal_ca_tensor = vector3OuterProduct(this%ss_normal_ca(:,ss_idx),this%ss_normal_ca(:,ss_idx))
          do j=1,n_loop
            loop_normal_ca_tensor = vector3OuterProduct(temporary_loop_normal_ca(:,j),temporary_loop_normal_ca(:,j))
            this%loop_crystallography_factor(ss_idx,i) = this%loop_crystallography_factor(ss_idx,i) + (1._k_real/n_loop) &
                                              *doubleContractionBetweenTwoMatrix(ss_normal_ca_tensor,loop_normal_ca_tensor)
          end do
        end do
      deallocate(temporary_loop_normal_ca)
      nullify(temporary_loop_normal_ca)
      end do
    endif

    call matf_reader%skipEmptyLine()

    call matf_reader%readLineAndCheckStringsAreEqual("--Integral-Formulation", dummy_string_array)
    call matf_reader%readParameter("n-gauss", this%n_gauss)
    call matf_reader%readParameter("n-standard-deviation-gaussian-integration", this%n_std_dev_gauss_integration)

    end associate
    if (this%n_gauss<=0) then
    write(*,*) &
              " the number of gauss points must be >=0, instead I have ", this%n_gauss
    stop
    endif
    if (this%n_std_dev_gauss_integration<0) then
    write(*,*) &
             " the gauss point standard deviation must be >0, instead I have ", this%n_std_dev_gauss_integration
    stop
    endif
  end subroutine

  subroutine setCrystalTypeString(this, crystal_type_str)
    use embedded_slip_systems_mod, only : stringToEnumCrystalType
    implicit none
    class(crystal_paremeters_type), intent(inout) :: this
    character(3), intent(in) :: crystal_type_str
  
    call stringToEnumCrystalType(crystal_type_str, this%crystal_type)
  
  end subroutine
  
  subroutine setCrystalTypeEnum(this, crystal_type)
    use embedded_slip_systems_mod, only : stringToEnumCrystalType
    implicit none
    class(crystal_paremeters_type), intent(inout) :: this
    integer(kind(crystal_type_enum)), intent(in) :: crystal_type
  
    this%crystal_type = crystal_type
  
  end subroutine

  function getCrystalType(this) result(crystal_type)
    implicit none
    class(crystal_paremeters_type), intent(in) :: this
    integer(kind(crystal_type_enum)) :: crystal_type
    crystal_type = this%crystal_type
  end function

  subroutine setCAratio(this, ca_ratio)
    implicit none
    class(crystal_paremeters_type), intent(inout) :: this
    real(k_real), intent(in) :: ca_ratio
  
    this%ca_ratio=ca_ratio
  
  end subroutine
  
  function getCAratio(this) result(ca_ratio)
    implicit none
    class(crystal_paremeters_type), intent(inout) :: this
    real(k_real) :: ca_ratio
  
    ca_ratio = this%ca_ratio
  
  end function

  subroutine addSSTypeString(this, ss_type_str)
    implicit none
    class(crystal_paremeters_type), intent(inout) :: this
    character(len=*), intent(in) :: ss_type_str
  
    call this%slip_system_types_str%addString(ss_type_str)
  
  end subroutine

end module
