module inelastic_strain_mod
use kinds
use material_base_mod, only : material_base
use string_module, only : string_type
use all_grid_data_mod, only : all_grid_data
use simumaltion_macro_data_mod, only : sim_all_macro_data_obj
use gauss_legendre_integration_mod, only : gauss_legendre_integration
use polymorphic_dtype_array_mod, only : dtype_array_ptr
use grid_data_types, only : griddata_vector6
use common_material_parameter_mod, only : common_material_parameter

#include "macro_debug.fpp"

implicit none

type, extends(material_base) :: inelastic_strain_base
  real(k_real), pointer, dimension(:,:,:,:,:) :: stress_grid=> null() ! a pointer to the stress vector
  real(k_real), pointer, dimension(:,:) :: stress_ptr=> null() ! a pointer to the stress vector
  real(k_real), pointer, dimension(:,:,:) :: porosity_grid=> null() ! a pointer to the porosity, only used for damage
  real(k_real), pointer, dimension(:,:,:,:) :: inelastic_strain_rate_grid=> null()

  real(k_real), pointer :: porosity_ptr=> null() ! a pointer to the porosity, only used for damage
  real(k_real),dimension(:), pointer :: inelastic_strain_rate_ptr => null(), &
                                        inelastic_strain_rate_avg => null()
  class(griddata_vector6), pointer :: inelastic_strain_rate_var => null()


  type(string_type) :: inelastic_strain_rate_name

  real(k_real) :: stress6(6) ! a the stress vector in B basis
  logical :: use_damage = .false. ! used to check we initialized the required variables

  integer :: n_gauss
  real(k_real) :: n_std_dev_gauss_integration
  type(gauss_legendre_integration), pointer :: gaussian => null()
  real(k_real), pointer, dimension(:) :: integration_weights=> null(), integration_normalized_x=> null()
  real(k_real), pointer, dimension(:) :: normalized_integration_weight => null(), normalized_pdf_coeff => null(), normalized_pdf_coeff2 => null()
  real(k_real) :: leblond_stress(6), dleblond_stress_dstress(6,6) !-> variables needed if damage is true
  integer :: igauss

  ! damage parameters
  real(k_real) :: n_exp_damage = 0 ! -> the exponent for the damage

  !! paramters
  real(k_real), pointer, dimension(:) :: E_formation_vacancy => null(), &
                                         E_migration_vacancy => null(), &
                                         E_activation_pipe_diffusion => null()
  real(k_real), pointer, dimension(:) :: vacancy_diffusivity_coefficient => null(), &
                                         pipe_diffusivity_coefficient => null()

  integer, pointer :: gb_id_ptr => null()
  integer :: gb_id
  integer, pointer, dimension(:,:,:) :: gb_id_grid => null()

contains

  procedure :: initParametersInelasticStrainBase !-> base function initializing grid_data pointers and allocating the required space
  procedure :: initGridPointers => initGridPointersInelasticStrainBase !-> base function initializing grid_data pointers and allocating the required space

  ! procedure :: computeEpsilonDot => computeEpsilonDotBase!-> compute and store the strain_rate. don't forget to override it
  ! ! procedure that one shall override to get the proper jacobian!
  ! procedure :: computeDEpsilonDotDStress => computeDEpsilonDotDStressBase!->  compute and store the strain_rate. don't forget to override it


  ! this is THE procedure that computes and provide the strain rate at a point in the grid
  procedure :: getStrainRateaAndStressJacobian => getStrainRateaAndStressJacobianBase

  procedure :: computeEpsilonDotAndDepsilonDotDStress => computeEpsilonDotAndDepsilonDotDStressInelasticStrainBase

  procedure :: getRotationRateContribution => getRotationRateContributionInelasticStrainBase

  ! setPointData must always be carryed over trough classes. This is used to set pointers for doing local calculations
  procedure :: setPointData => setPointersInelasticStrainBase
  procedure :: addFieldVariablesToGrid => addFieldVariablesToGridInelasticStrainBase

  procedure :: computeLeblondStress         !-> compute Leblond stress (if needed)

  procedure :: writeInelastiStrainRateToFile => writeInelastiStrainRateToFileInelasticStrain

  procedure :: acceptRejectSolution => acceptRejectSolutionInelasticStrainBase

  procedure :: computeVacancyDiffusivity !-> Compute Diffusivity
  procedure :: computePipeDiffusivity !-> Compute Diffusivity
  procedure :: computeSelfDiffusivity
  procedure :: thermal_equilibrium_concentration
  procedure :: getGaussianValueAndIntegrationWeights

end type inelastic_strain_base

type, extends(dtype_array_ptr) :: inelastic_strain_array
contains
  procedure :: addElement => addInelsticStrainType
  procedure :: addElementFirst => addInelsticStrainTypeFirst
  procedure :: getElementPtr => getInelsticStrainPtr
end type

contains

  subroutine addFieldVariablesToGridInelasticStrainBase(this)
    use grid_data_var_type_mod
    implicit none
    class(inelastic_strain_base), intent(inout) :: this

    if (.not.(this%macro_object_linked)) error stop "addFieldVariablesToGridInelasticStrainBase: you can't add field varaibles to the grid without first linking the global obejcts"
    if ( .not.(this%parameters_initialized)) error stop "addFieldVariablesToGridInelasticStrainBase: you can't add field varaibles to the grid without first initializing a material parameters"

    call this%material_base%addFieldVariablesToGrid()
    associate (all_grid_data_vars => this%grid_data)
    call all_grid_data_vars%addVar(this%inelastic_strain_rate_name%getString(), vector6)
    end associate

    this%grid_variables_provided =.true.
  end subroutine

  subroutine initParametersInelasticStrainBase(this, phase_id, common_material_parameter_ptr, use_damage, n_gauss, n_std_dev_gauss_integration)
    use probability_mod, only : gaussianProbabilityFromMeanAndStdDev
    implicit none
    class(inelastic_strain_base), intent(inout) :: this
    integer, intent(in) :: phase_id
    logical, intent(in) :: use_damage
    integer, intent(in) :: n_gauss
    real(k_real), intent(in) :: n_std_dev_gauss_integration
    class(common_material_parameter), pointer, intent(in) :: common_material_parameter_ptr

    if (this%parameters_initialized) error stop "you can initialize parameters only once"
    if (n_gauss < 1) error stop "initParametersInelasticStrainBase: n_gauss_points_in < 1"

    call this%initParametersMaterialBase(phase_id, common_material_parameter_ptr)
    if (n_gauss>1) then
      if (n_std_dev_gauss_integration <= 0._k_real) error stop "initParametersInelasticStrainBase: n_std_dev_gauss_integration <= 0 but n_gauss > 1"
      this%n_gauss = n_gauss
      this%n_std_dev_gauss_integration = n_std_dev_gauss_integration
      allocate(this%gaussian)
      call this%gaussian%init(n_gauss)
      this%integration_weights => this%gaussian%w_i;
      this%integration_normalized_x => this%gaussian%x_i
    else
      this%n_gauss = 1
    endif

    

    this%use_damage = use_damage
    this%vacancy_diffusivity_coefficient => this%common_material_parameter_ptr%vacancy_diffusivity_coefficient ! point to the bulk value
    this%E_formation_vacancy => this%common_material_parameter_ptr%vacancy_formation_energy ! point to the bulk value
    this%E_migration_vacancy => this%common_material_parameter_ptr%vacancy_migration_energy ! point to the bulk value
    this%E_activation_pipe_diffusion => this%common_material_parameter_ptr%pipe_diffusion_migration_energy ! point to the bulk value
    this%pipe_diffusivity_coefficient => this%common_material_parameter_ptr%pipe_diffusivity_coefficient ! point to the bulk value


  end subroutine

  subroutine getGaussianValueAndIntegrationWeights(this, mean, std_dev, gaussian_value, total_weights)
    use probability_mod, only : gaussianProbabilityFromMeanAndStdDev
    implicit none
    class(inelastic_strain_base), intent(inout) :: this
    real(k_real), intent(in) :: mean, std_dev
    real(k_real), intent(inout) :: gaussian_value(this%n_gauss), total_weights(this%n_gauss)
    real(k_real) :: v_min, v_max, pdf
    integer :: i
    if (this%n_gauss>1) then
    v_min = mean-this%n_std_dev_gauss_integration*std_dev
    v_max = mean+this%n_std_dev_gauss_integration*std_dev
    call this%gaussian%getXiRealSpace(v_min,  v_max, gaussian_value)
    do i = 1, this%n_gauss
      call gaussianProbabilityFromMeanAndStdDev(mean, std_dev, gaussian_value(i), pdf)
      total_weights(i) = std_dev*pdf*this%integration_weights(i)*this%n_std_dev_gauss_integration
    enddo
  else 
    gaussian_value(1) = mean
    total_weights(1) = 1._k_real
  endif
  end subroutine

  subroutine initGridPointersInelasticStrainBase(this)
    implicit none
    class(inelastic_strain_base), intent(inout) :: this

    if (.not.this%macro_object_linked) &
      error stop "initGridPointersInelasticStrainBase you cannot initialized GridPointers without first linking base objects. Abort!"
    if (.not.this%grid_variables_provided) &
      error stop "initGridPointersInelasticStrainBase you cannot initialized GridPointers without first providing the grid variables. Abort!"
    if (this%grid_pointers_linked) error stop "you can link grid pointers only once. Abort!"
    if (.not.(this%grid_data%initialized)) error stop "initGridPointersInelasticStrainBase you cannot init grid pointer without first initializing the grid!"

    call this%material_base%initGridPointers()

    call this%grid_data%getTensor2DataPointerByName("stress", this%stress_grid)
    call this%grid_data%getVector6DataPointerByName(this%inelastic_strain_rate_name%getString(), this%inelastic_strain_rate_grid)
    call this%grid_data%getAvgVector6DataPointerByName(this%inelastic_strain_rate_name%getString(), this%inelastic_strain_rate_avg)
    call this%grid_data%getVector6PointerToVariableByName(this%inelastic_strain_rate_name%getString() , this%inelastic_strain_rate_var)

    if (this%use_damage) call this%grid_data%getScalarDataPointerByName("effective_porosity", this%porosity_grid)

    call this%common_grid_data%getScalarIntegerDataPointerByName("gb_id", this%gb_id_grid)

    this%grid_pointers_linked =.true.
  end subroutine

subroutine setPointersInelasticStrainBase(this, ix, iy, iz)
  use change_tensor_basis, only : chg_basis_tensor2_to_vector6
  implicit none
  class(inelastic_strain_base), intent(inout) :: this
  integer, intent(in) ::  ix, iy, iz
  if (.not.(this%grid_data%initialized)) error stop "setPointersInelasticStrainBase you cannot init grid pointer without first initializing the grid!"

  call this%material_base%setPointData(ix, iy, iz)
  this%stress_ptr => this%stress_grid(:,:, ix, iy, iz)
  call chg_basis_tensor2_to_vector6(this%stress_ptr, this%stress6)
  this%inelastic_strain_rate_ptr => this%inelastic_strain_rate_grid(:, ix, iy, iz)

  if (this%use_damage) then
    this%porosity_ptr => this%porosity_grid(ix, iy, iz)
  else
    if (.not.(associated(this%porosity_ptr))) allocate(this%porosity_ptr)
    this%porosity_ptr = 0._k_real
  endif
  
  this%gb_id_ptr => this%gb_id_grid(ix, iy, iz)
  this%gb_id = min(this%gb_id_ptr, 2)

end subroutine

subroutine getStrainRateaAndStressJacobianBase(this, stress6, epsilon_dot, depsilon_dot_dstress, ix, iy, iz)
  use tensor_math_mod, only : mat66InnerProdct
  implicit none
  class(inelastic_strain_base), intent(inout) :: this
  real(k_real), target, intent(in) ::stress6(6)
  real(k_real), dimension(6), intent(out) :: epsilon_dot
  real(k_real), dimension(6,6), intent(out) :: depsilon_dot_dstress
  integer, intent(in) :: ix, iy, iz
  integer :: i, j

  call this%setPointData(ix, iy, iz)
  this%stress6 = stress6
  if (this%use_damage.and.this%porosity_ptr.ne.0._k_real) then
    call this%computeLeblondStress(stress6)
  else
    this%stress6(6) = 0._k_real
  end if

  if (this%n_gauss == 1) then
    call this%computeEpsilonDotAndDepsilonDotDStress(this%stress6, epsilon_dot, depsilon_dot_dstress)
    ! TODO we shall rewirte all models in terms of flow magnitude nad direction to avoid these nasty loops
  else
    error stop "getStrainRateaAndStressJacobianBase gaussian model not implemented! "
  endif

  if (this%use_damage.and.this%porosity_ptr.ne.0._k_real) then
    depsilon_dot_dstress = matmul(depsilon_dot_dstress, this%dleblond_stress_dstress)
  else
    do j=1,6; do i=1,6
      if (i.eq.6.or.j.eq.6) depsilon_dot_dstress = 0._k_real
    enddo; enddo
  endif

  this%inelastic_strain_rate_ptr(:) = epsilon_dot
end subroutine

subroutine getRotationRateContributionInelasticStrainBase(this, R_dot, ix, iy, iz)
  implicit none
  class(inelastic_strain_base), intent(inout) :: this
  real(k_real), dimension(3,3), intent(out) :: R_dot
  integer, intent(in) :: ix, iy, iz

  call this%setPointData(ix, iy, iz)

  R_dot = 0.
end subroutine

subroutine computeEpsilonDotAndDepsilonDotDStressInelasticStrainBase(this, stress6, epsilon_dot, depsilon_dot_dstress)
  class(inelastic_strain_base), intent(inout) :: this
  real(k_real), target, intent(in) ::stress6(6)
  real(k_real), target, intent(out) :: epsilon_dot(6), depsilon_dot_dstress(6,6)
  __DECL_CLASS_UNUSED_THIS__
  __DECL_UNUSED_VECTOR_PTR__
  __DECL_UNUSED_MATRIX_PTR__

  __SUPPRESS_CLASS_UNUSED_THIS__
  __SUPPRESS_UNUSED_VECTOR_WARNING__(stress6)
  __SUPPRESS_UNUSED_VECTOR_WARNING__(epsilon_dot)
  __SUPPRESS_UNUSED_MATRIX_WARNING__(depsilon_dot_dstress)
  error stop "If you end up here it means you did not override computeEpsilonDotAndDepsilonDotDStress in your class"
end subroutine

subroutine computeLeblondStress(this,stress6)
  implicit none
  class(inelastic_strain_base), intent(inout) :: this
  real(k_real), intent(in) :: stress6(6)
  __DECL_CLASS_UNUSED_THIS__
  __DECL_UNUSED_VECTOR__(6)
  ! this routine shall modify the value of the stress
  __SUPPRESS_CLASS_UNUSED_THIS__
  __SUPPRESS_UNUSED_VECTOR__(stress6)
  error stop "If you end up here it means you did not override computeLeblondStress in your class"
end subroutine

subroutine writeInelastiStrainRateToFileInelasticStrain(this, csv_writer_obj, write_headers)
  use write_to_file_utils_mod
  use change_tensor_basis, only : chg_basis_vector6_to_tensor2
  use csv_writer_mod, only : csv_writer
  use mpi_useful_routines_mod, only : MPIAverageVoxelWeightGridVector
  class(inelastic_strain_base), intent(inout) :: this
  class(csv_writer), intent(inout) :: csv_writer_obj
  logical, intent(in) :: write_headers
  real(k_real) :: strain_rate_tensor(3,3), strain_rate_vector6(6)
  


  call this%isInitialized()

  call MPIAverageVoxelWeightGridVector(this%inelastic_strain_rate_grid, this%phase_fraction_grid, strain_rate_vector6)
  call chg_basis_vector6_to_tensor2(strain_rate_vector6,  strain_rate_tensor)
  call csv_writer_obj%AppendTensor(strain_rate_tensor, "strain", this%inelastic_strain_rate_name%getString(), write_headers)
  call this%inelastic_strain_rate_var%computeAverage()

end subroutine

subroutine acceptRejectSolutionInelasticStrainBase(this, dt_max, accept_solution_flag)
  class(inelastic_strain_base), intent(inout) :: this
  real(k_real), intent(out) :: dt_max
  logical, intent(out) :: accept_solution_flag
  __DECL_CLASS_UNUSED_THIS__

  __SUPPRESS_CLASS_UNUSED_THIS__
  __SUPPRESS_UNUSED_REAL_OUT__(dt_max)
  __SUPPRESS_UNUSED_LOGICAL_OUT__(accept_solution_flag)
  error stop "If you end up here it means you did not override acceptRejectSolutionInelasticStrainBase in your class"
end subroutine

!------------------------------------------------------------------------------!
!                   INELASTIC STRAIN ARRAY SUBROUTINE                          !
!------------------------------------------------------------------------------!
subroutine addInelsticStrainType(this, new_element)
  class(inelastic_strain_array), intent(inout) :: this
  class(inelastic_strain_base), pointer, intent(inout) :: new_element

  call this%extend()
  this%all_pt(this%n_elements)%pt => new_element
  nullify(new_element)
end subroutine

subroutine addInelsticStrainTypeFirst(this, new_element)
  class(inelastic_strain_array), intent(inout) :: this
  class(inelastic_strain_base), pointer, intent(inout) :: new_element

  call this%extendFirstFree()
  this%all_pt(1)%pt => new_element
  nullify(new_element)
end subroutine

subroutine getInelsticStrainPtr(this, idx, ptr)
  class(inelastic_strain_array), intent(in) :: this
  integer, intent(in) :: idx
  class(inelastic_strain_base), pointer, intent(out) :: ptr

  nullify(ptr)
  call this%checkElementExist(idx)

  associate(elem => this%all_pt(idx)%pt)
  select type(elem)
  class is (inelastic_strain_base)
      ptr => elem
  class default
      error stop "getInelsticStrainPtr wrong class"
  end select
  end associate
  if (.not.(associated(ptr))) error stop "getInelsticStrainPtr ptr not associated "

end subroutine

subroutine computeVacancyDiffusivity(this, Dv)
  use math_constants, only : kBSI
  implicit none
  class(inelastic_strain_base), intent(in) :: this
  real(k_real), intent(out) :: Dv

  associate( Diff => this%vacancy_diffusivity_coefficient, &
             T => this%temperature, &
             Ev=> this%E_migration_vacancy )

  Dv = Diff(this%gb_id)*exp(-Ev(this%gb_id)/(kBSI*T))
  end associate
end subroutine

subroutine computeSelfDiffusivity(this, Dv)
  use math_constants, only : kBSI
  implicit none
  class(inelastic_strain_base), intent(in) :: this
  real(k_real), intent(out) :: Dv

  associate( Diff => this%vacancy_diffusivity_coefficient, &
             T => this%temperature, &
             Ev=> this%E_migration_vacancy, &
             Ef=> this%E_formation_vacancy  )

  Dv = Diff(this%gb_id)*exp(-(Ev(this%gb_id)+Ef(this%gb_id))/(kBSI*T))
  end associate
end subroutine

subroutine computePipeDiffusivity(this, Dv)
  use math_constants, only : kBSI
  implicit none
  class(inelastic_strain_base), intent(in) :: this
  real(k_real), intent(out) :: Dv

  associate( Diff => this%pipe_diffusivity_coefficient, &
             T => this%temperature, &
             Ev=> this%E_activation_pipe_diffusion ) !this is the total pipe activation energy

  Dv = Diff(this%gb_id)*exp(-Ev(this%gb_id)/(kBSI*T))
  end associate
end subroutine

subroutine thermal_equilibrium_concentration(this, c_th, dcth_dtauclimb)
  use math_constants, only : kBSI
  ! use print_utils_mod, only : printToScreen
  implicit none
  class(inelastic_strain_base), intent(in) :: this
  real(k_real), intent(out) :: c_th, dcth_dtauclimb

  associate( T => this%temperature, &
             G_f=> this%E_formation_vacancy )

  c_th = exp(-G_f(this%gb_id)/(kBSI*T))
  dcth_dtauclimb = 0._k_real
  end associate

end subroutine

end module inelastic_strain_mod
