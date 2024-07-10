module stiffness_base_mod
use kinds
use string_module, only : string_type
use all_grid_data_mod, only : all_grid_data
use simumaltion_macro_data_mod, only : sim_all_macro_data_obj
use polymorphic_dtype_array_mod, only : dtype_array_ptr
use material_base_mod, only : material_base
use matrix_inversion
use print_utils_mod, only : printToScreen
use common_material_parameter_mod, only : common_material_parameter

#include "macro_debug.fpp"

implicit none

type, extends(material_base)  :: stiffness_base
  real(k_real), pointer, dimension(:,:,:,:,:) ::  stiffness_B_basis_grid=> null() ! a pointer to the stiffness
  real(k_real), pointer, dimension(:,:,:,:,:) ::  R_crystal2sample_grid => null()
  real(k_real), pointer, dimension(:,:,:) ::  effective_porosity_grid => null()


  real(k_real), pointer, dimension(:,:) ::  stiffness_B_basis_ptr=> null() ! a pointer to the Stiffness, only used for damage
  real(k_real), pointer, dimension(:,:) ::  R_crystal2sample_ptr => null()
  real(k_real), pointer :: porosity_ptr => null()


  real(k_real), dimension(:,:,:,:), pointer :: stiffness_crystal_axes => null()
  real(k_real), dimension(:,:), pointer  :: stiffness_crystal_axes_voigt => null()
  real(k_real), pointer :: G_avg => null()
  real(k_real), pointer :: nu_avg => null()
contains

  procedure :: initGridPointers => initGridPointersStiffnessBase !-> base function initializing grid_data pointers and allocating the required space
  procedure :: initParametersStiffnessBase

  ! setPointData must always be carryed over trough classes. This is used to set pointers for doing local calculations
  procedure :: setPointData => setPointersStiffnessBase
  procedure :: addFieldVariablesToGrid => addFieldVariablesToGridStiffnessBase

  procedure :: initStateVariablesAtMaterialPoint => initStateVariablesAtMaterialPointStiffnessBase

  procedure :: updateStateVariablesAtMaterialPointStaggered => updateStateVariablesAtMaterialPointStaggeredStiffnessBase

  procedure :: computeStiffnessAtMaterialPoint
  procedure :: computeStiffnessCrystalAxes => computeStiffnessCrystalAxesStiffnessBase
  procedure :: computeCompliance => computeComplianceStiffnessBase
  procedure :: computeAveragePoisson => computeAveragePoissonStiffnessBase
  procedure :: computeAverageShearModulus => computeAverageShearModulusStiffnessBase
  procedure :: updatePhaseParameters => updatePhaseParametersStiffnessBase


end type stiffness_base


type, extends(stiffness_base)  :: stiffness_isotropic
  real(k_real) :: E, E_alpha
  real(k_real) :: nu, nu_alpha
contains
  procedure :: computeStiffnessCrystalAxes => computeStiffnessCrystalAxesIsotropic
  procedure :: readMaterialParametersFromFile => readMaterialParametersFromFileStiffnessIsotropic
  procedure :: initParameters => initParametersStiffnessIsotropic

end type

type, extends(stiffness_base)  :: stiffness_isotropic_shear
  real(k_real) :: G, G_alpha
  real(k_real) :: nu, nu_alpha
contains
  procedure :: computeStiffnessCrystalAxes => computeStiffnessCrystalAxesIsotropicShearPoisson
  procedure :: readMaterialParametersFromFile => readStiffnessIsotropicShearPoisson
  procedure :: initParameters => initParametersStiffnessIsotropicShearPoisson

end type

type, extends(stiffness_base)  :: stiffness_cubic
  real(k_real) :: C11, C11_alpha
  real(k_real) :: C12, C12_alpha
  real(k_real) :: C44, C44_alpha
contains
  procedure :: computeStiffnessCrystalAxes => computeStiffnessCrystalAxesCubic
  procedure :: readMaterialParametersFromFile => readMaterialParametersFromFileStiffnessCubic
  procedure :: initParameters => initParametersStiffnessCubic
  procedure :: computeAveragePoisson => computeAveragePoissonStiffnessCubic
  procedure :: computeAverageShearModulus => computeAverageShearModulusStiffnessCubic


end type

type, extends(stiffness_base)  :: stiffness_hexagonal
  real(k_real) :: C11, C11_alpha
  real(k_real) :: C12, C12_alpha
  real(k_real) :: C13, C13_alpha
  real(k_real) :: C33, C33_alpha
  real(k_real) :: C44, C44_alpha
contains
  procedure :: computeStiffnessCrystalAxes => computeStiffnessCrystalAxesHexagonal
  procedure :: readMaterialParametersFromFile => readMaterialParametersFromFileStiffnessHexagonal
  procedure :: initParameters => initParametersStiffnessHexagonal

end type

contains

  subroutine addFieldVariablesToGridStiffnessBase(this)
    use grid_data_var_type_mod
    implicit none
    class(stiffness_base), intent(inout) :: this

    call this%material_base%addFieldVariablesToGrid()

    associate (all_grid_data_vars => this%grid_data)
      call all_grid_data_vars%addVar("R_crystal2sample", tensor2)
      call all_grid_data_vars%addVar("stiffness", matrix66)
      call all_grid_data_vars%addVar("effective_porosity", scalar)
    end associate

  end subroutine

  subroutine initGridPointersStiffnessBase(this)
    implicit none
    class(stiffness_base), intent(inout) :: this

    call this%material_base%initGridPointers()

    call this%grid_data%getTensor2DataPointerByName("R_crystal2sample", this%R_crystal2sample_grid)
    call this%grid_data%getMatrix66DataPointerByName("stiffness", this%stiffness_B_basis_grid)
    call this%grid_data%getScalarDataPointerByName("effective_porosity", this%effective_porosity_grid)

  end subroutine

  subroutine initParametersStiffnessBase(this, phase_id, common_material_parameter_ptr)
    implicit none
    class(stiffness_base), intent(inout) :: this
    integer, intent(in) :: phase_id
    class(common_material_parameter), pointer, intent(in) :: common_material_parameter_ptr


    if (this%parameters_initialized) error stop "you can initialize parameters only once"
    call this%initParametersMaterialBase(phase_id, common_material_parameter_ptr)

    allocate(this%G_avg)
    this%G_avg = 0._k_real
    allocate(this%nu_avg)
    this%nu_avg = 0._k_real
    allocate(this%stiffness_crystal_axes(3,3,3,3))
    this%stiffness_crystal_axes = 0._k_real
    allocate(this%stiffness_crystal_axes_voigt(6,6))
    this%stiffness_crystal_axes_voigt = 0._k_real
  end subroutine

subroutine setPointersStiffnessBase(this, ix, iy, iz)
  use change_tensor_basis, only : chg_basis_tensor2_to_vector6
  implicit none
  class(stiffness_base), intent(inout) :: this
  integer, intent(in) ::  ix, iy, iz

  call this%material_base%setPointData(ix, iy, iz)

  this%R_crystal2sample_ptr => this%R_crystal2sample_grid(:,:,ix,iy,iz)
  this%stiffness_B_basis_ptr => this%stiffness_B_basis_grid(:,:,ix,iy,iz)
  this%porosity_ptr => this%effective_porosity_grid(ix,iy,iz)

end subroutine

subroutine updatePhaseParametersStiffnessBase(this)
  class(stiffness_base), intent(inout) :: this

  ! the stiffness in cystal axis can be updated only once
  call this%computeStiffnessCrystalAxes
  ! poisson ratio does not change with porosity
  call this%computeAveragePoisson(this%stiffness_crystal_axes_voigt, this%nu_avg)
  ! shear modules G, does cahnge with porosity but only linearly
  call this%computeAverageShearModulus(this%stiffness_crystal_axes_voigt, this%G_avg)
end subroutine

subroutine computeComplianceStiffnessBase(this, c, s)
  class(stiffness_base), intent(in) :: this
  real(k_real), intent(in) :: C(6,6)
  real(k_real), intent(out) :: S(6,6)
  real(k_real), parameter, dimension(6,6) :: S_identity =reshape((/ &
  1, 1, 1, 2, 2, 2,&
  1, 1, 1, 2, 2, 2,&
  1, 1, 1, 2, 2, 2,&
  2, 2, 2, 4, 4, 4,&
  2, 2, 2, 4, 4, 4,&
  2, 2, 2, 4, 4, 4/), (/6,6/))
  __DECL_CLASS_UNUSED_THIS__

  __SUPPRESS_CLASS_UNUSED_THIS__
  call lapackInverseSymmetric(C, S)
  S = S*S_identity
end subroutine

subroutine computeAveragePoissonStiffnessBase(this, c, nu)
  use print_utils_mod, only : printToScreen
  use mpi_variables_mod, only : mpi_rank, mpi_master_rank
  class(stiffness_base), intent(in) :: this
  real(k_real), intent(in) :: C(6,6)
  real(k_real), intent(out) :: nu
  real(k_real) :: S(6,6)

  call this%computeCompliance(C, S)
  nu = 0
  nu = nu - 0.5_k_real*( S(1,2)+S(1,3))/S(1,1)
  nu = nu - 0.5_k_real*( S(2,1)+S(2,3))/S(2,2)
  nu = nu - 0.5_k_real*( S(3,1)+S(3,2))/S(3,3)
  nu = nu/3._k_real
  if(mpi_rank.eq.mpi_master_rank) &
    call printToScreen(nu, "computeAveragePoissonStiffnessBase average poisson ratio")
end subroutine

subroutine computeAveragePoissonStiffnessCubic(this, c, nu)
  use print_utils_mod, only : printToScreen
  use mpi_variables_mod, only : mpi_rank, mpi_master_rank
  class(stiffness_cubic), intent(in) :: this
  real(k_real), intent(in) :: C(6,6)
  real(k_real), intent(out) :: nu
  real(k_real) :: S(6,6)

  call this%computeCompliance(C, S)
  nu = -S(1,2)/S(1,1)
  if(mpi_rank.eq.mpi_master_rank) &
    call printToScreen(nu, "computeAveragePoissonStiffnessCubic average poisson ratio")
end subroutine

subroutine computeAverageShearModulusStiffnessBase(this, C, G)
  use print_utils_mod, only : printToScreen
  use mpi_variables_mod, only : mpi_rank, mpi_master_rank
  class(stiffness_base), intent(in) :: this
  real(k_real), intent(in) :: C(6,6)
  real(k_real), intent(out) :: G
  integer :: i
  __DECL_CLASS_UNUSED_THIS__

  __SUPPRESS_CLASS_UNUSED_THIS__
  G=0
  do i =4,6
    G = G + C(i,i)
  enddo

  G = G/3._k_real
  if (mpi_rank.eq.mpi_master_rank) &
  call printToScreen(G, "computeAverageShearModulusStiffnessBase average shear moduls")

end subroutine

subroutine computeAverageShearModulusStiffnessCubic(this, C, G)
  use print_utils_mod, only : printToScreen
  use mpi_variables_mod, only : mpi_rank, mpi_master_rank
  class(stiffness_cubic), intent(in) :: this
  real(k_real), intent(in) :: C(6,6)
  real(k_real), intent(out) :: G
  __DECL_CLASS_UNUSED_THIS__

  __SUPPRESS_CLASS_UNUSED_THIS__

  G = C(4,4)

  if (mpi_rank.eq.mpi_master_rank) &
  call printToScreen(G, "computeAverageShearModulusStiffnessCubic average shear moduls")

end subroutine

subroutine initStateVariablesAtMaterialPointStiffnessBase(this)
  implicit none
  class(stiffness_base), intent(inout) :: this

  call this%computeStiffnessAtMaterialPoint()
end subroutine

subroutine updateStateVariablesAtMaterialPointStaggeredStiffnessBase(this)
  implicit none
  class(stiffness_base), intent(inout) :: this
  call this%computeStiffnessAtMaterialPoint()
end subroutine

subroutine computeStiffnessAtMaterialPoint(this)
  use tensor_math_mod, only : rotateTensor4
  use change_tensor_basis, only : chg_basis_tensor4_to_matrix66
  implicit none
  class(stiffness_base), intent(inout) :: this
  real(k_real), dimension(3,3,3,3) :: rotated_stiffness
  associate (stiffness_ca => this%stiffness_crystal_axes, &
            R_c2s => this%R_crystal2sample_ptr, &
            stiffness_66 => this%stiffness_B_basis_ptr, &
            P => this%porosity_ptr)

  call rotateTensor4(stiffness_ca, R_c2s, rotated_stiffness)
  call chg_basis_tensor4_to_matrix66(rotated_stiffness, stiffness_66)
  stiffness_66 = stiffness_66 * (1._k_real-P)
  end associate
end subroutine

subroutine computeStiffnessCrystalAxesStiffnessBase(this)
  implicit none
  class(stiffness_base), intent(inout) :: this
  __DECL_CLASS_UNUSED_THIS__

  __SUPPRESS_CLASS_UNUSED_THIS__
  error stop "if you end up here it means you didn't override computeStiffnessCrystalAxesStiffnessBase"
end subroutine

subroutine initParametersStiffnessIsotropic(this, phase_id, common_material_parameter_ptr)
  implicit none
  class(stiffness_isotropic), intent(inout) :: this
  integer, intent(in) :: phase_id
  class(common_material_parameter), pointer, intent(in) :: common_material_parameter_ptr

  call this%initParametersStiffnessBase(phase_id, common_material_parameter_ptr)
end subroutine

subroutine readMaterialParametersFromFileStiffnessIsotropic(this, matf_reader)
  use read_from_file_utils, only : file_reader
  class(stiffness_isotropic), intent(inout) :: this
  class(file_reader), intent(inout) :: matf_reader
  real(k_real), pointer, dimension(:) :: param_values => null()


  call matf_reader%readVectorParameter("E-and-E_alpha", 2, param_values)
  this%E = param_values(1)
  this%E_alpha = param_values(2)
  deallocate(param_values); nullify(param_values)
  call matf_reader%readVectorParameter("nu-and-nu_alpha", 2, param_values)
  this%nu = param_values(1)
  this%nu_alpha = param_values(2)
  deallocate(param_values); nullify(param_values)

end subroutine

subroutine computeStiffnessCrystalAxesIsotropic(this)
  use voigt_indicial_conversion_mod, only : MatrixVoigtToTensor4, Tensor4ToMatrixVoigt
  class(stiffness_isotropic), intent(inout) :: this
  integer :: i,j

  associate(stiffness_mtx_voigt => this%stiffness_crystal_axes_voigt, &
            stiffness_tensor => this%stiffness_crystal_axes, &
            E => this%E + this%E_alpha*this%temperature, &
            nu => this%nu + this%nu_alpha*this%temperature )

  stiffness_mtx_voigt = 0._k_real
  do j=1,6
    do  i=1,6
      if(i.le.3.and.j.le.3) then
        if (i==j) then
          stiffness_mtx_voigt(i,j) = 1._k_real-nu
        else
          stiffness_mtx_voigt(i,j) = nu
        endif
      else if(i.gt.3.and.j.gt.3.and.(i==j)) then
        stiffness_mtx_voigt(i,j) = (1._k_real-2._k_real*nu)/2._k_real
      endif
    enddo
  enddo
  stiffness_mtx_voigt = stiffness_mtx_voigt *E/((1._k_real+nu)*(1._k_real-2._k_real*nu))

  stiffness_tensor = MatrixVoigtToTensor4(stiffness_mtx_voigt)

  ! call printToScreen(Tensor4ToMatrixVoigt(stiffness_tensor), "stiffness isotropic new read")

  end associate
end subroutine

subroutine initParametersStiffnessIsotropicShearPoisson(this, phase_id, common_material_parameter_ptr)
  implicit none
  class(stiffness_isotropic_shear), intent(inout) :: this
  integer, intent(in) :: phase_id
  class(common_material_parameter), pointer, intent(in) :: common_material_parameter_ptr

  call this%initParametersStiffnessBase(phase_id, common_material_parameter_ptr)
end subroutine

subroutine readStiffnessIsotropicShearPoisson(this, matf_reader)
  use read_from_file_utils, only : file_reader
  class(stiffness_isotropic_shear), intent(inout) :: this
  class(file_reader), intent(inout) :: matf_reader
  real(k_real), pointer, dimension(:) :: param_values => null()


  call matf_reader%readVectorParameter("G-and-G_alpha", 2, param_values)
  this%G = param_values(1)
  this%G_alpha = param_values(2)
  deallocate(param_values); nullify(param_values)
  call matf_reader%readVectorParameter("nu-and-nu_alpha", 2, param_values)
  this%nu = param_values(1)
  this%nu_alpha = param_values(2)
  deallocate(param_values); nullify(param_values)

end subroutine

subroutine computeStiffnessCrystalAxesIsotropicShearPoisson(this)
  use voigt_indicial_conversion_mod, only : MatrixVoigtToTensor4, Tensor4ToMatrixVoigt
  class(stiffness_isotropic_shear), intent(inout) :: this
  real(k_real) :: E
  integer :: i,j

  associate(stiffness_mtx_voigt => this%stiffness_crystal_axes_voigt, &
            stiffness_tensor => this%stiffness_crystal_axes, &
            G => this%G + this%G_alpha*this%temperature, &
            nu => this%nu + this%nu_alpha*this%temperature )

  E = 2._k_real*G*(1._k_real+nu)
  stiffness_mtx_voigt = 0._k_real
  do j=1,6
    do  i=1,6
      if(i.le.3.and.j.le.3) then
        if (i==j) then
          stiffness_mtx_voigt(i,j) = 1._k_real-nu
        else
          stiffness_mtx_voigt(i,j) = nu
        endif
      else if(i.gt.3.and.j.gt.3.and.(i==j)) then
        stiffness_mtx_voigt(i,j) = (1._k_real-2._k_real*nu)/2._k_real
      endif
    enddo
  enddo
  stiffness_mtx_voigt = stiffness_mtx_voigt *E/((1._k_real+nu)*(1._k_real-2._k_real*nu))

  stiffness_tensor = MatrixVoigtToTensor4(stiffness_mtx_voigt)

  ! call printToScreen(Tensor4ToMatrixVoigt(stiffness_tensor), "stiffness isotropic new read")

  end associate
end subroutine

subroutine initParametersStiffnessCubic(this, phase_id, common_material_parameter_ptr)
  implicit none
  class(stiffness_cubic), intent(inout) :: this
  integer, intent(in) :: phase_id
  class(common_material_parameter), pointer, intent(in) :: common_material_parameter_ptr

  call this%initParametersStiffnessBase(phase_id, common_material_parameter_ptr)

end subroutine

subroutine readMaterialParametersFromFileStiffnessCubic(this, matf_reader)
  use read_from_file_utils, only : file_reader
  class(stiffness_cubic), intent(inout) :: this
  class(file_reader), intent(inout) :: matf_reader
  real(k_real), pointer, dimension(:) :: param_values => null()

  call matf_reader%readVectorParameter("C1111-and-C1111_alpha", 2, param_values)
  this%C11 = param_values(1)
  this%C11_alpha = param_values(2)
  deallocate(param_values); nullify(param_values)
  call matf_reader%readVectorParameter("C1122-and-C1122_alpha", 2, param_values)
  this%C12 = param_values(1)
  this%C12_alpha = param_values(2)
  deallocate(param_values); nullify(param_values)
  call matf_reader%readVectorParameter("C2323-and-C2323_alpha", 2, param_values)
  this%C44 = param_values(1)
  this%C44_alpha = param_values(2)
  deallocate(param_values); nullify(param_values)

end subroutine

subroutine computeStiffnessCrystalAxesCubic(this)
  use voigt_indicial_conversion_mod, only : MatrixVoigtToTensor4
  class(stiffness_cubic), intent(inout) :: this
  integer :: i,j

  associate(stiffness_mtx_voigt => this%stiffness_crystal_axes_voigt, &
            stiffness_tensor => this%stiffness_crystal_axes, &
            c11 => this%c11 + this%c11_alpha*this%temperature, &
            c12 => this%c12 + this%c12_alpha*this%temperature, &
            c44 => this%c44 + this%c44_alpha*this%temperature )

  stiffness_mtx_voigt = 0._k_real
  do j=1,6
    do  i=1,6
      if(i.le.3.and.j.le.3) then
        if (i==j) then
          stiffness_mtx_voigt(i,j) = c11
        else
          stiffness_mtx_voigt(i,j) = c12
        endif
      else if(i.gt.3.and.j.gt.3.and.(i==j)) then
        stiffness_mtx_voigt(i,j) = c44
      endif
    enddo
  enddo
  stiffness_tensor = MatrixVoigtToTensor4(stiffness_mtx_voigt)
  end associate
end subroutine

subroutine initParametersStiffnessHexagonal(this, phase_id, common_material_parameter_ptr)
  implicit none
  class(stiffness_hexagonal), intent(inout) :: this
  integer, intent(in) :: phase_id
  class(common_material_parameter), pointer, intent(in) :: common_material_parameter_ptr

  call this%initParametersStiffnessBase(phase_id, common_material_parameter_ptr)

end subroutine

subroutine readMaterialParametersFromFileStiffnessHexagonal(this, matf_reader)
  use read_from_file_utils, only : file_reader
  class(stiffness_hexagonal), intent(inout) :: this
  class(file_reader), intent(inout) :: matf_reader
  real(k_real), pointer, dimension(:) :: param_values => null()

  call matf_reader%readVectorParameter( "C1111-and-C1111_alpha", 2, param_values)
  this%C11 = param_values(1)
  this%C11_alpha = param_values(2)
  deallocate(param_values); nullify(param_values)
  call matf_reader%readVectorParameter( "C1122-and-C1122_alpha", 2, param_values)
  this%C12 = param_values(1)
  this%C12_alpha = param_values(2)
  deallocate(param_values); nullify(param_values)
  call matf_reader%readVectorParameter( "C1133-and-C1133_alpha", 2, param_values)
  this%C13 = param_values(1)
  this%C13_alpha = param_values(2)
  deallocate(param_values); nullify(param_values)
  call matf_reader%readVectorParameter( "C3333-and-C3333_alpha", 2, param_values)
  this%C33 = param_values(1)
  this%C33_alpha = param_values(2)
  deallocate(param_values); nullify(param_values)
  call matf_reader%readVectorParameter( "C2323-and-C2323_alpha", 2, param_values)
  this%C44 = param_values(1)
  this%C44_alpha = param_values(2)
  deallocate(param_values); nullify(param_values)

end subroutine

subroutine computeStiffnessCrystalAxesHexagonal(this)
  use voigt_indicial_conversion_mod, only : MatrixVoigtToTensor4
  class(stiffness_hexagonal), intent(inout) :: this

  associate(stiffness_mtx_voigt => this%stiffness_crystal_axes_voigt, &
            stiffness_tensor => this%stiffness_crystal_axes, &
            c11 => this%c11 + this%c11_alpha*this%temperature, &
            c12 => this%c12 + this%c12_alpha*this%temperature, &
            c13 => this%c13 + this%c13_alpha*this%temperature, &
            c33 => this%c33 + this%c33_alpha*this%temperature, &
            c44 => this%c44 + this%c44_alpha*this%temperature )

  stiffness_mtx_voigt = 0._k_real
  stiffness_mtx_voigt(1,1) = c11
  stiffness_mtx_voigt(2,2) = c11
  stiffness_mtx_voigt(3,3) = c33

  stiffness_mtx_voigt(1,2) = c12
  stiffness_mtx_voigt(2,1) = c12

  stiffness_mtx_voigt(1,3) = c13
  stiffness_mtx_voigt(2,3) = c13
  stiffness_mtx_voigt(3,1) = c13
  stiffness_mtx_voigt(3,2) = c13

  stiffness_mtx_voigt(4,4) = c44
  stiffness_mtx_voigt(5,5) = c44
  stiffness_mtx_voigt(6,6) = (c11-c12)/2._k_real


  stiffness_tensor = MatrixVoigtToTensor4(stiffness_mtx_voigt)
  end associate
end subroutine

subroutine readElasticityFromFile(matf_reader, phase_id, all_mighty_grid_in, sim_all_macro_data, the_bc_object, elasticity_base_ptr, common_material_parameter_ptr)
  use read_from_file_utils, only : file_reader
  use all_mighty_grid_mod, only : all_mighty_grid_type
  use simumaltion_macro_data_mod, only : sim_all_macro_data_obj
  use bc_objects_mod, only : boundary_condition_array_type
  use common_material_parameter_mod, only : common_material_parameter
  ! use macro
  implicit none
  type(file_reader), intent(inout) :: matf_reader
  class(stiffness_base), pointer, intent(inout) :: elasticity_base_ptr
  class(common_material_parameter), pointer, intent(inout) :: common_material_parameter_ptr
  type(string_type) :: elasticity_model_name
  integer, intent(in) :: phase_id
  class(all_mighty_grid_type), intent(in), target :: all_mighty_grid_in
  class(sim_all_macro_data_obj), intent(in), target :: sim_all_macro_data
  class(boundary_condition_array_type), intent(in), target :: the_bc_object

  class(stiffness_isotropic), pointer :: stiffness_model_isotropic_temp => null()
  class(stiffness_isotropic_shear), pointer :: stiffness_model_isotropic_shear_temp => null()
  class(stiffness_cubic), pointer :: stiffness_model_cubic_temp => null()
  class(stiffness_hexagonal), pointer :: stiffness_model_hexagonal_temp => null()

  call matf_reader%readParameter("--Elasticity", elasticity_model_name)
  select case(elasticity_model_name%getString())
  

  case ("isotropic-linear")
    allocate(stiffness_model_isotropic_temp)
    elasticity_base_ptr => stiffness_model_isotropic_temp
    call stiffness_model_isotropic_temp%linkBaseObjects(phase_id, all_mighty_grid_in, sim_all_macro_data, the_bc_object)
    call stiffness_model_isotropic_temp%initParameters(phase_id, common_material_parameter_ptr)

  case ("isotropic-shear-poisson-linear")
    allocate(stiffness_model_isotropic_shear_temp)
    elasticity_base_ptr => stiffness_model_isotropic_shear_temp
    call stiffness_model_isotropic_shear_temp%linkBaseObjects(phase_id, all_mighty_grid_in, sim_all_macro_data, the_bc_object)
    call stiffness_model_isotropic_shear_temp%initParameters(phase_id, common_material_parameter_ptr)

  case ("cubic-linear")
    allocate(stiffness_model_cubic_temp)
    elasticity_base_ptr => stiffness_model_cubic_temp
    call stiffness_model_cubic_temp%linkBaseObjects(phase_id, all_mighty_grid_in, sim_all_macro_data, the_bc_object)
    call stiffness_model_cubic_temp%initParameters(phase_id, common_material_parameter_ptr)

  case ("hexagonal-linear")
    allocate(stiffness_model_hexagonal_temp)
    elasticity_base_ptr => stiffness_model_hexagonal_temp
    call stiffness_model_hexagonal_temp%linkBaseObjects(phase_id, all_mighty_grid_in, sim_all_macro_data, the_bc_object)
    call stiffness_model_hexagonal_temp%initParameters(phase_id, common_material_parameter_ptr)

  end select

  call elasticity_base_ptr%readMaterialParametersFromFile(matf_reader)

  select case(elasticity_model_name%getString())
  case ("isotropic-linear")
    nullify(stiffness_model_isotropic_temp)
  case ("isotropic-shear-poisson-linear")
    nullify(stiffness_model_isotropic_shear_temp)
  case ("cubic-linear")
    nullify(stiffness_model_cubic_temp)
  case ("hexagonal-linear")
    nullify(stiffness_model_hexagonal_temp)
  case default
    write(*,*) "somthing is wrong when trying nullifying the pointer for ", elasticity_model_name%getString()
    error stop "abort"
  end select
end subroutine

end module
