module material_base_mod
use kinds
use string_module, only : string_type
use all_grid_data_mod, only : all_grid_data
! use all_mighty_grid_mod, only : all_mighty_grid_type ! using the all mighty grid will cause a common block generation (WHY???)
! we need to dig into this, for now the material_base does not need to access the all mighty
use bc_objects_mod, only : boundary_condition_array_type
use simumaltion_macro_data_mod, only : sim_all_macro_data_obj
use polymorphic_dtype_array_mod, only : dtype_array_ptr
use common_material_parameter_mod, only : common_material_parameter

#include "macro_debug.fpp"

implicit none

type :: material_base
  ! type(all_mighty_grid_type), pointer :: all_mighty_grid => null()
  class(all_grid_data), pointer :: grid_data => null(), &
                                   common_grid_data => null()
  class(sim_all_macro_data_obj), pointer :: sim_all_macro_data => null()

  logical :: macro_object_linked =.false.
  logical :: grid_variables_provided =.false.
  logical :: grid_pointers_linked =.false.
  logical :: parameters_initialized =.false.
  logical :: state_variables_initialized =.false.
  logical :: write_properties_header =.true.
  integer :: phase_id, ix, iy, iz
  integer :: xs_rank, xe_rank, ys_rank, ye_rank, zs_rank, ze_rank
  real(k_real), pointer, dimension(:,:,:) :: phase_fraction_grid => null() ! a pointer to the material phase_id

  __DECL_CLASS_UNUSED_THIS__

  ! simulation parameter
  real(k_real), pointer :: dt => null(), &
                           temperature => null(), &
                           time => null(), &
                           dt_max => null()

  ! common material parameter pointer
  class(common_material_parameter), pointer :: common_material_parameter_ptr => null()

contains

  procedure :: linkBaseObjects => linkBaseObjectsMaterialBase

  procedure :: initParametersMaterialBase !-> base function initializing grid_data pointers and allocating the required space
  procedure :: initGridPointers=>initGridPointersMaterialBase !-> base function initializing grid_data pointers and allocating the required space

  ! setPointData must always be carryed over trough classes. This is used to set pointers for doing local calculations
  procedure :: setPointData => setPointersMaterialBase
  procedure :: addFieldVariablesToGrid => addFieldVariablesToGridMaterialBase

  procedure :: initStateVariables
  procedure :: initStateVariablesAtMaterialPoint => initStateVariablesAtMaterialPointMaterialBase

  procedure :: updateStateVariablesInnerLoop
  procedure :: updateStateVariablesAtMaterialPointInnerLoop => updateStateVariablesAtMaterialPointInnerLoopMaterialBase

  procedure :: updateStateVariablesOuterLoop
  procedure :: updateStateVariablesAtMaterialPointOuterLoop => updateStateVariablesAtMaterialPointOuterLoopMaterialBase

  procedure :: updateStateVariablesStaggered
  procedure :: updateStateVariablesAtMaterialPointStaggered => updateStateVariablesAtMaterialPointStaggeredMaterialBase

  procedure :: isInitialized => isInitializedMaterialBase

  procedure :: readMaterialParametersFromFile => readMaterialParametersFromFileMaterialBase

  procedure :: updatePhaseParameters => updatePhaseParametersMaterialBase

  procedure :: acceptRejectSolution => acceptRejectSolutionMaterialBase

  procedure :: writeAverageQuantitiesCSV => writeAverageQuantitiesCSVMaterialBase

end type material_base

contains
  subroutine linkBaseObjectsMaterialBase(this, phase_id, all_mighty_grid_in, sim_all_macro_data, the_bc_object)
    use all_mighty_grid_mod, only : all_mighty_grid_type
    implicit none
    class(material_base), intent(inout) :: this
    integer, intent(in) :: phase_id
    class(all_mighty_grid_type), intent(in), pointer :: all_mighty_grid_in
    class(sim_all_macro_data_obj), intent(in), pointer :: sim_all_macro_data
    class(boundary_condition_array_type), intent(in) :: the_bc_object

    if (this%macro_object_linked) error stop "linkBaseObjects already intialized. abort! "
    this%phase_id = phase_id

    ! this%all_mighty_grid => all_mighty_grid_in
    call all_mighty_grid_in%AMGgetCommonGridPtr(this%common_grid_data)
    call all_mighty_grid_in%AMGgetPhaseGridPtrByIdx(this%phase_id, this%grid_data)
    call all_mighty_grid_in%AMGGetLoopLimitsRank(this%xs_rank, this%xe_rank, this%ys_rank, this%ye_rank, this%zs_rank, this%ze_rank )
    this%sim_all_macro_data => sim_all_macro_data

    call the_bc_object%getTemperaturePointer(this%temperature)
    call sim_all_macro_data%sim_time%getTimePointer(this%time)
    call sim_all_macro_data%sim_time%getDeltaTimePointer(this%dt)
    call sim_all_macro_data%sim_time%getDeltaTMaxPointer(this%dt_max)

    this%macro_object_linked =.true.

  end subroutine

  subroutine writeAverageQuantitiesCSVMaterialBase(this, csv_writer_obj, write_headers)
    use csv_writer_mod, only : csv_writer
    implicit none
    class(material_base), intent(inout) :: this
    class(csv_writer), intent(inout) :: csv_writer_obj
    logical, intent(in) :: write_headers
__DECL_UNUSED_LOGICAL__
__DECL_CLASS_UNUSED_THIS__

__SUPPRESS_CLASS_UNUSED__(csv_writer_obj)
__SUPPRESS_CLASS_UNUSED_THIS__
__SUPPRESS_UNUSED_LOGICAL__(write_headers)
  end subroutine

  subroutine addFieldVariablesToGridMaterialBase(this)
    use grid_data_var_type_mod
    implicit none
    class(material_base), intent(inout) :: this

    if (.not.(this%macro_object_linked)) error stop "addFieldVariablesToGridMaterialBase: you can't add field varaibles to the grid without first linking the global obejcts"
    if ( .not.(this%parameters_initialized)) error stop "addFieldVariablesToGridMaterialBase: you can't add field varaibles to the grid without first initializing a material parameters"

    associate (all_grid_data_vars => this%grid_data)
    call all_grid_data_vars%addVar("phase_id", scalar_integer)
    end associate

    this%grid_variables_provided =.true.
  end subroutine

  subroutine initParametersMaterialBase(this, phase_id, common_material_parameter_ptr)
    implicit none
    class(material_base), intent(inout) :: this
    integer, intent(in) :: phase_id
    class(common_material_parameter), pointer, intent(in) :: common_material_parameter_ptr

    if (.not.this%macro_object_linked) error stop "initParametersMaterialBase base objects have not been linked"
    if (this%parameters_initialized) error stop "initParametersMaterialBase you can initialize parameters only once"

    this%phase_id = phase_id

    this%parameters_initialized = .true.
    this%common_material_parameter_ptr => common_material_parameter_ptr
  end subroutine

  subroutine initGridPointersMaterialBase(this)
    implicit none
    class(material_base), intent(inout) :: this
    real(k_real), pointer, dimension(:,:,:,:) :: temp_vector_pointer => null()
    if (.not.this%macro_object_linked) &
      error stop "initGridPointersMaterialBase you cannot initialized GridPointers without first linking base objects. Abort!"
    if (.not.this%grid_variables_provided) &
      error stop "initGridPointersMaterialBase you cannot initialized GridPointers without first providing the grid variables. Abort!"
    if (this%grid_pointers_linked) error stop "you can link grid pointers only once. Abort!"
    if (.not.(this%grid_data%initialized)) error stop "initGridPointersMaterialBase you cannot init grid pointer without first initializing the grid!"

    call this%common_grid_data%getGenericVectorDataPointerByName("phase_fraction", temp_vector_pointer)
    this%phase_fraction_grid => temp_vector_pointer(this%phase_id,:,:,:)
    nullify(temp_vector_pointer)
    this%grid_pointers_linked =.true.
  end subroutine

subroutine setPointersMaterialBase(this, ix, iy, iz)
  use change_tensor_basis, only : chg_basis_tensor2_to_vector6
  implicit none
  class(material_base), intent(inout) :: this
  integer, intent(in) ::  ix, iy, iz
  if (.not.(this%grid_data%initialized)) error stop "setPointersMaterialBase you cannot init grid pointer without first initializing the grid!"
  this%ix = ix; this%iy = iy; this%iz = iz;

end subroutine

subroutine initStateVariables(this)
  implicit none
  class(material_base), intent(inout) :: this

  if (.not.(this%parameters_initialized)) error stop "cannot initialize state variables without first initializing parameters "
  if (.not.(this%macro_object_linked)) error stop "cannot initialize state variables without first linking macro objects "
  if (.not.(this%grid_variables_provided)) error stop "cannot initialize state variables without first providing state variables "
  if (.not.(this%grid_pointers_linked)) error stop "cannot initialize state variables without first linking grid_data pointers "

  associate (ix=>this%ix, xs=>this%xs_rank, xe=>this%xe_rank, &
             iy=>this%iy, ys=>this%ys_rank, ye=>this%ye_rank, &
             iz=>this%iz, zs=>this%zs_rank, ze=>this%ze_rank)

  call this%updatePhaseParameters
  
  do iz=zs,ze
    do iy=ys,ye
      do ix=xs,xe
        call this%setPointData(ix,iy,iz)
        call this%initStateVariablesAtMaterialPoint()
      enddo
    enddo
  enddo

  this%state_variables_initialized = .true.

  end associate
end subroutine

subroutine initStateVariablesAtMaterialPointMaterialBase(this)
  class(material_base), intent(inout) :: this
  __DECL_CLASS_UNUSED_THIS__

  __SUPPRESS_CLASS_UNUSED_THIS__
  error stop "If you end up here you did not override initStateVariablesAtMaterialPoint in your class"
end subroutine


subroutine updateStateVariablesAtMaterialPointInnerLoopMaterialBase(this)
  implicit none
  class(material_base), intent(inout) :: this
  __DECL_CLASS_UNUSED_THIS__

  __SUPPRESS_CLASS_UNUSED_THIS__
  error stop "If you end up here you did not override updateStateVariablesAtMaterialPointInnerLoopMaterialBase in your class"
end subroutine

subroutine updateStateVariablesAtMaterialPointOuterLoopMaterialBase(this)
  implicit none
  class(material_base), intent(inout) :: this
  __DECL_CLASS_UNUSED_THIS__

  __SUPPRESS_CLASS_UNUSED_THIS__
  ! error stop "If you end up here you did not override updateStateVariablesAtMaterialPointOuterLoopMaterialBase in your class"
end subroutine

subroutine updateStateVariablesInnerLoop(this, ix_in, iy_in, iz_in)
  implicit none
  class(material_base), intent(inout) :: this
  integer, intent(in) :: ix_in, iy_in, iz_in

  if (.not.(associated(this%phase_fraction_grid))) error stop "initStateVariables phase_fraction_grid not associated"

  call this%updatePhaseParameters


  this%ix = ix_in
  this%iy = iy_in
  this%iz = iz_in
  call this%setPointData(ix_in, iy_in, iz_in)
  if (this%phase_fraction_grid(ix_in, iy_in, iz_in).gt.0._k_real) then
  call this%updateStateVariablesAtMaterialPointInnerLoop()
  endif

end subroutine

subroutine updateStateVariablesOuterLoop(this)
  implicit none
  class(material_base), intent(inout) :: this

  associate (ix=>this%ix, xs=>this%xs_rank, xe=>this%xe_rank, &
            iy=>this%iy, ys=>this%ys_rank, ye=>this%ye_rank, &
            iz=>this%iz, zs=>this%zs_rank, ze=>this%ze_rank)

    do iz=zs,ze
      do iy=ys,ye
        do ix=xs,xe
          call this%setPointData(ix,iy,iz)
          if (this%phase_fraction_grid(ix,iy,iz).gt.0._k_real) then
          call this%updateStateVariablesAtMaterialPointOuterLoop()
          endif
        enddo
      enddo
    enddo

  end associate

end subroutine

subroutine updateStateVariablesAtMaterialPointStaggeredMaterialBase(this)
  implicit none
  class(material_base), intent(inout) :: this
  __DECL_CLASS_UNUSED_THIS__

  __SUPPRESS_CLASS_UNUSED_THIS__
  error stop "If you end up here you did not override updateStateVariablesAtMaterialPointStaggeredMaterialBase in your class"
end subroutine


subroutine updateStateVariablesStaggered(this)
  implicit none
  class(material_base), intent(inout) :: this

  call this%updatePhaseParameters


  associate (ix=>this%ix, xs=>this%xs_rank, xe=>this%xe_rank, &
            iy=>this%iy, ys=>this%ys_rank, ye=>this%ye_rank, &
            iz=>this%iz, zs=>this%zs_rank, ze=>this%ze_rank)


    do iz=zs,ze
      do iy=ys,ye
        do ix=xs,xe
          call this%setPointData(ix,iy,iz)
          if (this%phase_fraction_grid(ix,iy,iz).gt.0._k_real) then
          call this%updateStateVariablesAtMaterialPointStaggered()
          endif
        enddo
      enddo
    enddo

  end associate

end subroutine

subroutine updatePhaseParametersMaterialBase(this)
  class(material_base), intent(inout) :: this
  __DECL_CLASS_UNUSED_THIS__

  __SUPPRESS_CLASS_UNUSED_THIS__
end subroutine

subroutine readMaterialParametersFromFileMaterialBase(this, matf_reader)
  use read_from_file_utils, only : file_reader
  implicit none
  class(material_base), intent(inout) :: this
  class(file_reader), intent(inout) :: matf_reader
  __DECL_CLASS_UNUSED_THIS__
  __DECL_UNUSED_FILE_READER__
 
  __SUPPRESS_CLASS_UNUSED_THIS__
  __SUPPRESS_UNUSED_FILE_READER__(matf_reader)
  error stop "If you end up here you did not override readMaterialParametersFromFileMaterialBase in your class"
end

subroutine isInitializedMaterialBase(this)
  class(material_base), intent(inout) :: this

  if (this%macro_object_linked.and. &
      this%grid_pointers_linked.and. &
      this%parameters_initialized.and. &
      this%state_variables_initialized) then
    return
  else
    if (.not.(this%macro_object_linked)) error stop "Macro Object not linked"
    if (.not.(this%grid_variables_provided)) error stop "Grid variables have not provided"
    if (.not.(this%grid_pointers_linked)) error stop "Grid Pointers not linked"
    if (.not.(this%parameters_initialized)) error stop "Parameters not initialized"
    if (.not.(this%state_variables_initialized)) error stop "State Variables not initialized"
  endif
end subroutine

subroutine acceptRejectSolutionMaterialBase(this, dt_max, accept_solution_flag)
  implicit none
  class(material_base), intent(inout) :: this
  real(k_real), intent(out) :: dt_max
  logical, intent(out) :: accept_solution_flag
  __DECL_CLASS_UNUSED_THIS__

  __SUPPRESS_CLASS_UNUSED_THIS__
  __SUPPRESS_UNUSED_REAL_OUT__(dt_max)
  __SUPPRESS_UNUSED_LOGICAL_OUT__(accept_solution_flag)
  error stop "If you end up here it means you did not override acceptRejectSolutionMaterialBase in your class"
end subroutine

end module material_base_mod
