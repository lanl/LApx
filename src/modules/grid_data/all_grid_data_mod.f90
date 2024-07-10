module all_grid_data_mod
#include "macro_debug.fpp"
  use grid_data_types
  use grid_data_var_type_mod
  use string_module
  use kinds
  use polymorphic_dtype_array_mod
  use grid_data_types, only : var_dimension_type
  implicit none


  type all_grid_data
    logical :: initialized = .FALSE.
    type(griddata_scalar), pointer, dimension(:) :: all_scalar =>null()
    type(griddata_scalar_integer), pointer, dimension(:) :: all_scalar_integer =>null()
    type(griddata_real_vector), pointer, dimension(:) :: all_real_vector =>null()
    type(griddata_vector5), pointer, dimension(:) :: all_vector5 =>null()
    type(griddata_vector6), pointer, dimension(:) :: all_vector6 =>null()
    type(griddata_generic_vector), pointer, dimension(:) :: all_generic_vector =>null()
    type(griddata_ssscalar), pointer, dimension(:) :: all_ss_scalar =>null()
    type(griddata_tensor2), pointer, dimension(:) :: all_tensor2 =>null()
    type(griddata_matrix66), pointer, dimension(:) :: all_matrix66 =>null()
    type(griddata_ss_generic_vector), pointer, dimension(:) :: all_ss_generic_vector =>null()
    type(griddata_ss_vector5), pointer, dimension(:) :: all_ss_vector5 =>null()
    type(griddata_ss_vector6), pointer, dimension(:) :: all_ss_vector6 =>null()
    type(griddata_ss_generic_matrix), pointer, dimension(:) :: all_ss_generic_matrix =>null()
    type(griddata_generic_matrix), pointer, dimension(:) :: all_generic_matrix =>null()
    class(var_dimension_type), pointer :: dimension_obj => null()
    integer, dimension(3) :: nx_ny_nz=(/0,0,0/), &
                             nx_ny_nz_rank =(/0,0,0/), &
                             xyz_start_rank=(/0,0,0/), &
                             xyz_end_rank=(/0,0,0/), &
                             xyz_offset_rank=(/0,0,0/)
    integer :: x_start_rank, x_end_rank, x_offset_rank
    integer :: y_start_rank, y_end_rank, y_offset_rank
    integer :: z_start_rank, z_end_rank, z_offset_rank
    integer :: n_points = 0
    integer :: n_points_rank = 0
    integer :: n_ss = 0
    integer :: phase_id
    real(k_real) :: voxel_weight = 0._k_real
    logical :: all_scalar_has_data
    logical :: all_scalar_integer_has_data
    logical :: all_real_vector_has_data
    logical :: all_vector5_has_data
    logical :: all_vector6_has_data
    logical :: all_generic_vector_has_data
    logical :: all_ss_scalar_has_data
    logical :: all_tensor2_has_data
    logical :: all_matrix66_has_data
    logical :: all_ss_generic_vector_has_data
    logical :: all_ss_vector5_has_data
    logical :: all_ss_vector6_has_data
    logical :: all_ss_matrix_has_data
    logical :: all_generic_matrix_has_data
    logical :: num_ss_was_set = .false.
    logical :: dimensions_were_set = .false.
    type(grid_data_var_type_array) :: all_grid_vars_list

    !! dump for restart parameters
    integer :: n_file_to_keep = 0
    integer :: dump_every_n_steps = 100 
    type(string_array) :: restart_file_names, xmf_restart_file_names
    type(string_type) :: restart_file_base_name
    type(string_type) :: phase_group_name
    type(string_type) :: field_absolute_path
    integer :: write_fields_every_n_steps = 0
    type(string_type) :: field_file_base_name
    __DECL_CLASS_UNUSED_THIS__


  contains
    procedure :: init => initAGD
    procedure :: setDimensions => setDimensionsAGD
    procedure :: getGlobalGridDimension
    procedure :: getGlobalGridNpoints
    procedure :: getRankGridDimension
    procedure :: getRankZStart
    procedure :: getRankZEnd
    procedure :: getZRankFromZGlobal
    procedure :: getZGlobalFromZRank
    procedure :: AGDupdateAllStatefulVariables
    procedure :: AGDResetAllStatefulVariables
    procedure :: setDumpForRestartOptions
    procedure :: setWriteFieldOptions
    procedure :: AGDDumpForRestart
    procedure :: AGDReloadFromDump
    procedure :: DumpMaterialPointValuesToTextFile
    procedure :: addVar
    procedure :: setNumSlipSystems => setNumSlipSystemsAGD

    !************************************************************!
    !***********************SCALAR DATA API**********************!
    !************************************************************!
    procedure, private :: AGDGetGridDataScalarPointerByType
    procedure, private :: AGDGetScalarDataPointerByNameTypeAndStatefulIndex
    procedure, private :: AGDGetAvgScalarDataPointerByNameTypeAndStatefulIndex
    procedure, private :: AGDGetGridDataScalarVariablePointerByName

    procedure :: getScalarDataPointerByName => AGDGetScalarDataPointerByName
    procedure :: getScalarDataPointerByNameOld => AGDGetScalarDataPointerByNameOld
    procedure :: getScalarDataPointerByNameOlder => AGDGetScalarDataPointerByNameOlder
    procedure :: getAvgScalarDataPointerByName => AGDGetAvgScalarDataPointerByName
    procedure :: getAvgScalarDataPointerByNameOld => AGDGetAvgScalarDataPointerByNameOld
    procedure :: getAvgScalarDataPointerByNameOlder => AGDGetAvgScalarDataPointerByNameOlder
    procedure :: getScalarPointerToVariableByName => AGDGetScalarPointerToVariableByName

    ! this procedure can pick a pointer from any integer scalar field
    procedure, private :: AGDGetScalarIntegerDataPointerByNameTypeAndStatefulIndex
    procedure :: getScalarIntegerDataPointerByName => AGDGetScalarIntegerDataPointerByName
    procedure :: getScalarIntegerDataPointerByNameOld => AGDGetScalarIntegerDataPointerByNameOld
    procedure :: getScalarIntegerDataPointerByNameOlder => AGDGetScalarIntegerDataPointerByNameOlder

    !************************************************************!
    !***********************VECTOR DATA API**********************!
    !************************************************************!
    procedure, private :: AGDGetGridDataVectorPointerByType
    procedure, private :: AGDGetVectorDataPointerByNameTypeAndStatefulIndex
    procedure, private :: AGDGetAvgVectorDataPointerByNameTypeAndStatefulIndex
    procedure, private :: AGDGetGridDataVectorVariablePointerByName

    procedure :: getRealVectorDataPointerByName => AGDGetRealVectorDataPointerByName
    procedure :: getRealVectorDataPointerByNameOld => AGDGetRealVectorDataPointerByNameOld
    procedure :: getRealVectorDataPointerByNameOlder => AGDGetRealVectorDataPointerByNameOlder
    procedure :: getRealAvgVectorDataPointerByName => AGDGetAvgRealVectorDataPointerByName
    procedure :: getRealAvgVectorDataPointerByNameOld => AGDGetAvgRealVectorDataPointerByNameOld
    procedure :: getRealAvgVectorDataPointerByNameOlder => AGDGetAvgRealVectorDataPointerByNameOlder

    procedure :: getVector5DataPointerByName => AGDGetVector5DataPointerByName
    procedure :: getVector5DataPointerByNameOld => AGDGetVector5DataPointerByNameOld
    procedure :: getVector5DataPointerByNameOlder => AGDGetVector5DataPointerByNameOlder
    procedure :: getAvgVector5DataPointerByName => AGDGetAvgVector5DataPointerByName
    procedure :: getAvgVector5DataPointerByNameOld => AGDGetAvgVector5DataPointerByNameOld
    procedure :: getAvgVector5DataPointerByNameOlder => AGDGetAvgVector5DataPointerByNameOlder

    procedure :: getVector6DataPointerByName => AGDGetVector6DataPointerByName
    procedure :: getVector6DataPointerByNameOld => AGDGetVector6DataPointerByNameOld
    procedure :: getVector6DataPointerByNameOlder => AGDGetVector6DataPointerByNameOlder
    procedure :: getAvgVector6DataPointerByName => AGDGetAvgVector6DataPointerByName
    procedure :: getAvgVector6DataPointerByNameOld => AGDGetAvgVector6DataPointerByNameOld
    procedure :: getAvgVector6DataPointerByNameOlder => AGDGetAvgVector6DataPointerByNameOlder
    procedure :: getVector6PointerToVariableByName => AGDGetVector6PointerToVariableByName

    procedure :: getGenericVectorDataPointerByName => AGDGetGenericVectorDataPointerByName
    procedure :: getGenericVectorDataPointerByNameOld => AGDGetGenericVectorDataPointerByNameOld
    procedure :: getGenericVectorDataPointerByNameOlder => AGDGetGenericVectorDataPointerByNameOlder
    procedure :: getAvgGenericVectorDataPointerByName => AGDGetAvgGenericVectorDataPointerByName
    procedure :: getAvgGenericVectorDataPointerByNameOld => AGDGetAvgGenericVectorDataPointerByNameOld
    procedure :: getAvgGenericVectorDataPointerByNameOlder => AGDGetAvgGenericVectorDataPointerByNameOlder

    procedure :: getSSScalarDataPointerByName => AGDgetSSScalarDataPointerByName
    procedure :: getSSScalarDataPointerByNameOld => AGDgetSSScalarDataPointerByNameOld
    procedure :: getSSScalarDataPointerByNameOlder => AGDgetSSScalarDataPointerByNameOlder
    procedure :: getAvgSSScalarDataPointerByName => AGDgetAvgSSScalarDataPointerByName
    procedure :: getAvgSSScalarDataPointerByNameOld => AGDgetAvgSSScalarDataPointerByNameOld
    procedure :: getAvgSSScalarDataPointerByNameOlder => AGDgetAvgSSScalarDataPointerByNameOlder

    !************************************************************!
    !***********************MATRIX DATA API**********************!
    !************************************************************!
    procedure, private :: AGDGetGridDataMatrixPointerByType
    procedure, private :: AGDGetMatrixDataPointerByNameTypeAndStatefulIndex
    procedure, private :: AGDGetAvgMatrixDataPointerByNameTypeAndStatefulIndex
    procedure, private :: AGDGetGridDataMatrixVariablePointerByName

    procedure :: getGenericMatrixDataPointerByName => AGDGetGenericMatrixDataPointerByName
    procedure :: getGenericMatrixDataPointerByNameOld => AGDGetGenericMatrixDataPointerByNameOld
    procedure :: getGenericMatrixDataPointerByNameOlder => AGDGetGenericMatrixDataPointerByNameOlder
    procedure :: getAvgGenericMatrixDataPointerByName => AGDGetAvgGenericMatrixDataPointerByName
    procedure :: getAvgGenericMatrixDataPointerByNameOld => AGDGetAvgGenericMatrixDataPointerByNameOld
    procedure :: getAvgGenericMatrixDataPointerByNameOlder => AGDGetAvgGenericMatrixDataPointerByNameOlder

    procedure :: getMatrix66DataPointerByName => AGDGetMatrix66DataPointerByName
    procedure :: getMatrix66DataPointerByNameOld => AGDGetMatrix66DataPointerByNameOld
    procedure :: getMatrix66DataPointerByNameOlder => AGDGetMatrix66DataPointerByNameOlder
    procedure :: getAvgMatrix66DataPointerByName => AGDGetAvgMatrix66DataPointerByName
    procedure :: getAvgMatrix66DataPointerByNameOld => AGDGetAvgMatrix66DataPointerByNameOld
    procedure :: getAvgMatrix66DataPointerByNameOlder => AGDGetAvgMatrix66DataPointerByNameOlder

    procedure :: getTensor2DataPointerByName => AGDGetTensor2DataPointerByName
    procedure :: getTensor2DataPointerByNameOld => AGDGetTensor2DataPointerByNameOld
    procedure :: getTensor2DataPointerByNameOlder => AGDGetTensor2DataPointerByNameOlder
    procedure :: getAvgTensor2DataPointerByName => AGDGetAvgTensor2DataPointerByName
    procedure :: getAvgTensor2DataPointerByNameOld => AGDGetAvgTensor2DataPointerByNameOld
    procedure :: getAvgTensor2DataPointerByNameOlder => AGDGetAvgTensor2DataPointerByNameOlder
    procedure :: getTensor2PointerToVariableByName => AGDGetTensor2PointerToVariableByName

    procedure :: getSSGenericVectorDataPointerByName => AGDGetSSGenericVectorDataPointerByName
    procedure :: getSSGenericVectorDataPointerByNameOld => AGDGetSSGenericVectorDataPointerByNameOld
    procedure :: getSSGenericVectorDataPointerByNameOlder => AGDGetSSGenericVectorDataPointerByNameOlder
    procedure :: getAvgSSGenericVectorDataPointerByName => AGDGetAvgSSGenericVectorDataPointerByName
    procedure :: getAvgSSGenericVectorDataPointerByNameOld => AGDGetAvgSSGenericVectorDataPointerByNameOld
    procedure :: getAvgSSGenericVectorDataPointerByNameOlder => AGDGetAvgSSGenericVectorDataPointerByNameOlder

    procedure :: getSSVector5DataPointerByName => AGDGetSSVector5DataPointerByName
    procedure :: getSSVector5DataPointerByNameOld => AGDGetSSVector5DataPointerByNameOld
    procedure :: getSSVector5DataPointerByNameOlder => AGDGetSSVector5DataPointerByNameOlder
    procedure :: getAvgSSVector5DataPointerByName => AGDGetAvgSSVector5DataPointerByName
    procedure :: getAvgSSVector5DataPointerByNameOld => AGDGetAvgSSVector5DataPointerByNameOld
    procedure :: getAvgSSVector5DataPointerByNameOlder => AGDGetAvgSSVector5DataPointerByNameOlder

    procedure :: getSSVector6DataPointerByName => AGDGetSSVector6DataPointerByName
    procedure :: getSSVector6DataPointerByNameOld => AGDGetSSVector6DataPointerByNameOld
    procedure :: getSSVector6DataPointerByNameOlder => AGDGetSSVector6DataPointerByNameOlder
    procedure :: getSSAvgVector6DataPointerByName => AGDGetAvgSSVector6DataPointerByName
    procedure :: getSSAvgVector6DataPointerByNameOld => AGDGetAvgSSVector6DataPointerByNameOld
    procedure :: getSSAvgVector6DataPointerByNameOlder => AGDGetAvgSSVector6DataPointerByNameOlder

    !************************************************************!
    !***********************SS MATRIX DATA API**********************!
    !************************************************************!

    procedure, private :: AGDGetGridDataSSMatrixPointerByType
    procedure, private :: AGDGetSSMatrixDataPointerByNameTypeAndStatefulIndex
    procedure, private :: AGDGetAvgSSMatrixDataPointerByNameTypeAndStatefulIndex

    procedure :: getSSMatrixDataPointerByName => AGDGetSSMatrixDataPointerByName
    procedure :: getSSMatrixDataPointerByNameOld => AGDGetSSMatrixDataPointerByNameOld
    procedure :: getSSMatrixDataPointerByNameOlder => AGDGetSSMatrixDataPointerByNameOlder
    procedure :: getSSAvgMatrixDataPointerByName => AGDGetAvgSSMatrixDataPointerByName
    procedure :: getSSAvgMatrixDataPointerByNameOld => AGDGetAvgSSMatrixDataPointerByNameOld
    procedure :: getSSAvgMatrixDataPointerByNameOlder => AGDGetAvgSSMatrixDataPointerByNameOlder
  end type

  type, extends(dtype_array_ptr) :: all_grid_data_multi_phase
  contains
    procedure :: addElement => addGridData
    procedure :: getElementPtr => getGridDataPtr
    ! procedure :: setDimensions => GridDatasetDimensions
  end type

interface

  module subroutine DumpMaterialPointValuesToTextFile(this, ix,iy,iz, file_id)
    implicit none
    class(all_grid_data), intent(inout) :: this
    integer, intent(in) :: ix,iy,iz, file_id
  end subroutine

  module subroutine initAGD(this, phase_id)
    implicit none
    class(all_grid_data), intent(inout) :: this
    integer, intent(in), optional :: phase_id
  end subroutine

  module subroutine  setDimensionsAGD(this, dimension_obj)
    implicit none
    class(all_grid_data), intent(inout) :: this
    class(var_dimension_type), intent(in), pointer :: dimension_obj
    ! integer, intent(in) :: n_points_global(3), &
    !                         xyz_start_rank(3), &
    !                         xyz_end_rank(3), &
    !                         xyz_offset_rank(3)
  end subroutine

  module subroutine  setNumSlipSystemsAGD(this, nss)
    class(all_grid_data), intent(inout) :: this
    integer, intent(in) :: nss
  end subroutine

  module subroutine AGDDumpForRestart(this, file_name, xmf_writer, write_only_current_value)
    use xmf_writer_mod, only :xmf_writer_type
    implicit none
    class(all_grid_data), intent(inout) :: this
    character(len=*), intent(in) :: file_name
    type(xmf_writer_type), intent(inout) :: xmf_writer
    logical, intent(in) :: write_only_current_value
  end subroutine

  module subroutine AGDReloadFromDump(this, file_name, time, dt, step, XLSEC, XLSEC_old)
    implicit none
    class(all_grid_data), intent(inout) :: this
    character(len=*), intent(in) :: file_name
    real(k_real), intent(out) :: time, dt
    integer, intent(out) :: step
    real(k_real), intent(out) :: XLSEC(6,6), XLSEC_old(6,6)
  end subroutine

  module function getGlobalGridNpoints(this) result(npoints)
    implicit none
    class(all_grid_data), intent(in) :: this
    integer :: npoints
  end function

  module function getGlobalGridDimension(this) result(nx_ny_nz)
    implicit none
    class(all_grid_data), intent(in) :: this
    integer, dimension(3) :: nx_ny_nz
  end function

  module function getRankGridDimension(this) result(nx_ny_nz)
    implicit none
    class(all_grid_data), intent(in) :: this
    integer, dimension(3) :: nx_ny_nz
  end function

  module function getRankZStart(this) result(z_start_rank)
    implicit none
    class(all_grid_data), intent(in) :: this
    integer :: z_start_rank
  end function

  module function getRankZEnd(this) result(z_end_rank)
    implicit none
    class(all_grid_data), intent(in) :: this
    integer :: z_end_rank
  end function

  module function getZRankFromZGlobal(this, z_global) result(z_local)
    implicit none
    class(all_grid_data), intent(in) :: this
    integer, intent(in):: z_global
    integer :: z_local
  end function

  module function getZGlobalFromZRank(this, z_local) result(z_global)
    implicit none
    class(all_grid_data), intent(in) :: this
    integer, intent(in):: z_local
    integer :: z_global
  end function

  module subroutine AGDupdateAllStatefulVariables(this)
    implicit none
    class(all_grid_data), intent(inout) :: this
  end subroutine

  module subroutine AGDResetAllStatefulVariables(this)
    implicit none
    class(all_grid_data), intent(inout) :: this
  end subroutine

  !************************************************************!
  !**************SCALAR DATA ROUTINE PROTOTYPES ***************!
  !************************************************************!

  module subroutine AGDGetGridDataScalarPointerByType(this, scalar_type, generic_griddata_scalar_ptr, has_data)
      implicit none
      class(all_grid_data), intent(in) :: this
      integer(kind(grid_data_var_type_enum)), intent(in) :: scalar_type
      class(griddata_scalar), dimension(:), pointer, intent(out):: generic_griddata_scalar_ptr
      logical, intent(out) :: has_data
  end subroutine

  module subroutine AGDGetScalarDataPointerByNameTypeAndStatefulIndex(this, scalar_name, scalar_type, stateful_index, pointer_2_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: scalar_name
    integer(kind(grid_data_var_type_enum)), intent(in) :: scalar_type
    integer, intent(in) :: stateful_index
    real(k_real), dimension(:,:,:), pointer, intent(inout) :: pointer_2_data
  end subroutine

  module subroutine AGDGetAvgScalarDataPointerByNameTypeAndStatefulIndex(this, scalar_name, scalar_type, stateful_index, pointer_2_avg_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: scalar_name
    integer(kind(grid_data_var_type_enum)), intent(in) :: scalar_type
    integer, intent(in) :: stateful_index
    real(k_real), pointer, intent(inout) :: pointer_2_avg_data
  end subroutine

  module subroutine AGDGetGridDataScalarVariablePointerByName(this, scalar_name, scalar_type, pointer_2_variable)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: scalar_name
    integer(kind(grid_data_var_type_enum)), intent(in) :: scalar_type
    class(griddata_scalar), pointer, intent(inout) :: pointer_2_variable
  end subroutine

  !************************************************************!
  !*********SCALAR INTEGER DATA ROUTINE PROTOTYPES ************!
  !************************************************************!

  module subroutine AGDGetScalarIntegerDataPointerByNameTypeAndStatefulIndex(this, scalar_integer_name, scalar_integer_type, stateful_index, pointer_2_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: scalar_integer_name
    integer(kind(grid_data_var_type_enum)), intent(in) :: scalar_integer_type
    integer, intent(in) :: stateful_index
    integer(k_int), dimension(:,:,:), pointer, intent(inout) :: pointer_2_data
  end subroutine

  !************************************************************!
  !**************VECTOR DATA ROUTINE PROTOTYPES ***************!
  !************************************************************!

  module subroutine AGDGetGridDataVectorPointerByType(this, vector_type, generic_griddata_vector_ptr, has_data)
      implicit none
      class(all_grid_data), intent(in) :: this
      integer(kind(grid_data_var_type_enum)), intent(in) :: vector_type
      class(griddata_vector), dimension(:), pointer, intent(out):: generic_griddata_vector_ptr
      logical, intent(out) :: has_data
  end subroutine

  module subroutine AGDGetVectorDataPointerByNameTypeAndStatefulIndex(this, vector_name, vector_type, stateful_index, pointer_2_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: vector_name
    integer(kind(grid_data_var_type_enum)), intent(in) :: vector_type
    integer, intent(in) :: stateful_index
    real(k_real), dimension(:,:,:,:), pointer, intent(inout) :: pointer_2_data
  end subroutine

  module subroutine AGDGetAvgVectorDataPointerByNameTypeAndStatefulIndex(this, vector_name, vector_type, stateful_index, pointer_2_avg_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: vector_name
    integer(kind(grid_data_var_type_enum)), intent(in) :: vector_type
    integer, intent(in) :: stateful_index
    real(k_real), dimension(:), pointer, intent(inout) :: pointer_2_avg_data
  end subroutine

  module subroutine AGDGetGridDataVectorVariablePointerByName(this, vector_name, vector_type, pointer_2_variable)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: vector_name
    integer(kind(grid_data_var_type_enum)), intent(in) :: vector_type
    class(griddata_vector), pointer, intent(inout) :: pointer_2_variable
  end subroutine

  !************************************************************!
  !**************MATRIX DATA ROUTINE PROTOTYPES ***************!
  !************************************************************!

  module subroutine AGDGetGridDataMatrixPointerByType(this, matrix_type, generic_griddata_matrix_ptr, has_data)
      implicit none
      class(all_grid_data), intent(in) :: this
      integer(kind(grid_data_var_type_enum)), intent(in) :: matrix_type
      class(griddata_matrix), dimension(:), pointer, intent(out):: generic_griddata_matrix_ptr
      logical, intent(out) :: has_data
  end subroutine

  module subroutine AGDGetMatrixDataPointerByNameTypeAndStatefulIndex(this, matrix_name, matrix_type, stateful_index, pointer_2_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: matrix_name
    integer(kind(grid_data_var_type_enum)), intent(in) :: matrix_type
    integer, intent(in) :: stateful_index
    real(k_real), dimension(:,:,:,:,:), pointer, intent(inout) :: pointer_2_data
  end subroutine

  module subroutine AGDGetAvgMatrixDataPointerByNameTypeAndStatefulIndex(this, matrix_name, matrix_type, stateful_index, pointer_2_avg_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: matrix_name
    integer(kind(grid_data_var_type_enum)), intent(in) :: matrix_type
    integer, intent(in) :: stateful_index
    real(k_real), dimension(:,:), pointer, intent(inout) :: pointer_2_avg_data
  end subroutine

  module subroutine AGDGetGridDataMatrixVariablePointerByName(this, matrix_name, matrix_type, pointer_2_variable)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: matrix_name
    integer(kind(grid_data_var_type_enum)), intent(in) :: matrix_type
    class(griddata_matrix), pointer, intent(inout) :: pointer_2_variable
  end subroutine

  !***************************************************************!
  !**************SS MATRIX DATA ROUTINE PROTOTYPES ***************!
  !***************************************************************!

  module subroutine AGDGetGridDataSSMatrixPointerByType(this, ss_matrix_type, generic_griddata_ss_matrix_ptr, has_data)
      implicit none
      class(all_grid_data), intent(in) :: this
      integer(kind(grid_data_var_type_enum)), intent(in) :: ss_matrix_type
      class(griddata_ss_generic_matrix), dimension(:), pointer, intent(out):: generic_griddata_ss_matrix_ptr
      logical, intent(out) :: has_data
  end subroutine

  module subroutine AGDGetSSMatrixDataPointerByNameTypeAndStatefulIndex(this, ss_matrix_name, ss_matrix_type, stateful_index, pointer_2_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: ss_matrix_name
    integer(kind(grid_data_var_type_enum)), intent(in) :: ss_matrix_type
    integer, intent(in) :: stateful_index
    real(k_real), dimension(:,:,:,:,:,:), pointer, intent(inout) :: pointer_2_data
  end subroutine

  module subroutine AGDGetAvgSSMatrixDataPointerByNameTypeAndStatefulIndex(this, ss_matrix_name, ss_matrix_type, stateful_index, pointer_2_avg_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: ss_matrix_name
    integer(kind(grid_data_var_type_enum)), intent(in) :: ss_matrix_type
    integer, intent(in) :: stateful_index
    real(k_real), dimension(:,:,:), pointer, intent(inout) :: pointer_2_avg_data
  end subroutine
end interface

contains

  subroutine setDumpForRestartOptions(this, n_file_to_keep, dump_every_n_steps, restart_file_base_name )
    implicit none
    class(all_grid_data), intent(inout) :: this
    integer, intent(in) :: n_file_to_keep, dump_every_n_steps
    character(len=*) :: restart_file_base_name

    this%n_file_to_keep=n_file_to_keep
    this%dump_every_n_steps= dump_every_n_steps
    call this%restart_file_base_name%setString(restart_file_base_name)

  end subroutine

  subroutine setWriteFieldOptions(this, write_fields_every_n_steps, field_file_base_name )
    use number_to_string_mod, only : int2string
    implicit none
    class(all_grid_data), intent(inout) :: this
    integer, intent(in) :: write_fields_every_n_steps
    character(len=*) :: field_file_base_name

    this%write_fields_every_n_steps= write_fields_every_n_steps
    call this%field_absolute_path%setString("/")
    if (this%phase_id.gt.0) then
      call this%phase_group_name%setString("phase_"//trim(adjustl(int2string(this%phase_id))))
      call this%field_absolute_path%setString("/"//trim(adjustl((this%phase_group_name%getString())))//"/")
    else 
      call this%phase_group_name%resetString()
    endif 
    call this%field_file_base_name%setString(field_file_base_name)

  end subroutine

  !************************************************************!
  !***********************SCALAR DATA API *********************!
  !************************************************************!

  ! get pointer to data
  subroutine AGDGetScalarDataPointerByName(this, scalar_name, pointer_2_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: scalar_name
    real(k_real), dimension(:,:,:), pointer, intent(inout) :: pointer_2_data
    call this%AGDGetScalarDataPointerByNameTypeAndStatefulIndex(scalar_name, scalar, 1, pointer_2_data)
  end subroutine

  subroutine AGDGetScalarDataPointerByNameOld(this, scalar_name, pointer_2_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: scalar_name
    real(k_real), dimension(:,:,:), pointer, intent(inout) :: pointer_2_data
    call this%AGDGetScalarDataPointerByNameTypeAndStatefulIndex(scalar_name, scalar, 2, pointer_2_data)
  end subroutine

  subroutine AGDGetScalarDataPointerByNameOlder(this, scalar_name, pointer_2_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: scalar_name
    real(k_real), dimension(:,:,:), pointer, intent(inout) :: pointer_2_data
    call this%AGDGetScalarDataPointerByNameTypeAndStatefulIndex(scalar_name, scalar, 3, pointer_2_data)
  end subroutine

  ! get pointer to average data
  subroutine AGDGetAvgScalarDataPointerByName(this, scalar_name, pointer_2_avg_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: scalar_name
    real(k_real), pointer, intent(inout) :: pointer_2_avg_data
    call this%AGDGetAvgScalarDataPointerByNameTypeAndStatefulIndex(scalar_name, scalar, 1, pointer_2_avg_data)
  end subroutine

  subroutine AGDGetAvgScalarDataPointerByNameOld(this, scalar_name, pointer_2_avg_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: scalar_name
    real(k_real), pointer, intent(inout) :: pointer_2_avg_data
    call this%AGDGetAvgScalarDataPointerByNameTypeAndStatefulIndex(scalar_name, scalar, 2, pointer_2_avg_data)
  end subroutine

  subroutine AGDGetAvgScalarDataPointerByNameOlder(this, scalar_name, pointer_2_avg_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: scalar_name
    real(k_real), pointer, intent(inout) :: pointer_2_avg_data
    call this%AGDGetAvgScalarDataPointerByNameTypeAndStatefulIndex(scalar_name, scalar, 3, pointer_2_avg_data)
  end subroutine

  ! get pointer to variable
  subroutine AGDGetScalarPointerToVariableByName(this, scalar_name, pointer_2_variable)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: scalar_name
    class(griddata_scalar), pointer, intent(inout) :: pointer_2_variable
    class(griddata_scalar), pointer :: temp_pointer_to_variable => null()

    nullify(temp_pointer_to_variable)

    call this%AGDGetGridDataScalarVariablePointerByName(scalar_name, scalar, temp_pointer_to_variable)

    select type(temp_pointer_to_variable)
    class is (griddata_scalar)
        pointer_2_variable => temp_pointer_to_variable
    class default
        error stop "AGDGetScalarPointerToVariableByName wrong class"
    end select

  end subroutine

  !************************************************************!
  !*****************SCALAR INTEGER DATA API *******************!
  !************************************************************!

  ! get pointer to data
  subroutine AGDGetScalarIntegerDataPointerByName(this, scalar_integer_name, pointer_2_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: scalar_integer_name
    integer(k_int), dimension(:,:,:), pointer, intent(inout) :: pointer_2_data
    call this%AGDGetScalarIntegerDataPointerByNameTypeAndStatefulIndex(scalar_integer_name, scalar_integer, 1, pointer_2_data)
  end subroutine

  subroutine AGDGetScalarIntegerDataPointerByNameOld(this, scalar_integer_name, pointer_2_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: scalar_integer_name
    integer(k_int), dimension(:,:,:), pointer, intent(inout) :: pointer_2_data
    call this%AGDGetScalarIntegerDataPointerByNameTypeAndStatefulIndex(scalar_integer_name, scalar_integer,  2, pointer_2_data)
  end subroutine

  subroutine AGDGetScalarIntegerDataPointerByNameOlder(this, scalar_integer_name, pointer_2_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: scalar_integer_name
    integer(k_int), dimension(:,:,:), pointer, intent(inout) :: pointer_2_data
    call this%AGDGetScalarIntegerDataPointerByNameTypeAndStatefulIndex(scalar_integer_name, scalar_integer,  3, pointer_2_data)
  end subroutine

  !************************************************************!
  !***************REAL VECTOR DATA API ************************!
  !************************************************************!

  ! get pointer to data
  subroutine AGDGetRealVectorDataPointerByName(this, real_vector_name, pointer_2_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: real_vector_name
    real(k_real), dimension(:,:,:,:), pointer, intent(inout) :: pointer_2_data
    call this%AGDGetVectorDataPointerByNameTypeAndStatefulIndex(real_vector_name, real_vector, 1, pointer_2_data)
  end subroutine

  subroutine AGDGetRealVectorDataPointerByNameOld(this, real_vector_name, pointer_2_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: real_vector_name
    real(k_real), dimension(:,:,:,:), pointer, intent(inout) :: pointer_2_data
    call this%AGDGetVectorDataPointerByNameTypeAndStatefulIndex(real_vector_name, real_vector, 2, pointer_2_data)
  end subroutine

  subroutine AGDGetRealVectorDataPointerByNameOlder(this, real_vector_name, pointer_2_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: real_vector_name
    real(k_real), dimension(:,:,:,:), pointer, intent(inout) :: pointer_2_data
    call this%AGDGetVectorDataPointerByNameTypeAndStatefulIndex(real_vector_name, real_vector, 3, pointer_2_data)
  end subroutine

  ! get pointer to avergae data
  subroutine AGDGetAvgRealVectorDataPointerByName(this, real_vector_name, pointer_2_avg_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: real_vector_name
    real(k_real), dimension(:), pointer, intent(inout) :: pointer_2_avg_data
    call this%AGDGetAvgVectorDataPointerByNameTypeAndStatefulIndex(real_vector_name, real_vector, 1, pointer_2_avg_data)
  end subroutine

  subroutine AGDGetAvgRealVectorDataPointerByNameOld(this, real_vector_name, pointer_2_avg_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: real_vector_name
    real(k_real), dimension(:), pointer, intent(inout) :: pointer_2_avg_data
    call this%AGDGetAvgVectorDataPointerByNameTypeAndStatefulIndex(real_vector_name, real_vector, 2, pointer_2_avg_data)
  end subroutine

  subroutine AGDGetAvgRealVectorDataPointerByNameOlder(this, real_vector_name, pointer_2_avg_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: real_vector_name
    real(k_real), dimension(:), pointer, intent(inout) :: pointer_2_avg_data
    call this%AGDGetAvgVectorDataPointerByNameTypeAndStatefulIndex(real_vector_name, real_vector, 3, pointer_2_avg_data)
  end subroutine

  !************************************************************!
  !*******************VECTOR5 DATA API ************************!
  !************************************************************!

  ! get pointer to data
  subroutine AGDGetVector5DataPointerByName(this, vector5_name, pointer_2_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: vector5_name
    real(k_real), dimension(:,:,:,:), pointer, intent(inout) :: pointer_2_data
    call this%AGDGetVectorDataPointerByNameTypeAndStatefulIndex(vector5_name, vector5, 1, pointer_2_data)
  end subroutine

  subroutine AGDGetVector5DataPointerByNameOld(this, vector5_name, pointer_2_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: vector5_name
    real(k_real), dimension(:,:,:,:), pointer, intent(inout) :: pointer_2_data
    call this%AGDGetVectorDataPointerByNameTypeAndStatefulIndex(vector5_name, vector5, 2, pointer_2_data)
  end subroutine

  subroutine AGDGetVector5DataPointerByNameOlder(this, vector5_name, pointer_2_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: vector5_name
    real(k_real), dimension(:,:,:,:), pointer, intent(inout) :: pointer_2_data
    call this%AGDGetVectorDataPointerByNameTypeAndStatefulIndex(vector5_name, vector5, 3, pointer_2_data)
  end subroutine

  ! get pointer to average data
  subroutine AGDGetAvgVector5DataPointerByName(this, vector5_name, pointer_2_avg_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: vector5_name
    real(k_real), dimension(:), pointer, intent(inout) :: pointer_2_avg_data
    call this%AGDGetAvgVectorDataPointerByNameTypeAndStatefulIndex(vector5_name, vector5, 1, pointer_2_avg_data)
  end subroutine

  subroutine AGDGetAvgVector5DataPointerByNameOld(this, vector5_name, pointer_2_avg_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: vector5_name
    real(k_real), dimension(:), pointer, intent(inout) :: pointer_2_avg_data
    call this%AGDGetAvgVectorDataPointerByNameTypeAndStatefulIndex(vector5_name, vector5, 2, pointer_2_avg_data)
  end subroutine

  subroutine AGDGetAvgVector5DataPointerByNameOlder(this, vector5_name, pointer_2_avg_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: vector5_name
    real(k_real), dimension(:), pointer, intent(inout) :: pointer_2_avg_data
    call this%AGDGetAvgVectorDataPointerByNameTypeAndStatefulIndex(vector5_name, vector5, 3, pointer_2_avg_data)
  end subroutine

  !************************************************************!
  !*******************VECTOR6 DATA API ************************!
  !************************************************************!

  ! get pointer to data
  subroutine AGDGetVector6DataPointerByName(this, vector6_name, pointer_2_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: vector6_name
    real(k_real), dimension(:,:,:,:), pointer, intent(inout) :: pointer_2_data
    call this%AGDGetVectorDataPointerByNameTypeAndStatefulIndex(vector6_name, vector6, 1, pointer_2_data)
  end subroutine

  subroutine AGDGetVector6DataPointerByNameOld(this, vector6_name, pointer_2_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: vector6_name
    real(k_real), dimension(:,:,:,:), pointer, intent(inout) :: pointer_2_data
    call this%AGDGetVectorDataPointerByNameTypeAndStatefulIndex(vector6_name, vector6, 2, pointer_2_data)
  end subroutine

  subroutine AGDGetVector6DataPointerByNameOlder(this, vector6_name, pointer_2_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: vector6_name
    real(k_real), dimension(:,:,:,:), pointer, intent(inout) :: pointer_2_data
    call this%AGDGetVectorDataPointerByNameTypeAndStatefulIndex(vector6_name, vector6, 3, pointer_2_data)
  end subroutine

  ! get pointer to average data
  subroutine AGDGetAvgVector6DataPointerByName(this, vector6_name, pointer_2_avg_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: vector6_name
    real(k_real), dimension(:), pointer, intent(inout) :: pointer_2_avg_data
    call this%AGDGetAvgVectorDataPointerByNameTypeAndStatefulIndex(vector6_name, vector6, 1, pointer_2_avg_data)
  end subroutine

  subroutine AGDGetAvgVector6DataPointerByNameOld(this, vector6_name, pointer_2_avg_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: vector6_name
    real(k_real), dimension(:), pointer, intent(inout) :: pointer_2_avg_data
    call this%AGDGetAvgVectorDataPointerByNameTypeAndStatefulIndex(vector6_name, vector6, 2, pointer_2_avg_data)
  end subroutine

  subroutine AGDGetAvgVector6DataPointerByNameOlder(this, vector6_name, pointer_2_avg_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: vector6_name
    real(k_real), dimension(:), pointer, intent(inout) :: pointer_2_avg_data
    call this%AGDGetAvgVectorDataPointerByNameTypeAndStatefulIndex(vector6_name, vector6, 3, pointer_2_avg_data)
  end subroutine

  ! get pointer to variable
  subroutine AGDGetVector6PointerToVariableByName(this, vector6_name, pointer_2_variable)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: vector6_name
    class(griddata_vector6), pointer, intent(inout) :: pointer_2_variable
    class(griddata_vector), pointer :: temp_pointer_to_variable => null()

    nullify(temp_pointer_to_variable)

    call this%AGDGetGridDataVectorVariablePointerByName(vector6_name, vector6, temp_pointer_to_variable)

    select type(temp_pointer_to_variable)
    class is (griddata_vector6)
        pointer_2_variable => temp_pointer_to_variable
    class default
        error stop "AGDGetVector6PointerToVariableByName wrong class"
    end select

  end subroutine

  !************************************************************!
  !************GENERIC VECTOR DATA API ************************!
  !************************************************************!

  ! get pointer to data
  subroutine AGDGetGenericVectorDataPointerByName(this, generic_vector_name, pointer_2_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: generic_vector_name
    real(k_real), dimension(:,:,:,:), pointer, intent(inout) :: pointer_2_data
    call this%AGDGetVectorDataPointerByNameTypeAndStatefulIndex(generic_vector_name, generic_vector, 1, pointer_2_data)
  end subroutine

  subroutine AGDGetGenericVectorDataPointerByNameOld(this, generic_vector_name, pointer_2_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: generic_vector_name
    real(k_real), dimension(:,:,:,:), pointer, intent(inout) :: pointer_2_data
    call this%AGDGetVectorDataPointerByNameTypeAndStatefulIndex(generic_vector_name, generic_vector, 2, pointer_2_data)
  end subroutine

  subroutine AGDGetGenericVectorDataPointerByNameOlder(this, generic_vector_name, pointer_2_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: generic_vector_name
    real(k_real), dimension(:,:,:,:), pointer, intent(inout) :: pointer_2_data
    call this%AGDGetVectorDataPointerByNameTypeAndStatefulIndex(generic_vector_name, generic_vector, 3, pointer_2_data)
  end subroutine

  ! get pointer to average data
  subroutine AGDGetAvgGenericVectorDataPointerByName(this, generic_vector_name, pointer_2_avg_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: generic_vector_name
    real(k_real), dimension(:), pointer, intent(inout) :: pointer_2_avg_data
    call this%AGDGetAvgVectorDataPointerByNameTypeAndStatefulIndex(generic_vector_name, generic_vector, 1, pointer_2_avg_data)
  end subroutine

  subroutine AGDGetAvgGenericVectorDataPointerByNameOld(this, generic_vector_name, pointer_2_avg_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: generic_vector_name
    real(k_real), dimension(:), pointer, intent(inout) :: pointer_2_avg_data
    call this%AGDGetAvgVectorDataPointerByNameTypeAndStatefulIndex(generic_vector_name, generic_vector, 2, pointer_2_avg_data)
  end subroutine

  subroutine AGDGetAvgGenericVectorDataPointerByNameOlder(this, generic_vector_name, pointer_2_avg_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: generic_vector_name
    real(k_real), dimension(:), pointer, intent(inout) :: pointer_2_avg_data
    call this%AGDGetAvgVectorDataPointerByNameTypeAndStatefulIndex(generic_vector_name, generic_vector, 3, pointer_2_avg_data)
  end subroutine

  !************************************************************!
  !*********SLIPSSYTEM VECTOR SCALAR DATA API *****************!
  !************************************************************!

  ! get pointer to data
  subroutine AGDgetSSScalarDataPointerByName(this, ss_vector_name, pointer_2_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: ss_vector_name
    real(k_real), dimension(:,:,:,:), pointer, intent(inout) :: pointer_2_data
    call this%AGDGetVectorDataPointerByNameTypeAndStatefulIndex(ss_vector_name, ss_scalar, 1, pointer_2_data)
  end subroutine

  subroutine AGDgetSSScalarDataPointerByNameOld(this, ss_vector_name, pointer_2_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: ss_vector_name
    real(k_real), dimension(:,:,:,:), pointer, intent(inout) :: pointer_2_data
    call this%AGDGetVectorDataPointerByNameTypeAndStatefulIndex(ss_vector_name, ss_scalar, 2, pointer_2_data)
  end subroutine

  subroutine AGDgetSSScalarDataPointerByNameOlder(this, ss_vector_name, pointer_2_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: ss_vector_name
    real(k_real), dimension(:,:,:,:), pointer, intent(inout) :: pointer_2_data
    call this%AGDGetVectorDataPointerByNameTypeAndStatefulIndex(ss_vector_name, ss_scalar, 3, pointer_2_data)
  end subroutine

  ! get pointer to average data
  subroutine AGDgetAvgSSScalarDataPointerByName(this, ss_vector_name, pointer_2_avg_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: ss_vector_name
    real(k_real), dimension(:), pointer, intent(inout) :: pointer_2_avg_data
    call this%AGDGetAvgVectorDataPointerByNameTypeAndStatefulIndex(ss_vector_name, ss_scalar, 1, pointer_2_avg_data)
  end subroutine

  subroutine AGDgetAvgSSScalarDataPointerByNameOld(this, ss_vector_name, pointer_2_avg_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: ss_vector_name
    real(k_real), dimension(:), pointer, intent(inout) :: pointer_2_avg_data
    call this%AGDGetAvgVectorDataPointerByNameTypeAndStatefulIndex(ss_vector_name, ss_scalar, 2, pointer_2_avg_data)
  end subroutine

  subroutine AGDgetAvgSSScalarDataPointerByNameOlder(this, ss_vector_name, pointer_2_avg_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: ss_vector_name
    real(k_real), dimension(:), pointer, intent(inout) :: pointer_2_avg_data
    call this%AGDGetAvgVectorDataPointerByNameTypeAndStatefulIndex(ss_vector_name, ss_scalar, 3, pointer_2_avg_data)
  end subroutine

  !************************************************************!
  !******************TENSOR2 SCALAR DATA API ******************!
  !************************************************************!

  ! get pointer to data
  subroutine AGDGetTensor2DataPointerByName(this, tensor2_name, pointer_2_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: tensor2_name
    real(k_real), dimension(:,:,:,:,:), pointer, intent(inout) :: pointer_2_data
    call this%AGDGetMatrixDataPointerByNameTypeAndStatefulIndex(tensor2_name, tensor2, 1, pointer_2_data)
  end subroutine

  subroutine AGDGetTensor2DataPointerByNameOld(this, tensor2_name, pointer_2_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: tensor2_name
    real(k_real), dimension(:,:,:,:,:), pointer, intent(inout) :: pointer_2_data
    call this%AGDGetMatrixDataPointerByNameTypeAndStatefulIndex(tensor2_name, tensor2, 2, pointer_2_data)
  end subroutine

  subroutine AGDGetTensor2DataPointerByNameOlder(this, tensor2_name, pointer_2_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: tensor2_name
    real(k_real), dimension(:,:,:,:,:), pointer, intent(inout) :: pointer_2_data
    call this%AGDGetMatrixDataPointerByNameTypeAndStatefulIndex(tensor2_name, tensor2, 3, pointer_2_data)
  end subroutine

  subroutine AGDGetAvgTensor2DataPointerByName(this, tensor2_name, pointer_2_avg_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: tensor2_name
    real(k_real), dimension(:,:), pointer, intent(inout) :: pointer_2_avg_data
    call this%AGDGetAvgMatrixDataPointerByNameTypeAndStatefulIndex(tensor2_name, tensor2, 1, pointer_2_avg_data)
  end subroutine

  subroutine AGDGetAvgTensor2DataPointerByNameOld(this, tensor2_name, pointer_2_avg_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: tensor2_name
    real(k_real), dimension(:,:), pointer, intent(inout) :: pointer_2_avg_data
    call this%AGDGetAvgMatrixDataPointerByNameTypeAndStatefulIndex(tensor2_name, tensor2, 2, pointer_2_avg_data)
  end subroutine

  subroutine AGDGetAvgTensor2DataPointerByNameOlder(this, tensor2_name, pointer_2_avg_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: tensor2_name
    real(k_real), dimension(:,:), pointer, intent(inout) :: pointer_2_avg_data
    call this%AGDGetAvgMatrixDataPointerByNameTypeAndStatefulIndex(tensor2_name, tensor2, 3, pointer_2_avg_data)
  end subroutine

  ! get pointer to variable
  subroutine AGDGetTensor2PointerToVariableByName(this, tensor2_name, pointer_2_variable)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: tensor2_name
    class(griddata_tensor2), pointer, intent(inout) :: pointer_2_variable
    class(griddata_matrix), pointer :: temp_pointer_to_variable => null()

    nullify(temp_pointer_to_variable)

    call this%AGDGetGridDataMatrixVariablePointerByName(tensor2_name, tensor2, temp_pointer_to_variable)

    select type(temp_pointer_to_variable)
    class is (griddata_tensor2)
        pointer_2_variable => temp_pointer_to_variable
    class default
        error stop "AGDGetTensor2PointerToVariableByName wrong class"
    end select

  end subroutine

  !************************************************************!
  !*****************GENERIC MATRIX DATA API *******************!
  !************************************************************! 

  ! get pointer to data
  subroutine AGDGetGenericMatrixDataPointerByName(this, generic_matrix_name, pointer_2_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: generic_matrix_name
    real(k_real), dimension(:,:,:,:,:), pointer, intent(inout) :: pointer_2_data
    call this%AGDGetMatrixDataPointerByNameTypeAndStatefulIndex(generic_matrix_name, generic_matrix, 1, pointer_2_data)
  end subroutine

  subroutine AGDGetGenericMatrixDataPointerByNameOld(this, generic_matrix_name, pointer_2_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: generic_matrix_name
    real(k_real), dimension(:,:,:,:,:), pointer, intent(inout) :: pointer_2_data
    call this%AGDGetMatrixDataPointerByNameTypeAndStatefulIndex(generic_matrix_name, generic_matrix, 2, pointer_2_data)
  end subroutine

  subroutine AGDGetGenericMatrixDataPointerByNameOlder(this, generic_matrix_name, pointer_2_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: generic_matrix_name
    real(k_real), dimension(:,:,:,:,:), pointer, intent(inout) :: pointer_2_data
    call this%AGDGetMatrixDataPointerByNameTypeAndStatefulIndex(generic_matrix_name, generic_matrix, 3, pointer_2_data)
  end subroutine

  ! get pointer to average data
  subroutine AGDGetAvgGenericMatrixDataPointerByName(this, generic_matrix_name, pointer_2_avg_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: generic_matrix_name
    real(k_real), dimension(:,:), pointer, intent(inout) :: pointer_2_avg_data
    call this%AGDGetAvgMatrixDataPointerByNameTypeAndStatefulIndex(generic_matrix_name, generic_matrix, 1, pointer_2_avg_data)
  end subroutine

  subroutine AGDGetAvgGenericMatrixDataPointerByNameOld(this, generic_matrix_name, pointer_2_avg_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: generic_matrix_name
    real(k_real), dimension(:,:), pointer, intent(inout) :: pointer_2_avg_data
    call this%AGDGetAvgMatrixDataPointerByNameTypeAndStatefulIndex(generic_matrix_name, generic_matrix, 2, pointer_2_avg_data)
  end subroutine

  subroutine AGDGetAvgGenericMatrixDataPointerByNameOlder(this, generic_matrix_name, pointer_2_avg_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: generic_matrix_name
    real(k_real), dimension(:,:), pointer, intent(inout) :: pointer_2_avg_data
    call this%AGDGetAvgMatrixDataPointerByNameTypeAndStatefulIndex(generic_matrix_name, generic_matrix, 3, pointer_2_avg_data)
  end subroutine

  !************************************************************!
  !*****************MATRIX66 SCALAR DATA API ******************!
  !************************************************************!

  ! get pointer to data
  subroutine AGDGetMatrix66DataPointerByName(this, matrix66_name, pointer_2_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: matrix66_name
    real(k_real), dimension(:,:,:,:,:), pointer, intent(inout) :: pointer_2_data
    call this%AGDGetMatrixDataPointerByNameTypeAndStatefulIndex(matrix66_name, matrix66, 1, pointer_2_data)
  end subroutine

  subroutine AGDGetMatrix66DataPointerByNameOld(this, matrix66_name, pointer_2_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: matrix66_name
    real(k_real), dimension(:,:,:,:,:), pointer, intent(inout) :: pointer_2_data
    call this%AGDGetMatrixDataPointerByNameTypeAndStatefulIndex(matrix66_name, matrix66, 2, pointer_2_data)
  end subroutine

  subroutine AGDGetMatrix66DataPointerByNameOlder(this, matrix66_name, pointer_2_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: matrix66_name
    real(k_real), dimension(:,:,:,:,:), pointer, intent(inout) :: pointer_2_data
    call this%AGDGetMatrixDataPointerByNameTypeAndStatefulIndex(matrix66_name, matrix66, 3, pointer_2_data)
  end subroutine

  ! get pointer to average data
  subroutine AGDGetAvgMatrix66DataPointerByName(this, matrix66_name, pointer_2_avg_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: matrix66_name
    real(k_real), dimension(:,:), pointer, intent(inout) :: pointer_2_avg_data
    call this%AGDGetAvgMatrixDataPointerByNameTypeAndStatefulIndex(matrix66_name, matrix66, 1, pointer_2_avg_data)
  end subroutine

  subroutine AGDGetAvgMatrix66DataPointerByNameOld(this, matrix66_name, pointer_2_avg_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: matrix66_name
    real(k_real), dimension(:,:), pointer, intent(inout) :: pointer_2_avg_data
    call this%AGDGetAvgMatrixDataPointerByNameTypeAndStatefulIndex(matrix66_name, matrix66, 2, pointer_2_avg_data)
  end subroutine

  subroutine AGDGetAvgMatrix66DataPointerByNameOlder(this, matrix66_name, pointer_2_avg_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: matrix66_name
    real(k_real), dimension(:,:), pointer, intent(inout) :: pointer_2_avg_data
    call this%AGDGetAvgMatrixDataPointerByNameTypeAndStatefulIndex(matrix66_name, matrix66, 3, pointer_2_avg_data)
  end subroutine

  !*************************************************************!
  !*********SLIPSSYTEM GENERIC VECTOR DATA API *****************!
  !*************************************************************!

  ! get pointer to data
  subroutine AGDGetSSGenericVectorDataPointerByName(this, ss_generic_vector_name, pointer_2_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: ss_generic_vector_name
    real(k_real), dimension(:,:,:,:,:), pointer, intent(inout) :: pointer_2_data
    call this%AGDGetMatrixDataPointerByNameTypeAndStatefulIndex(ss_generic_vector_name, ss_generic_vector, 1, pointer_2_data)
  end subroutine

  subroutine AGDGetSSGenericVectorDataPointerByNameOld(this, ss_generic_vector_name, pointer_2_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: ss_generic_vector_name
    real(k_real), dimension(:,:,:,:,:), pointer, intent(inout) :: pointer_2_data
    call this%AGDGetMatrixDataPointerByNameTypeAndStatefulIndex(ss_generic_vector_name, ss_generic_vector, 2, pointer_2_data)
  end subroutine

  subroutine AGDGetSSGenericVectorDataPointerByNameOlder(this, ss_generic_vector_name, pointer_2_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: ss_generic_vector_name
    real(k_real), dimension(:,:,:,:,:), pointer, intent(inout) :: pointer_2_data
    call this%AGDGetMatrixDataPointerByNameTypeAndStatefulIndex(ss_generic_vector_name, ss_generic_vector, 3, pointer_2_data)
  end subroutine

  ! get pointer to average data
  subroutine AGDGetAvgSSGenericVectorDataPointerByName(this, ss_generic_vector_name, pointer_2_avg_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: ss_generic_vector_name
    real(k_real), dimension(:,:), pointer, intent(inout) :: pointer_2_avg_data
    call this%AGDGetAvgMatrixDataPointerByNameTypeAndStatefulIndex(ss_generic_vector_name, ss_generic_vector, 1, pointer_2_avg_data)
  end subroutine

  subroutine AGDGetAvgSSGenericVectorDataPointerByNameOld(this, ss_generic_vector_name, pointer_2_avg_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: ss_generic_vector_name
    real(k_real), dimension(:,:), pointer, intent(inout) :: pointer_2_avg_data
    call this%AGDGetAvgMatrixDataPointerByNameTypeAndStatefulIndex(ss_generic_vector_name, ss_generic_vector, 2, pointer_2_avg_data)
  end subroutine

  subroutine AGDGetAvgSSGenericVectorDataPointerByNameOlder(this, ss_generic_vector_name, pointer_2_avg_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: ss_generic_vector_name
    real(k_real), dimension(:,:), pointer, intent(inout) :: pointer_2_avg_data
    call this%AGDGetAvgMatrixDataPointerByNameTypeAndStatefulIndex(ss_generic_vector_name, ss_generic_vector, 3, pointer_2_avg_data)
  end subroutine

  !*************************************************************!
  !************SLIPSSYTEM VECTOR5 DATA API *********************!
  !*************************************************************!

  ! get pointer to data
  subroutine AGDGetSSVector5DataPointerByName(this, ss_vector5_name, pointer_2_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: ss_vector5_name
    real(k_real), dimension(:,:,:,:,:), pointer, intent(inout) :: pointer_2_data
    call this%AGDGetMatrixDataPointerByNameTypeAndStatefulIndex(ss_vector5_name, ss_vector5, 1, pointer_2_data)
  end subroutine

  subroutine AGDGetSSVector5DataPointerByNameOld(this, ss_vector5_name, pointer_2_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: ss_vector5_name
    real(k_real), dimension(:,:,:,:,:), pointer, intent(inout) :: pointer_2_data
    call this%AGDGetMatrixDataPointerByNameTypeAndStatefulIndex(ss_vector5_name, ss_vector5, 2, pointer_2_data)
  end subroutine

  subroutine AGDGetSSVector5DataPointerByNameOlder(this, ss_vector5_name, pointer_2_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: ss_vector5_name
    real(k_real), dimension(:,:,:,:,:), pointer, intent(inout) :: pointer_2_data
    call this%AGDGetMatrixDataPointerByNameTypeAndStatefulIndex(ss_vector5_name, ss_vector5, 3, pointer_2_data)
  end subroutine

  ! get pointer to average data
  subroutine AGDGetAvgSSVector5DataPointerByName(this, ss_vector5_name, pointer_2_avg_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: ss_vector5_name
    real(k_real), dimension(:,:), pointer, intent(inout) :: pointer_2_avg_data
    call this%AGDGetAvgMatrixDataPointerByNameTypeAndStatefulIndex(ss_vector5_name, ss_vector5, 1, pointer_2_avg_data)
  end subroutine

  subroutine AGDGetAvgSSVector5DataPointerByNameOld(this, ss_vector5_name, pointer_2_avg_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: ss_vector5_name
    real(k_real), dimension(:,:), pointer, intent(inout) :: pointer_2_avg_data
    call this%AGDGetAvgMatrixDataPointerByNameTypeAndStatefulIndex(ss_vector5_name, ss_vector5, 2, pointer_2_avg_data)
  end subroutine

  subroutine AGDGetAvgSSVector5DataPointerByNameOlder(this, ss_vector5_name, pointer_2_avg_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: ss_vector5_name
    real(k_real), dimension(:,:), pointer, intent(inout) :: pointer_2_avg_data
    call this%AGDGetAvgMatrixDataPointerByNameTypeAndStatefulIndex(ss_vector5_name, ss_vector5, 3, pointer_2_avg_data)
  end subroutine

  !*************************************************************!
  !************SLIPSSYTEM VECTOR6 DATA API *********************!
  !*************************************************************!

  ! get pointer to data
  subroutine AGDGetSSVector6DataPointerByName(this, ss_vector6_name, pointer_2_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: ss_vector6_name
    real(k_real), dimension(:,:,:,:,:), pointer, intent(inout) :: pointer_2_data
    call this%AGDGetMatrixDataPointerByNameTypeAndStatefulIndex(ss_vector6_name, ss_vector6, 1, pointer_2_data)
  end subroutine

  subroutine AGDGetSSVector6DataPointerByNameOld(this, ss_vector6_name, pointer_2_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: ss_vector6_name
    real(k_real), dimension(:,:,:,:,:), pointer, intent(inout) :: pointer_2_data
    call this%AGDGetMatrixDataPointerByNameTypeAndStatefulIndex(ss_vector6_name, ss_vector6, 2, pointer_2_data)
  end subroutine

  subroutine AGDGetSSVector6DataPointerByNameOlder(this, ss_vector6_name, pointer_2_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: ss_vector6_name
    real(k_real), dimension(:,:,:,:,:), pointer, intent(inout) :: pointer_2_data
    call this%AGDGetMatrixDataPointerByNameTypeAndStatefulIndex(ss_vector6_name, ss_vector6, 3, pointer_2_data)
  end subroutine

  ! get pointer to average data
  subroutine AGDGetAvgSSVector6DataPointerByName(this, ss_vector6_name, pointer_2_avg_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: ss_vector6_name
    real(k_real), dimension(:,:), pointer, intent(inout) :: pointer_2_avg_data
    call this%AGDGetAvgMatrixDataPointerByNameTypeAndStatefulIndex(ss_vector6_name, ss_vector6, 1, pointer_2_avg_data)
  end subroutine

  subroutine AGDGetAvgSSVector6DataPointerByNameOld(this, ss_vector6_name, pointer_2_avg_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: ss_vector6_name
    real(k_real), dimension(:,:), pointer, intent(inout) :: pointer_2_avg_data
    call this%AGDGetAvgMatrixDataPointerByNameTypeAndStatefulIndex(ss_vector6_name, ss_vector6, 2, pointer_2_avg_data)
  end subroutine

  subroutine AGDGetAvgSSVector6DataPointerByNameOlder(this, ss_vector6_name, pointer_2_avg_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: ss_vector6_name
    real(k_real), dimension(:,:), pointer, intent(inout) :: pointer_2_avg_data
    call this%AGDGetAvgMatrixDataPointerByNameTypeAndStatefulIndex(ss_vector6_name, ss_vector6, 3, pointer_2_avg_data)
  end subroutine

  !*************************************************************!
  !************SLIPSSYTEM MATRIX DATA API *********************!
  !*************************************************************!

  ! get pointer to data
  subroutine AGDGetSSMatrixDataPointerByName(this, ss_vector6_name, pointer_2_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: ss_vector6_name
    real(k_real), dimension(:,:,:,:,:,:), pointer, intent(inout) :: pointer_2_data
    call this%AGDGetSSMatrixDataPointerByNameTypeAndStatefulIndex(ss_vector6_name, ss_vector6, 1, pointer_2_data)
  end subroutine

  subroutine AGDGetSSMatrixDataPointerByNameOld(this, ss_vector6_name, pointer_2_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: ss_vector6_name
    real(k_real), dimension(:,:,:,:,:,:), pointer, intent(inout) :: pointer_2_data
    call this%AGDGetSSMatrixDataPointerByNameTypeAndStatefulIndex(ss_vector6_name, ss_vector6, 2, pointer_2_data)
  end subroutine

  subroutine AGDGetSSMatrixDataPointerByNameOlder(this, ss_vector6_name, pointer_2_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: ss_vector6_name
    real(k_real), dimension(:,:,:,:,:,:), pointer, intent(inout) :: pointer_2_data
    call this%AGDGetSSMatrixDataPointerByNameTypeAndStatefulIndex(ss_vector6_name, ss_vector6, 3, pointer_2_data)
  end subroutine

  ! get pointer to average data
  subroutine AGDGetAvgSSMatrixDataPointerByName(this, ss_vector6_name, pointer_2_avg_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: ss_vector6_name
    real(k_real), dimension(:,:,:), pointer, intent(inout) :: pointer_2_avg_data
    call this%AGDGetAvgSSMatrixDataPointerByNameTypeAndStatefulIndex(ss_vector6_name, ss_vector6, 1, pointer_2_avg_data)
  end subroutine

  subroutine AGDGetAvgSSMatrixDataPointerByNameOld(this, ss_vector6_name, pointer_2_avg_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: ss_vector6_name
    real(k_real), dimension(:,:,:), pointer, intent(inout) :: pointer_2_avg_data
    call this%AGDGetAvgSSMatrixDataPointerByNameTypeAndStatefulIndex(ss_vector6_name, ss_vector6, 2, pointer_2_avg_data)
  end subroutine

  subroutine AGDGetAvgSSMatrixDataPointerByNameOlder(this, ss_vector6_name, pointer_2_avg_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: ss_vector6_name
    real(k_real), dimension(:,:,:), pointer, intent(inout) :: pointer_2_avg_data
    call this%AGDGetAvgSSMatrixDataPointerByNameTypeAndStatefulIndex(ss_vector6_name, ss_vector6, 3, pointer_2_avg_data)
  end subroutine

  subroutine addVar(this, var_name, var_type, additional_var_dimensions, stateful_level)
    implicit none
    class(all_grid_data), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    integer(kind(grid_data_var_type_enum)), intent(in) :: var_type
    integer, dimension(:), optional, intent(in) :: additional_var_dimensions
    integer, intent(in), optional :: stateful_level

    call this%all_grid_vars_list%addVar(var_name, var_type, additional_var_dimensions, stateful_level)
  end subroutine

  subroutine addGridData(this, new_element)
    class(all_grid_data_multi_phase), intent(inout) :: this
    class(all_grid_data), pointer, intent(inout) :: new_element

    call this%extend()
    this%all_pt(this%n_elements)%pt => new_element
    nullify(new_element)
  end subroutine


  subroutine getGridDataPtr(this, idx, ptr)
    class(all_grid_data_multi_phase), intent(inout) :: this
    integer, intent(in) :: idx
    class(all_grid_data), pointer, intent(out) :: ptr

    call this%checkElementExist(idx)

    associate(elem => this%all_pt(idx)%pt)
    select type(elem)
      class is (all_grid_data)
        ptr => elem
      class default
        error stop "getGridDataPtr wrong class"
    end select
    end associate

  end subroutine

end module all_grid_data_mod
