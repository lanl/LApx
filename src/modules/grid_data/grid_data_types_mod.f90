module grid_data_types
  use kinds
  use grid_data_var_type_mod
#include "macro_debug.fpp"

  implicit none

  type var_dimension_type
  integer, dimension(3) :: nx_ny_nz=(/0,0,0/), &
                           nx_ny_nz_rank =(/0,0,0/), &
                           xyz_start_rank=(/0,0,0/), &
                           xyz_end_rank=(/0,0,0/), &
                           xyz_offset_rank=(/0,0,0/)
  contains
    procedure :: init => initVarDimensionType
  end type

  type, abstract :: abstract_grid_data_type
    character(len=:), allocatable :: name
    class(var_dimension_type), pointer :: dimension_obj => null()
    integer :: nx, ny, nz, npoints, nz_total, z_offset_rank
    integer, dimension(3) :: nx_ny_nz, &
                             nx_ny_nz_rank, &
                             xyz_start_rank, &
                             xyz_end_rank, &
                             xyz_offset_rank
    logical :: initialized = .FALSE.
    logical :: stateful = .FALSE.
    integer :: stateful_level = 1 !-> the number of time step we want to keep in memory
    logical :: compute_field_average=.FALSE. !-> if true computes the field average
    integer(kind(grid_data_var_type_enum)) :: data_type
    logical :: dummy_logical
  contains
    procedure, public :: checkDimensionBase
    procedure, public :: checkGridDataDimension => checkDimensionBase
    procedure, public :: setStatefulLevel
    procedure, private :: errorOnGetOldOrOlder
    procedure, public  :: writePointToFile
    procedure, public  :: writeDataPointToFile => writeDataPointToFileBase
    procedure, private :: setVarDimsFromVarDimsType
  end type

  type, extends(abstract_grid_data_type) :: griddata_scalar
      real(k_real), dimension(:,:,:,:), pointer :: data  => null()
      real(k_real), dimension(:), pointer :: avg_data  => null()
    contains
      procedure, public :: init => initScalar
      procedure, public  :: getDataPointerByStatefulLevelIndex => getDataPointerByStatefulLevelIndexScalar
      procedure, public  :: getAvgDataPointerByStatefulLevelIndex => getAvgDataPointerByStatefulLevelIndexScalar
      procedure, public  :: computeAverage => computeAverageScalar
      procedure, public  :: updateHistory => updateHistoryScalar
      procedure, public  :: resetHistory => resetHistoryScalar
      procedure, public  :: setToValue => setScalarToValue
      procedure, public  :: setAvgToValue => setScalarAvgToValue
      procedure, public  :: dumpForRestart => dumpForRestartGridScalarStateful
      procedure, public  :: reloadFromDump => reloadFromDumpGridScalarStateful
      procedure, public  :: writeDataPointToFile => writeDataPointToFileScalar
      procedure, public  :: writeXDMFAttribute => writeXDMFAttributeScalarStateful
  end type

  type, extends(abstract_grid_data_type) :: griddata_scalar_integer
      integer(k_int), dimension(:,:,:,:), pointer :: data  => null()
    contains
      procedure, public :: init => initScalarInteger
      procedure, public  :: getDataPointerByStatefulLevelIndex => getDataPointerByStatefulLevelIndexScalarInteger
      procedure, public  :: updateHistory => updateHistoryScalarInteger
      procedure, public  :: resetHistory => resetHistoryScalarInteger
      procedure, public  :: setToValue => setScalarIntegerToValue
      procedure, public  :: dumpForRestart => dumpForRestartGridScalarIntegerStateful
      procedure, public  :: reloadFromDump => reloadFromDumpGridScalarIntegerStateful
      procedure, public  :: writeDataPointToFile => writeDataPointToFileIntegerScalar
      procedure, public  :: writeXDMFAttribute => writeXDMFAttributeScalarIntegerStateful
  end type

  type, abstract, extends(abstract_grid_data_type) :: griddata_vector
      integer ::  n_components
      real(k_real), dimension(:,:,:,:,:), pointer :: data  => null()
      real(k_real), dimension(:,:), pointer :: avg_data  => null()
    contains
      procedure, public  :: getDataPointerByStatefulLevelIndex => getDataPointerByStatefulLevelIndexVector
      procedure, public  :: getAvgDataPointerByStatefulLevelIndex => getAvgDataPointerByStatefulLevelIndexVector
      procedure, public  :: computeAverage => computeAverageVector
      procedure, public  :: computeAverageMasked => computeAverageVectorMasked
      procedure, public  :: updateHistory => updateHistoryVector
      procedure, public  :: resetHistory => resetHistoryVector
      procedure, public  :: setToValue => setVectorToValue
      procedure, public  :: setAvgToValue => setVectorAvgToValue
      procedure, public  :: dumpForRestart => dumpForRestartGridVectorStateful
      procedure, public  :: reloadFromDump => reloadFromDumpGridVectorStateful
      procedure, public  :: writeDataPointToFile => writeDataPointToFileVector
      procedure, public  :: writeXDMFAttribute => writeXDMFAttributeGridVectorStateful
  end type

  type, extends(griddata_vector) :: griddata_real_vector
    contains
      procedure, public :: init => initRealVector
  end type

  type, extends(griddata_vector) :: griddata_vector5
    contains
      procedure, public :: init => initVector5
  end type

  type, extends(griddata_vector) :: griddata_vector6
    contains
      procedure, public :: init => initVector6
  end type

  type, extends(griddata_vector) :: griddata_ssscalar
    integer ::  n_ss
    contains
      procedure, public :: init => initSSScalar
  end type

  type, extends(griddata_vector) :: griddata_generic_vector

    contains
      procedure, public :: init => initGenericVector
      procedure :: checkGridDataDimension => checkGridDataDimensionGenericVector
  end type

  type, abstract, extends(abstract_grid_data_type) :: griddata_matrix
    real(k_real), dimension(:,:,:,:,:,:), pointer :: data  => null()
    real(k_real), dimension(:,:,:), pointer :: avg_data  => null()
    integer :: n_comp_i, n_comp_j
  contains
    procedure, public  :: getDataPointerByStatefulLevelIndex => getDataPointerByStatefulLevelIndexMatrix
    procedure, public  :: getAvgDataPointerByStatefulLevelIndex => getAvgDataPointerByStatefulLevelIndexMatrix
    procedure, public  :: computeAverage => computeAverageMatrix
    procedure, public  :: computeAverageMasked => computeAverageMatrixMasked
    procedure, public  :: updateHistory => updateHistoryMatrix
    procedure, public  :: resetHistory => resetHistoryMatrix
    procedure, public  :: setToValue => setMatrixToValue
    procedure, public  :: setAvgToValue => setMatrixAvgToValue
    procedure, public  :: dumpForRestart => dumpForRestartGridMatrixStateful
    procedure, public  :: reloadFromDump => reloadFromDumpGridMatrixStateful
    procedure, public  :: writeDataPointToFile => writeDataPointToFileMatrix
    procedure, public  :: writeXDMFAttribute => writeXDMFAttributeGridMatrixStateful

  end type

  type, extends(griddata_matrix) :: griddata_tensor2
    contains
      procedure, public :: init => initTensor2
  end type

  type, extends(griddata_matrix) :: griddata_generic_matrix
    contains
      procedure, public :: init => initGenericMatrix
  end type

  type, extends(griddata_matrix) :: griddata_matrix66
    contains
      procedure, public :: init => initMatrix66
  end type

  type, extends(griddata_matrix) :: griddata_ss_generic_vector
    integer ::  n_components
    integer ::  n_ss
    contains
    procedure, public :: init => initSSGenericVector
  end type

  type, extends(griddata_matrix) :: griddata_ss_vector5
    integer ::  n_ss
    contains
    procedure, public :: init => initSSVector5
  end type

  type, extends(griddata_matrix) :: griddata_ss_vector6
    integer ::  n_ss
    contains
    procedure, public :: init => initSSVector6
  end type

  type, extends(abstract_grid_data_type) :: griddata_ss_generic_matrix
    integer ::  n_components
    integer ::  m_components
    integer ::  n_ss
    real(k_real), dimension(:,:,:,:,:,:,:), pointer :: data  => null()
    real(k_real), dimension(:,:,:,:), pointer :: avg_data  => null()
  contains
    procedure, public  :: init => initSSMatrix
    procedure, public  :: getDataPointerByStatefulLevelIndex => getDataPointerByStatefulLevelIndexSSMatrix
    procedure, public  :: getAvgDataPointerByStatefulLevelIndex => getAvgDataPointerByStatefulLevelIndexSSMatrix
    procedure, public  :: computeAverage => computeAverageSSMatrix
    procedure, public  :: updateHistory => updateHistorySSMatrix
    procedure, public  :: resetHistory => resetHistorySSMatrix
    procedure, public  :: setToValue => setSSMatrixToValue
    procedure, public  :: setAvgToValue => setSSMatrixAvgToValue
    procedure, public  :: dumpForRestart => dumpForRestartGridSSMatrixStateful
    procedure, public  :: reloadFromDump => reloadFromDumpGridSSMatrixStateful
    procedure, public  :: writeDataPointToFile => writeDataPointToFileSSMatrix
    procedure, public  :: writeXDMFAttribute => writeXDMFAttributeGridSSMatrix
  end type

interface
  module subroutine checkDimensionBase(this)
    implicit none
    class(abstract_grid_data_type), intent(in) :: this
  end subroutine

end interface

contains



  subroutine initVarDimensionType(this, nx_ny_nz, nx_ny_nz_rank, &
                                  xyz_start_rank, xyz_end_rank, xyz_offset_rank)
    implicit none 
    class(var_dimension_type), intent(inout) :: this
    integer, dimension(3), intent(in) :: nx_ny_nz, &
                            nx_ny_nz_rank, &
                            xyz_start_rank, &
                            xyz_end_rank, &
                            xyz_offset_rank

    this%nx_ny_nz = nx_ny_nz
    this%nx_ny_nz_rank = nx_ny_nz_rank
    this%xyz_start_rank = xyz_start_rank
    this%xyz_end_rank = xyz_end_rank
    this%xyz_offset_rank = xyz_offset_rank
  end subroutine

  subroutine setVarDimsFromVarDimsType(this, var_dims)
    class(abstract_grid_data_type), intent(inout) :: this 
    class(var_dimension_type), intent(in), pointer :: var_dims

    if (.not.(associated(var_dims))) error stop "setVarDimsFromVarDimsType: var_dims not associated"

    this%dimension_obj => var_dims
    this%nx_ny_nz = var_dims%nx_ny_nz
    this%nx_ny_nz_rank = var_dims%nx_ny_nz_rank
    this%xyz_start_rank = var_dims%xyz_start_rank
    this%xyz_end_rank = var_dims%xyz_end_rank
    this%xyz_offset_rank = var_dims%xyz_offset_rank 

    this%nx = this%nx_ny_nz_rank(1)
    this%ny = this%nx_ny_nz_rank(2)
    this%nz = this%nx_ny_nz_rank(3)
    this%npoints = product(this%nx_ny_nz_rank)
    this%nz_total = this%nx_ny_nz(3)
    this%z_offset_rank = this%xyz_offset_rank(3)
  end subroutine

  subroutine writePointToFile(this, ix,iy,iz, file_id)
    implicit none
    class(abstract_grid_data_type), intent(in) :: this
    integer, intent(in) :: ix, iy, iz, file_id
    integer :: i
    write(file_id,*) "variable name: ", trim(adjustl(this%name))
    write(file_id,*) "stateful level: ", this%stateful_level
    do i=1,this%stateful_level
      if (i==1) write(file_id,*) "current"
      if (i==2) write(file_id,*) "old"
      if (i==3) write(file_id,*) "older"
      call this%writeDataPointToFile(ix,iy,iz, i, file_id)
    enddo
    write(file_id,*) ""
  end subroutine

  subroutine writeDataPointToFileBase(this, ix,iy,iz, s_level, file_id)
    implicit none
    class(abstract_grid_data_type), intent(in) :: this
    integer, intent(in) :: ix, iy, iz, s_level, file_id
    __DECL_CLASS_UNUSED_THIS__
    __DECL_UNUSED_INTEGER__

    __SUPPRESS_CLASS_UNUSED_THIS__
    __SUPPRESS_UNUSED_INTEGER__(ix)
    __SUPPRESS_UNUSED_INTEGER__(iy)
    __SUPPRESS_UNUSED_INTEGER__(iz)
    __SUPPRESS_UNUSED_INTEGER__(s_level)
    __SUPPRESS_UNUSED_INTEGER__(file_id)

    error stop "if you end up here you didn't override writeDataPointToFile in your subclass"
  end subroutine

  subroutine errorOnGetOldOrOlder(this, requested_stateful_level)
    implicit none
    class(abstract_grid_data_type), intent(in) :: this
    integer, intent(in) :: requested_stateful_level

    if (.not.(this%stateful)) then
      write(*,*) "the variable named ", this%name, " is not stateful, you can't request older values"
      stop
    else
      if (requested_stateful_level > this%stateful_level) then
        write(*,*) "the variable named ", this%name, " only keeps track of ", this%stateful_level, " timestep, you are asking for ", requested_stateful_level
        stop
      end if
    end if
  end subroutine

  subroutine setStatefulLevel(this, stateful_level)
    implicit none
    class(abstract_grid_data_type), intent(inout) :: this
    integer, intent(in) :: stateful_level
    this%stateful = .TRUE.
    if (stateful_level > 3) then
      write(*,*) "the maximum number of timestep allowed in memory is 3, you provided  ", stateful_level
      stop
    else if (stateful_level < 1) then
      write(*,*) "stateful_level must be >= 1, you provided  ", stateful_level
      stop
    end if
    this%stateful_level = stateful_level
  end subroutine

  subroutine initScalar(this, name, dimension_obj, stateful_level, compute_field_average)
    implicit none
    class(griddata_scalar), intent(inout) :: this
    class(var_dimension_type), pointer :: dimension_obj 
    character(len=:), allocatable :: name
    integer, intent(in), optional :: stateful_level
    logical, intent(in), optional :: compute_field_average
    if (present(stateful_level)) call this%setStatefulLevel(stateful_level)
    if (present(compute_field_average)) this%compute_field_average = compute_field_average
    this%name = name
    call this%setVarDimsFromVarDimsType(dimension_obj)
    allocate(this%data(this%nx,this%ny,this%nz, this%stateful_level))
    allocate(this%avg_data(this%stateful_level))

    call this%setToValue(0._k_real, .true.)
    call this%setAvgToValue(0._k_real, .true.)

    this%data_type = scalar
    this%initialized = .TRUE.
  end subroutine

  subroutine dumpForRestartGridScalarStateful(this, file_name, absolute_path, only_current_value)
    use hdf5, only : HID_T
    use read_write_parallel_h5_mod, only : readOrWriteDatasetStateful, openexistingh5readwrite, closeh5file
    implicit none
    class(griddata_scalar), intent(in) :: this
    character(len=*), intent(in) :: file_name, absolute_path
    INTEGER(HID_T) :: file_id, facc_plist
    logical, intent(in) :: only_current_value

    call openExistingH5ReadWrite(file_name, file_id, facc_plist) 
    if (.not. only_current_value) then
    call readOrWriteDatasetStateful(this%data, absolute_path//this%name, this%nx_ny_nz, this%xyz_offset_rank, file_id, "write")
    else 
    call readOrWriteDatasetStateful(this%data(:,:,:,1:1), absolute_path//this%name, this%nx_ny_nz, this%xyz_offset_rank, file_id, "write")
    endif
    call CloseH5File(file_id, facc_plist)
  end subroutine

  subroutine writeXDMFAttributeScalarStateful(this, xmf_writer, hdf5_fname, field_full_path, write_only_current, suffix)
    use xmf_writer_mod, only : xmf_writer_type
    implicit none
    class(griddata_scalar), intent(in) :: this
    type(xmf_writer_type), intent(inout) :: xmf_writer
    CHARACTER(len=*), intent(in) :: hdf5_fname, field_full_path, suffix
    logical, intent(in) :: write_only_current
    integer :: n_stateful_level 
    __DECL_UNUSED_CHARACTER__

    __SUPPRESS_UNUSED_STRING__(hdf5_fname)
    if (write_only_current) then
      n_stateful_level = 1
    else 
      n_stateful_level = this%stateful_level
    endif

    call xmf_writer%write3Dscalar( hdf5_fname, field_full_path, trim(adjustl(this%name)), n_stateful_level, suffix)
  end subroutine

  subroutine reloadFromDumpGridScalarStateful(this, file_id)
    use hdf5, only : HID_T
    use read_write_parallel_h5_mod, only : readOrWriteDatasetStateful
    implicit none
    class(griddata_scalar), intent(inout) :: this
    INTEGER(HID_T), intent(in) :: file_id

    call readOrWriteDatasetStateful(this%data, this%name, this%nx_ny_nz, this%xyz_offset_rank, file_id, "read")

  end subroutine

  subroutine writeDataPointToFileScalar(this, ix,iy,iz, s_level, file_id)
    implicit none
    class(griddata_scalar), intent(in) :: this
    integer, intent(in) :: ix, iy, iz, s_level, file_id
      write(file_id,*) "dimensions: ", 1
      write(file_id,*) this%data(ix,iy,iz, s_level)
  end subroutine

  subroutine setScalarToValue(this, val, setAll)
    implicit none
    class(griddata_scalar), intent(inout) :: this
    real(k_real), intent(in) :: val
    logical, optional, intent(in) :: setAll
    integer :: ix,iy,iz,istateful, max_index

    max_index = 1
    if (present(setAll)) then
      if (setAll) max_index=this%stateful_level
    endif

    do istateful=1,max_index
    do iz=1,this%nz
    do iy=1,this%ny
    do ix=1,this%nx
      this%data(ix,iy,iz,istateful) = val
    enddo
    enddo
    enddo
    enddo

  end subroutine

  subroutine setScalarAvgToValue(this, val, setAll)
    implicit none
    class(griddata_scalar), intent(inout) :: this
    real(k_real), intent(in) :: val
    logical, optional, intent(in) :: setAll
    integer :: istateful, max_index

    max_index = 1
    if (present(setAll)) then
      if (setAll) max_index=this%stateful_level
    endif

    do istateful=1,max_index
      this%avg_data(istateful) = val
    enddo

  end subroutine

  subroutine getDataPointerByStatefulLevelIndexScalar(this, stateful_idx, pt_2_data)
    implicit none
    class(griddata_scalar), intent(in) :: this
    integer, intent(in) :: stateful_idx
    real(k_real), dimension(:,:,:), pointer, intent(inout) :: pt_2_data
    call this%errorOnGetOldOrOlder(stateful_idx)
    pt_2_data => this%data(:,:,:,stateful_idx)
  end subroutine

  subroutine getAvgDataPointerByStatefulLevelIndexScalar(this, stateful_idx, pt_2_avg_data)
    implicit none
    class(griddata_scalar), intent(in) :: this
    integer, intent(in) :: stateful_idx
    real(k_real), pointer, intent(inout) :: pt_2_avg_data
    call this%errorOnGetOldOrOlder(stateful_idx)
    pt_2_avg_data => this%avg_data(stateful_idx)
  end subroutine

  subroutine computeAverageScalar(this)
    use mpi_useful_routines_mod, only : MPIAverageGridScalar
    implicit none
    class(griddata_scalar), intent(inout) :: this

    call MPIAverageGridScalar(this%data(:,:,:,1), this%avg_data(1))
  end subroutine

  subroutine updateHistoryScalar(this)
    implicit none
    class(griddata_scalar), intent(inout) :: this
    integer :: i
    if (this%stateful) then
      if (this%stateful_level > 1) then
        do i= this%stateful_level,2,-1
          this%data(:,:,:,i) = this%data(:,:,:,i-1)
          this%avg_data(i) = this%avg_data(i-1)
        end do
      end if
    end if
  end subroutine

  subroutine resetHistoryScalar(this)
    implicit none
    class(griddata_scalar), intent(inout) :: this
    if (this%stateful) then
      if (this%stateful_level > 1) then
        this%data(:,:,:,1) = this%data(:,:,:,2)
        this%avg_data(1) = this%avg_data(2)
      end if
    end if
  end subroutine

  subroutine dumpForRestartGridScalarIntegerStateful(this, file_name, absolute_path, only_current_value)
    use hdf5, only : HID_T
    use read_write_parallel_h5_mod, only : readOrWriteDatasetStateful, openexistingh5readwrite, closeh5file
    implicit none
    class(griddata_scalar_integer), intent(in) :: this
    character(len=*), intent(in) :: file_name, absolute_path
    INTEGER(HID_T) :: file_id, facc_plist
    logical, intent(in) :: only_current_value

    call openExistingH5ReadWrite(file_name, file_id, facc_plist) 
    if (.not. only_current_value) then
      call readOrWriteDatasetStateful(this%data, absolute_path//this%name, this%nx_ny_nz, this%xyz_offset_rank, file_id, "write")
      else 
      call readOrWriteDatasetStateful(this%data(:,:,:,1:1), absolute_path//this%name, this%nx_ny_nz, this%xyz_offset_rank, file_id, "write")
    endif
    call CloseH5File(file_id, facc_plist)

  end subroutine

  subroutine writeXDMFAttributeScalarIntegerStateful(this, xmf_writer, hdf5_fname, field_full_path, write_only_current, suffix)
    use hdf5, only : HID_T
    use read_write_parallel_h5_mod, only : readOrWriteDatasetStateful
    use xmf_writer_mod, only : xmf_writer_type
    implicit none
    class(griddata_scalar_integer), intent(in) :: this
    type(xmf_writer_type), intent(inout) :: xmf_writer
    CHARACTER(len=*), intent(in) :: hdf5_fname, field_full_path, suffix
    logical, intent(in) :: write_only_current
    integer :: n_stateful_level 

    if (write_only_current) then
      n_stateful_level = 1
    else 
      n_stateful_level = this%stateful_level
    endif

    call xmf_writer%write3Dscalar( hdf5_fname, field_full_path, trim(adjustl(this%name)), n_stateful_level, suffix)

  end subroutine

  subroutine reloadFromDumpGridScalarIntegerStateful(this, file_id)
    use hdf5, only : HID_T
    use read_write_parallel_h5_mod, only : readOrWriteDatasetStateful
    implicit none
    class(griddata_scalar_integer), intent(inout) :: this
    INTEGER(HID_T), intent(in) :: file_id

    call readOrWriteDatasetStateful(this%data, this%name, this%nx_ny_nz, this%xyz_offset_rank, file_id, "read")

  end subroutine

  subroutine initScalarInteger(this, name, dimension_obj, stateful_level)
    implicit none
    class(griddata_scalar_integer), intent(inout) :: this
    class(var_dimension_type), intent(in), pointer :: dimension_obj

    character(len=:), allocatable :: name
    integer, intent(in), optional :: stateful_level
    this%data_type = scalar_integer

    if (present(stateful_level)) call this%setStatefulLevel(stateful_level)
    this%name = name
    call this%setVarDimsFromVarDimsType(dimension_obj)
    allocate(this%data(this%nx,this%ny,this%nz, this%stateful_level))

    call this%setToValue(0_k_int, .true.)

    this%initialized = .TRUE.
  end subroutine

  subroutine writeDataPointToFileIntegerScalar(this, ix,iy,iz, s_level, file_id)
    implicit none
    class(griddata_scalar_integer), intent(in) :: this
    integer, intent(in) :: ix, iy, iz, s_level, file_id
    write(file_id,*) "dimensions: ", 1
    write(file_id,*) this%data(ix,iy,iz,s_level)

  end subroutine

  subroutine setScalarIntegerToValue(this, val, setAll)
    implicit none
    class(griddata_scalar_integer), intent(inout) :: this
    integer(k_int), intent(in) :: val
    logical, optional, intent(in) :: setAll
    integer :: ix,iy,iz,istateful, max_index

    max_index = 1
    if (present(setAll)) then
      if (setAll) max_index=this%stateful_level
    endif

    do istateful=1,max_index
    do iz=1,this%nz
    do iy=1,this%ny
    do ix=1,this%nx
      this%data(ix,iy,iz,istateful) = val
    enddo
    enddo
    enddo
    enddo

  end subroutine

  subroutine getDataPointerByStatefulLevelIndexScalarInteger(this, stateful_idx, pt_2_data)
    implicit none
    class(griddata_scalar_integer), intent(in) :: this
    integer, intent(in) :: stateful_idx
    integer(k_int), dimension(:,:,:), pointer, intent(inout) :: pt_2_data
    call this%errorOnGetOldOrOlder(stateful_idx)
    pt_2_data => this%data(:,:,:,stateful_idx)
  end subroutine

  subroutine updateHistoryScalarInteger(this)
    implicit none
    class(griddata_scalar_integer), intent(inout) :: this
    integer :: i
    if (this%stateful) then
      if (this%stateful_level > 1) then
        do i= this%stateful_level,2,-1
          this%data(:,:,:,i) = this%data(:,:,:,i-1)
        end do
      end if
    end if
  end subroutine

  subroutine resetHistoryScalarInteger(this)
    implicit none
    class(griddata_scalar_integer), intent(inout) :: this
    if (this%stateful) then
      if (this%stateful_level > 1) then
        this%data(:,:,:,1) = this%data(:,:,:,2)
      end if
    end if
  end subroutine

  subroutine writeDataPointToFileVector(this, ix,iy,iz, s_level, file_id)
    implicit none
    class(griddata_vector), intent(in) :: this
    integer, intent(in) :: ix, iy, iz, s_level, file_id
    integer :: dims(5)
    dims = shape(this%data)
    write(file_id,*) "dimensions: ", dims(1)
    write(file_id,*) this%data(:,ix,iy,iz,s_level)

  end subroutine

  subroutine dumpForRestartGridVectorStateful(this,  file_name, absolute_path, only_current_value)
    use hdf5, only : HID_T
    use read_write_parallel_h5_mod, only : readOrWriteDatasetStateful, openexistingh5readwrite, closeh5file
    implicit none
    class(griddata_vector), intent(in) :: this
    character(len=*), intent(in) :: file_name, absolute_path
    INTEGER(HID_T) :: file_id, facc_plist
    logical, intent(in) :: only_current_value

    call openExistingH5ReadWrite(file_name, file_id, facc_plist) 
    if (.not. only_current_value) then
      call readOrWriteDatasetStateful(this%data, absolute_path//this%name, this%nx_ny_nz, this%xyz_offset_rank, file_id, "write")
    else 
      call readOrWriteDatasetStateful(this%data(:,:,:,:,1:1), absolute_path//this%name, this%nx_ny_nz, this%xyz_offset_rank, file_id, "write")
    endif
    call CloseH5File(file_id, facc_plist)

  end subroutine

  subroutine writeXDMFAttributeGridVectorStateful(this, xmf_writer, hdf5_fname, field_full_path, write_only_current, suffix)
    use hdf5, only : HID_T
    use read_write_parallel_h5_mod, only : readOrWriteDatasetStateful
    use number_to_string_mod
    use xmf_writer_mod, only : xmf_writer_type
    implicit none
    class(griddata_vector), intent(in) :: this
    type(xmf_writer_type), intent(inout) :: xmf_writer
    CHARACTER(len=*), intent(in) :: hdf5_fname, field_full_path, suffix
    logical, intent(in) :: write_only_current
    integer :: n_stateful_level 

    if (write_only_current) then
      n_stateful_level = 1
    else 
      n_stateful_level = this%stateful_level
    endif
    
    call xmf_writer%write3DVector( hdf5_fname, field_full_path, trim(adjustl(this%name)), this%n_components, n_stateful_level, suffix)
  end subroutine

  subroutine reloadFromDumpGridVectorStateful(this, file_id)
    use hdf5, only : HID_T
    use read_write_parallel_h5_mod, only : readOrWriteDatasetStateful
    implicit none
    class(griddata_vector), intent(inout) :: this
    INTEGER(HID_T), intent(in) :: file_id

    call readOrWriteDatasetStateful(this%data, this%name, this%nx_ny_nz, this%xyz_offset_rank, file_id, "read")

  end subroutine

  subroutine setVectorToValue(this, vec_val, setAll)
    implicit none
    class(griddata_vector), intent(inout) :: this
    real(k_real), dimension(:), intent(in) :: vec_val
    logical, optional, intent(in) :: setAll
    integer :: ix,iy,iz,istateful, max_index

    max_index = 1
    if (present(setAll)) then
      if (setAll) max_index=this%stateful_level
    endif

    do istateful=1,max_index
    do iz=1,this%nz
    do iy=1,this%ny
    do ix=1,this%nx
      this%data(:,ix,iy,iz,istateful) = vec_val
    enddo
    enddo
    enddo
    enddo

  end subroutine

  subroutine setVectorAvgToValue(this, vec_val, setAll)
    implicit none
    class(griddata_vector), intent(inout) :: this
    real(k_real), dimension(:), intent(in) :: vec_val
    logical, optional, intent(in) :: setAll
    integer :: istateful, max_index

    max_index = 1
    if (present(setAll)) then
      if (setAll) max_index=this%stateful_level
    endif

    do istateful=1,max_index
      this%avg_data(:,istateful) = vec_val
    enddo
  end subroutine

  subroutine initRealVector(this, name, dimension_obj, stateful_level)
    implicit none
    class(griddata_real_vector), intent(inout) :: this
    class(var_dimension_type), pointer :: dimension_obj 
    character(len=:), allocatable :: name
    integer, intent(in), optional :: stateful_level
    real(k_real), parameter, dimension(3) :: init_zero=(/0,0,0/)

    this%data_type = real_vector
    this%n_components = 3
    if (present(stateful_level)) call this%setStatefulLevel(stateful_level)
    this%name = name
    call this%setVarDimsFromVarDimsType(dimension_obj)
    allocate(this%data(3, this%nx,this%ny,this%nz, this%stateful_level))
    allocate(this%avg_data(3,this%stateful_level))
    call this%setToValue(init_zero, .true.)
    call this%setAvgToValue(init_zero, .true.)
    this%initialized = .TRUE.
  end subroutine

  subroutine initVector5(this, name, dimension_obj, stateful_level)
    implicit none
    class(griddata_vector5), intent(inout) :: this
    class(var_dimension_type), pointer :: dimension_obj 
    character(len=:), allocatable :: name
    integer, intent(in), optional :: stateful_level
    real(k_real), parameter, dimension(5) :: init_zero=(/0,0,0,0,0/)
    this%data_type = vector5
    this%n_components = 5
    if (present(stateful_level)) call this%setStatefulLevel(stateful_level)
    this%name = name
    call this%setVarDimsFromVarDimsType(dimension_obj)
    allocate(this%data(5, this%nx,this%ny,this%nz, this%stateful_level))
    allocate(this%avg_data(5,this%stateful_level))
    call this%setToValue(init_zero, .true.)
    call this%setAvgToValue(init_zero, .true.)
    this%initialized = .TRUE.
  end subroutine

  subroutine initVector6(this, name, dimension_obj, stateful_level)
    implicit none
    class(griddata_vector6), intent(inout) :: this
    class(var_dimension_type), pointer :: dimension_obj 
    character(len=:), allocatable :: name
    integer, intent(in), optional :: stateful_level
    real(k_real), parameter, dimension(6) :: init_zero=(/0,0,0,0,0,0/)
    this%data_type = vector6
    this%n_components = 6
    if (present(stateful_level)) call this%setStatefulLevel(stateful_level)
    this%name = name
    call this%setVarDimsFromVarDimsType(dimension_obj)
    allocate(this%data(6, this%nx,this%ny,this%nz, this%stateful_level))
    allocate(this%avg_data(6,this%stateful_level))
    call this%setToValue(init_zero, .true.)
    call this%setAvgToValue(init_zero, .true.)
    this%initialized = .TRUE.
  end subroutine

  subroutine initSSScalar(this, name, dimension_obj, n_ss, stateful_level)
    implicit none
    class(griddata_ssscalar), intent(inout) :: this
    class(var_dimension_type), pointer :: dimension_obj 
    character(len=:), allocatable :: name
    integer, intent(in) :: n_ss
    integer, intent(in), optional :: stateful_level
    real(k_real), allocatable, dimension(:) :: init_zero

    this%data_type = ss_scalar
    this%n_components = n_ss
    if (present(stateful_level)) call this%setStatefulLevel(stateful_level)
    this%name = name
    call this%setVarDimsFromVarDimsType(dimension_obj)
    this%n_ss = n_ss
    allocate(this%data(this%n_ss, this%nx, this%ny, this%nz, this%stateful_level))
    allocate(this%avg_data(this%n_ss, this%stateful_level))
    allocate(init_zero(n_ss))
    init_zero = 0
    call this%setToValue(init_zero, .true.)
    call this%setAvgToValue(init_zero, .true.)
    this%initialized = .TRUE.
  end subroutine

  subroutine dumpForRestartGridMatrixStateful(this, file_name, absolute_path, only_current_value)
    use hdf5, only : HID_T
    use read_write_parallel_h5_mod, only : readOrWriteDatasetStateful, openexistingh5readwrite, closeh5file
    implicit none
    class(griddata_matrix), intent(in) :: this
    character(len=*), intent(in) :: file_name, absolute_path
    INTEGER(HID_T) :: file_id, facc_plist
    logical, intent(in) :: only_current_value

    call openExistingH5ReadWrite(file_name, file_id, facc_plist) 
    if (.not. only_current_value) then
      call readOrWriteDatasetStateful(this%data, absolute_path//this%name, this%nx_ny_nz, this%xyz_offset_rank, file_id, "write")
    else 
      call readOrWriteDatasetStateful(this%data(:,:,:,:,:,1:1), absolute_path//this%name, this%nx_ny_nz, this%xyz_offset_rank, file_id, "write")
    endif
    call CloseH5File(file_id, facc_plist)

  end subroutine

  subroutine writeXDMFAttributeGridMatrixStateful(this, xmf_writer, hdf5_fname, field_full_path, write_only_current, suffix)
    use hdf5, only : HID_T
    use xmf_writer_mod, only : xmf_writer_type
    implicit none
    class(griddata_matrix), intent(in) :: this
    type(xmf_writer_type), intent(inout) :: xmf_writer
    CHARACTER(len=*), intent(in) :: hdf5_fname, field_full_path, suffix
    logical, intent(in) :: write_only_current
    integer :: n_stateful_level 

    if (write_only_current) then
      n_stateful_level = 1
    else 
      n_stateful_level = this%stateful_level
    endif

    call xmf_writer%write3DMatrix( hdf5_fname, field_full_path, trim(adjustl(this%name)), this%n_comp_i, this%n_comp_j, n_stateful_level, suffix)

  end subroutine


  subroutine reloadFromDumpGridMatrixStateful(this, file_id)
    use hdf5, only : HID_T
    use read_write_parallel_h5_mod, only : readOrWriteDatasetStateful
    implicit none
    class(griddata_matrix), intent(inout) :: this
    INTEGER(HID_T), intent(in) :: file_id

    call readOrWriteDatasetStateful(this%data, this%name, this%nx_ny_nz, this%xyz_offset_rank, file_id, "read")

  end subroutine

  subroutine writeDataPointToFileMatrix(this, ix,iy,iz, s_level, file_id)
    implicit none
    class(griddata_matrix), intent(in) :: this
    integer, intent(in) :: ix, iy, iz, s_level, file_id
    integer :: j, dims(6)
    dims = shape(this%data)
    write(file_id,*) "dimensions: ", dims(1:2)
    do j = 1,size(this%data,2)
    write(file_id,*) this%data(:,j,ix,iy,iz,s_level)
    enddo
  end subroutine


  subroutine setMatrixToValue(this, mtx_val, setAll)
    implicit none
    class(griddata_matrix), intent(inout) :: this
    real(k_real), dimension(:,:), intent(in) :: mtx_val
    logical, optional, intent(in) :: setAll
    integer :: ix,iy,iz,istateful, max_index

    max_index = 1
    if (present(setAll)) then
      if (setAll) max_index=this%stateful_level
    endif

    do istateful=1,max_index
    do iz=1,this%nz
    do iy=1,this%ny
    do ix=1,this%nx
      this%data(:,:,ix,iy,iz,istateful) = mtx_val
    enddo
    enddo
    enddo
    enddo

  end subroutine

  subroutine setMatrixAvgToValue(this, mtx_val, setAll)
    implicit none
    class(griddata_matrix), intent(inout) :: this
    real(k_real), dimension(:,:), intent(in) :: mtx_val
    logical, optional, intent(in) :: setAll
    integer :: istateful, max_index

    max_index = 1
    if (present(setAll)) then
      if (setAll) max_index=this%stateful_level
    endif

    do istateful=1,max_index
      this%avg_data(:,:,istateful) = mtx_val
    enddo
  end subroutine

  subroutine initGenericMatrix(this, name, dimension_obj, n_comp_i, n_comp_j, stateful_level)
    implicit none
    class(griddata_generic_matrix), intent(inout) :: this
    class(var_dimension_type), pointer :: dimension_obj 
    character(len=:), allocatable :: name
    integer, intent(in) :: n_comp_i, n_comp_j
    integer, intent(in), optional :: stateful_level
    real(k_real), allocatable, dimension(:,:) :: init_zero
    this%data_type = ss_scalar

    if (present(stateful_level)) call this%setStatefulLevel(stateful_level)
    this%name = name
    call this%setVarDimsFromVarDimsType(dimension_obj)
    this%n_comp_i=n_comp_i
    this%n_comp_j=n_comp_j
    allocate(this%data(this%n_comp_i, this%n_comp_j, this%nx, this%ny, this%nz, this%stateful_level))
    allocate(this%avg_data(this%n_comp_i, this%n_comp_j, this%stateful_level))
    allocate(init_zero(n_comp_i,n_comp_j))
    init_zero = 0
    call this%setToValue(init_zero, .true.)
    call this%setAvgToValue(init_zero, .true.)
    this%initialized = .TRUE.
    deallocate(init_zero)
    write(*,*) " initialized ", name, " of type griddata_generic_matrix"
  end subroutine

  subroutine initSSGenericVector(this, name, dimension_obj, n_ss, n_components, stateful_level)
    implicit none
    class(griddata_ss_generic_vector), intent(inout) :: this
    class(var_dimension_type), pointer :: dimension_obj 
    character(len=:), allocatable :: name
    integer, intent(in) :: n_ss, n_components
    integer, intent(in), optional :: stateful_level
    real(k_real), allocatable, dimension(:,:) :: init_zero
    this%data_type = ss_scalar

    if (present(stateful_level)) call this%setStatefulLevel(stateful_level)
    this%name = name
    call this%setVarDimsFromVarDimsType(dimension_obj)
    this%n_ss = n_ss
    this%n_components = n_components
    this%n_comp_i=n_ss
    this%n_comp_j=n_components
    allocate(this%data(this%n_components, this%n_ss, this%nx, this%ny, this%nz, this%stateful_level))
    allocate(this%avg_data(this%n_components, this%n_ss, this%stateful_level))
    allocate(init_zero(n_components,n_ss))
    init_zero = 0
    call this%setToValue(init_zero, .true.)
    call this%setAvgToValue(init_zero, .true.)
    this%initialized = .TRUE.
  end subroutine

  subroutine initSSVector5(this, name, dimension_obj, n_ss, stateful_level)
    implicit none
    class(griddata_ss_vector5), intent(inout) :: this
    class(var_dimension_type), pointer :: dimension_obj 
    character(len=:), allocatable :: name
    integer, intent(in) :: n_ss
    integer, intent(in), optional :: stateful_level
    real(k_real), allocatable, dimension(:,:) :: init_zero

    this%data_type = ss_scalar

    if (present(stateful_level)) call this%setStatefulLevel(stateful_level)
    this%name = name
    call this%setVarDimsFromVarDimsType(dimension_obj)
    this%n_ss = n_ss
    this%n_comp_i=n_ss
    this%n_comp_j=5
    allocate(this%data(5, this%n_ss, this%nx, this%ny, this%nz, this%stateful_level))
    allocate(this%avg_data(5,this%n_ss,this%stateful_level))
    allocate(init_zero(5,n_ss))
    init_zero = 0
    call this%setToValue(init_zero, .true.)
    call this%setAvgToValue(init_zero, .true.)
    this%initialized = .TRUE.
  end subroutine

  subroutine initSSVector6(this, name, dimension_obj, n_ss, stateful_level)
    implicit none
    class(griddata_ss_vector6), intent(inout) :: this
    class(var_dimension_type), pointer :: dimension_obj 
    character(len=:), allocatable :: name
    integer, intent(in) :: n_ss
    integer, intent(in), optional :: stateful_level
    real(k_real), allocatable, dimension(:,:) :: init_zero

    this%data_type = ss_scalar

    if (present(stateful_level)) call this%setStatefulLevel(stateful_level)
    this%name = name
    call this%setVarDimsFromVarDimsType(dimension_obj)
    this%n_ss = n_ss
    this%n_comp_i=n_ss
    this%n_comp_j=6
    allocate(this%data(6, this%n_ss, this%nx, this%ny, this%nz, this%stateful_level))
    allocate(this%avg_data(6,this%n_ss,this%stateful_level))
    allocate(init_zero(6,n_ss))
    init_zero = 0
    call this%setToValue(init_zero, .true.)
    call this%setAvgToValue(init_zero, .true.)
    this%initialized = .TRUE.
  end subroutine

  subroutine initGenericVector(this, name, dimension_obj, n_components, stateful_level)
    implicit none
    class(griddata_generic_vector), intent(inout) :: this
    class(var_dimension_type), pointer :: dimension_obj
    character(len=:), allocatable :: name
    integer, intent(in) :: n_components
    integer, intent(in), optional :: stateful_level
    real(k_real), allocatable, dimension(:) :: init_zero
    this%data_type = generic_vector

    if (present(stateful_level)) call this%setStatefulLevel(stateful_level)
    this%name = name
    call this%setVarDimsFromVarDimsType(dimension_obj)
    this%n_components = n_components
    allocate(this%data(this%n_components, this%nx, this%ny, this%nz, this%stateful_level))
    allocate(this%avg_data(this%n_components,this%stateful_level))
    allocate(init_zero(n_components))
    init_zero = 0
    call this%setToValue(init_zero, .true.)
    call this%setAvgToValue(init_zero, .true.)
    this%initialized = .TRUE.
  end subroutine

  subroutine getDataPointerByStatefulLevelIndexVector(this, stateful_idx, pt_2_data)
    implicit none
    class(griddata_vector), intent(in) :: this
    integer, intent(in) :: stateful_idx
    real(k_real), dimension(:,:,:,:), pointer, intent(inout) :: pt_2_data
    call this%errorOnGetOldOrOlder(stateful_idx)
    pt_2_data => this%data(:,:,:,:,stateful_idx)
  end subroutine

  subroutine getAvgDataPointerByStatefulLevelIndexVector(this, stateful_idx, pt_2_avg_data)
    implicit none
    class(griddata_vector), intent(in) :: this
    integer, intent(in) :: stateful_idx
    real(k_real), dimension(:), pointer, intent(inout) :: pt_2_avg_data
    call this%errorOnGetOldOrOlder(stateful_idx)
    pt_2_avg_data => this%avg_data(:,stateful_idx)
  end subroutine

  subroutine computeAverageVector(this)
    use mpi_useful_routines_mod, only : MPIAverageGridVector
    implicit none
    class(griddata_vector), intent(inout) :: this

    call MPIAverageGridVector(this%data(:,:,:,:,1), this%avg_data(:,1))

  end subroutine


  subroutine computeAverageVectorMasked(this, grid_mask)
    use mpi_useful_routines_mod, only : MPIAverageGridVectorMasked
    implicit none
    class(griddata_vector), intent(inout) :: this
    logical, intent(in), dimension(:,:,:) :: grid_mask

    call MPIAverageGridVectorMasked(this%data(:,:,:,:,1), grid_mask, this%avg_data(:,1))

  end subroutine

  subroutine updateHistoryVector(this)
    implicit none
    class(griddata_vector), intent(inout) :: this
    integer :: i
    if (this%stateful) then
      if (this%stateful_level > 1) then
        do i= this%stateful_level,2,-1
          this%data(:,:,:,:,i) = this%data(:,:,:,:,i-1)
          this%avg_data(:,i) = this%avg_data(:,i-1)
        end do
      end if
    end if
  end subroutine

  subroutine resetHistoryVector(this)
    implicit none
    class(griddata_vector), intent(inout) :: this
    if (this%stateful) then
      if (this%stateful_level > 1) then
        this%data(:,:,:,:,1) = this%data(:,:,:,:,2)
        this%avg_data(:,1) = this%avg_data(:,2)
      end if
    end if
  end subroutine

  subroutine checkGridDataDimensionGenericVector(this)
    implicit none
    class(griddata_generic_vector), intent(in) :: this
    call this%checkDimensionBase()
    if (this%n_components <= 0) then
      write(*,*) "n_components <= 0 for griddata named ", this%name
      stop
    end if
  end subroutine

  subroutine initTensor2(this, name, dimension_obj, stateful_level)
    implicit none
    class(griddata_tensor2), intent(inout) :: this
    class(var_dimension_type), pointer :: dimension_obj 
    character(len=:), allocatable :: name
    integer, intent(in), optional :: stateful_level
    real(k_real), allocatable, dimension(:,:) :: init_zero
    this%data_type = tensor2

    if (present(stateful_level)) call this%setStatefulLevel(stateful_level)
    this%name = name
    call this%setVarDimsFromVarDimsType(dimension_obj)
    this%n_comp_i = 3
    this%n_comp_j=3
    allocate(this%data(3, 3, this%nx,this%ny,this%nz, this%stateful_level))
    allocate(this%avg_data(3, 3,this%stateful_level))
    allocate(init_zero(3,3))
    init_zero = 0
    call this%setToValue(init_zero, .true.)
    call this%setAvgToValue(init_zero, .true.)
    this%initialized = .TRUE.
  end subroutine

  subroutine initMatrix66(this, name, dimension_obj, stateful_level)
    implicit none
    class(griddata_matrix66), intent(inout) :: this
    class(var_dimension_type), pointer :: dimension_obj 
    character(len=:), allocatable :: name
    integer, intent(in), optional :: stateful_level
    real(k_real), allocatable, dimension(:,:) :: init_zero

    this%data_type = matrix66

    if (present(stateful_level)) call this%setStatefulLevel(stateful_level)

    this%name = name
    call this%setVarDimsFromVarDimsType(dimension_obj)
    this%n_comp_i = 6
    this%n_comp_j=6
    allocate(this%data(6, 6, this%nx,this%ny,this%nz, this%stateful_level))
    allocate(this%avg_data(6, 6,this%stateful_level))
    allocate(init_zero(6,6))
    init_zero = 0
    call this%setToValue(init_zero, .true.)
    call this%setAvgToValue(init_zero, .true.)
    this%initialized = .TRUE.
  end subroutine

  subroutine getDataPointerByStatefulLevelIndexMatrix(this, stateful_idx, pt_2_data)
    implicit none
    class(griddata_matrix), intent(in) :: this
    integer, intent(in) :: stateful_idx
    real(k_real), dimension(:,:,:,:,:), pointer, intent(inout) :: pt_2_data
    call this%errorOnGetOldOrOlder(stateful_idx)
    pt_2_data => this%data(:,:,:,:,:,stateful_idx)
  end subroutine

  subroutine getAvgDataPointerByStatefulLevelIndexMatrix(this, stateful_idx, pt_2_avg_data)
    implicit none
    class(griddata_matrix), intent(in) :: this
    integer, intent(in) :: stateful_idx
    real(k_real), dimension(:,:), pointer, intent(inout) :: pt_2_avg_data
    call this%errorOnGetOldOrOlder(stateful_idx)
    pt_2_avg_data => this%avg_data(:,:,stateful_idx)
  end subroutine

  subroutine computeAverageMatrix(this)
    use mpi_useful_routines_mod, only : MPIAverageGridMatrix
    implicit none
    class(griddata_matrix), intent(inout) :: this

    call MPIAverageGridMatrix(this%data(:,:,:,:,:,1), this%avg_data(:,:,1))
  end subroutine

  subroutine computeAverageMatrixMasked(this, grid_mask)
    use mpi_useful_routines_mod, only : MPIAverageGridMatrixMasked
    implicit none
    class(griddata_matrix), intent(inout) :: this
    logical, intent(in), dimension(:,:,:) :: grid_mask

    call MPIAverageGridMatrixMasked(this%data(:,:,:,:,:,1), grid_mask, this%avg_data(:,:,1))
  end subroutine

  subroutine updateHistoryMatrix(this)
    implicit none
    class(griddata_matrix), intent(inout) :: this
    integer :: i
    if (this%stateful) then
      if (this%stateful_level > 1) then
        do i= this%stateful_level,2,-1
          this%data(:,:,:,:,:,i) = this%data(:,:,:,:,:,i-1)
          this%avg_data(:,:,i) = this%avg_data(:,:,i-1)
        end do
      end if
    end if
  end subroutine

  subroutine resetHistoryMatrix(this)
    implicit none
    class(griddata_matrix), intent(inout) :: this
    if (this%stateful) then
      if (this%stateful_level > 1) then
          this%data(:,:,:,:,:,1) = this%data(:,:,:,:,:,2)
          this%avg_data(:,:,1) = this%avg_data(:,:,2)
      end if
    end if
  end subroutine

  subroutine dumpForRestartGridSSMatrixStateful(this, file_name, absolute_path, only_current_value)
    use hdf5, only : HID_T
    use read_write_parallel_h5_mod, only : readOrWriteDatasetStateful, openexistingh5readwrite, closeh5file
    implicit none
    class(griddata_ss_generic_matrix), intent(in) :: this
    character(len=*), intent(in) :: file_name, absolute_path
    INTEGER(HID_T) :: file_id, facc_plist
    logical, intent(in) :: only_current_value

    call openExistingH5ReadWrite(file_name, file_id, facc_plist) 
    if (.not. only_current_value) then
      call readOrWriteDatasetStateful(this%data, absolute_path//this%name, this%nx_ny_nz, this%xyz_offset_rank, file_id, "write")
    else 
      call readOrWriteDatasetStateful(this%data(:,:,:,:,:,:,1:1), absolute_path//this%name, this%nx_ny_nz, this%xyz_offset_rank, file_id, "write")
    endif
    call CloseH5File(file_id, facc_plist)

  end subroutine

  subroutine writeXDMFAttributeGridSSMatrix(this, xmf_writer, hdf5_fname, field_full_path, write_only_current, suffix)
    use xmf_writer_mod, only : xmf_writer_type
    implicit none
    class(griddata_ss_generic_matrix), intent(in) :: this
    type(xmf_writer_type), intent(inout) :: xmf_writer
    CHARACTER(len=*), intent(in) :: hdf5_fname, field_full_path, suffix
    logical, intent(in) :: write_only_current
    integer :: n_stateful_level 
    __DECL_UNUSED_CHARACTER__
   

    __SUPPRESS_UNUSED_STRING__(hdf5_fname)
    __SUPPRESS_UNUSED_STRING__(field_full_path)
    __SUPPRESS_UNUSED_STRING__(suffix)


    if (write_only_current) then
      n_stateful_level = 1
    else 
      n_stateful_level = this%stateful_level
    endif

    write(*,*) "skikpping for now"
  end subroutine


  subroutine reloadFromDumpGridSSMatrixStateful(this, file_id)
    use hdf5, only : HID_T
    use read_write_parallel_h5_mod, only : readOrWriteDatasetStateful
    implicit none
    class(griddata_ss_generic_matrix), intent(inout) :: this
    INTEGER(HID_T), intent(in) :: file_id

    call readOrWriteDatasetStateful(this%data, this%name, this%nx_ny_nz, this%xyz_offset_rank, file_id, "read")

  end subroutine

  subroutine initSSMatrix(this, name, dimension_obj, n_ss, n_components, m_components, stateful_level)
    implicit none
    class(griddata_ss_generic_matrix), intent(inout) :: this
    class(var_dimension_type), pointer :: dimension_obj 
    character(len=:), allocatable :: name
    integer, intent(in) :: n_ss, n_components, m_components
    integer, intent(in), optional :: stateful_level
    real(k_real), allocatable, dimension(:,:,:) :: init_zero
    this%data_type = ss_scalar

    if (present(stateful_level)) call this%setStatefulLevel(stateful_level)
    this%name = name
    call this%setVarDimsFromVarDimsType(dimension_obj)
    this%n_ss = n_ss
    this%n_components = n_components
    this%m_components = m_components
    allocate(this%data(this%n_components, this%m_components, this%n_ss, this%nx, this%ny, this%nz, this%stateful_level))
    allocate(this%avg_data(this%n_components, this%m_components, this%n_ss, this%stateful_level))
    allocate(init_zero(n_components, m_components, n_ss))
    init_zero = 0
    call this%setToValue(init_zero, .true.)
    call this%setAvgToValue(init_zero, .true.)
    this%initialized = .TRUE.
  end subroutine

  subroutine writeDataPointToFileSSMatrix(this, ix,iy,iz, s_level, file_id)
    implicit none
    class(griddata_ss_generic_matrix), intent(in) :: this
    integer, intent(in) :: ix, iy, iz, s_level, file_id
    integer :: j,k, dims(7)
    dims = shape(this%data)
    write(file_id,*) "dimensions: ", dims(1:3)
    do k = 1,size(this%data,3)
    do j = 1,size(this%data,2)
    write(file_id,*) this%data(:,j,k,ix,iy,iz,s_level)
    enddo
    enddo
  end subroutine

  subroutine getDataPointerByStatefulLevelIndexSSMatrix(this, stateful_idx, pt_2_data)
    implicit none
    class(griddata_ss_generic_matrix), intent(in) :: this
    integer, intent(in) :: stateful_idx
    real(k_real), dimension(:,:,:,:,:,:), pointer, intent(inout) :: pt_2_data
    call this%errorOnGetOldOrOlder(stateful_idx)
    pt_2_data => this%data(:,:,:,:,:,:,stateful_idx)
  end subroutine

  subroutine getAvgDataPointerByStatefulLevelIndexSSMatrix(this, stateful_idx, pt_2_avg_data)
    implicit none
    class(griddata_ss_generic_matrix), intent(in) :: this
    integer, intent(in) :: stateful_idx
    real(k_real), dimension(:,:,:), pointer, intent(inout) :: pt_2_avg_data
    call this%errorOnGetOldOrOlder(stateful_idx)
    pt_2_avg_data => this%avg_data(:,:,:,stateful_idx)
  end subroutine

  subroutine computeAverageSSMatrix(this)
    use mpi_useful_routines_mod, only : MPIAverageGridSSMatrix
    implicit none
    class(griddata_ss_generic_matrix), intent(inout) :: this

    call MPIAverageGridSSMatrix(this%data(:,:,:,:,:,:,1),this%avg_data(:,:,:,1))

  end subroutine

  subroutine updateHistorySSMatrix(this)
    implicit none
    class(griddata_ss_generic_matrix), intent(inout) :: this
    integer :: i
    if (this%stateful) then
      if (this%stateful_level > 1) then
        do i= this%stateful_level,2,-1
          this%data(:,:,:,:,:,:,i) = this%data(:,:,:,:,:,:,i-1)
          this%avg_data(:,:,:,i) = this%avg_data(:,:,:,i-1)
        end do
      end if
    end if
  end subroutine

  subroutine resetHistorySSMatrix(this)
    implicit none
    class(griddata_ss_generic_matrix), intent(inout) :: this
    if (this%stateful) then
      if (this%stateful_level > 1) then
        this%data(:,:,:,:,:,:,1) = this%data(:,:,:,:,:,:,2)
        this%avg_data(:,:,:,1) = this%avg_data(:,:,:,2)
      end if
    end if
  end subroutine

  subroutine setSSMatrixToValue(this, val, setAll)
    implicit none
    class(griddata_ss_generic_matrix), intent(inout) :: this
    real(k_real), intent(in) :: val(:,:,:)
    logical, optional, intent(in) :: setAll
    integer :: ix,iy,iz,istateful, max_index

    max_index = 1
    if (present(setAll)) then
      if (setAll) max_index=this%stateful_level
    endif

    do istateful=1,max_index
    do iz=1,this%nz
    do iy=1,this%ny
    do ix=1,this%nx
      this%data(:,:,:,ix,iy,iz,istateful) = val
    enddo
    enddo
    enddo
    enddo

  end subroutine

  subroutine setSSMatrixAvgToValue(this, val, setAll)
    implicit none
    class(griddata_ss_generic_matrix), intent(inout) :: this
    real(k_real), intent(in) :: val(:,:,:)
    logical, optional, intent(in) :: setAll
    integer :: istateful, max_index

    max_index = 1
    if (present(setAll)) then
      if (setAll) max_index=this%stateful_level
    endif

    do istateful=1,max_index
      this%avg_data(:,:,:,istateful) = val
    enddo

  end subroutine

end module
