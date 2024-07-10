module read_write_parallel_h5_mod
  use kinds
  use hdf5
  implicit none

  private :: readOrwrite3dstateful, &
              readOrWrite3DStatefulReal, &
              readOrWrite3DStatefulInteger, &
              readOrwrite4dstateful, &
              readOrWrite4DStatefulReal, &
              readOrWrite4DStatefulInteger, &
              readOrwrite5dstateful, &
              readOrWrite5DStatefulReal, &
              readOrWrite5DStatefulInteger, &
              readOrwrite6dstateful, &
              readOrWrite6DStatefulReal, &
              readOrWrite6DStatefulInteger, &
              readOrWriteStateful, &
              findArraySizesStateful, &
              getDataType, &
              readOrWriteDatasetPreamble, &
              readOrWriteDatasetClosing

  public :: readOrWriteDatasetStateful, &
            createAndOpenH5EmptyFile, &
            openExistingH5ReadOnly, &
            openExistingH5ReadWrite, &
            CloseH5File

  interface readOrWrite3DStateful
    procedure readOrWrite3DStatefulReal
    procedure readOrWrite3DStatefulInteger
  end interface

  interface readOrWrite4DStateful
    procedure readOrWrite4DStatefulReal
    procedure readOrWrite4DStatefulInteger
  end interface

  interface readOrWrite5DStateful
    procedure readOrWrite5DStatefulReal
    procedure readOrWrite5DStatefulInteger
  end interface

  interface readOrWrite6DStateful
    procedure readOrWrite6DStatefulReal
    procedure readOrWrite6DStatefulInteger
  end interface

  interface readOrWriteStateful
    procedure readOrWrite3DStatefulReal
    procedure readOrWrite3DStatefulInteger
    procedure readOrWrite4DStatefulReal
    procedure readOrWrite4DStatefulInteger
    procedure readOrWrite5DStatefulReal
    procedure readOrWrite5DStatefulInteger
    procedure readOrWrite6DStatefulReal
    procedure readOrWrite6DStatefulInteger
  end interface

  interface readOrWriteDatasetStateful
    procedure readWrite3DRealDatasetStateful
    procedure readWrite3DIntegerDatasetStateful
    procedure readWrite4DRealDatasetStateful
    procedure readWrite4DIntegerDatasetStateful
    procedure readWrite5DRealDatasetStateful
    procedure readWrite5DIntegerDatasetStateful
    procedure readWrite6DRealDatasetStateful
    procedure readWrite6DIntegerDatasetStateful
  end interface

contains
  
  subroutine createAndOpenH5EmptyFile(file_name, file_id, file_acces_prop_list)
    use mpi, only : MPI_COMM_WORLD, MPI_INFO_NULL
    implicit none
    character(len=*), intent(in) :: file_name
    INTEGER(HID_T), intent(out) :: file_id, file_acces_prop_list
    integer :: ierr

    CALL h5pcreate_f(H5P_FILE_ACCESS_F, file_acces_prop_list, ierr)
    if (ierr.ne.0) error stop "createAndOpenH5EmptyFile: error creating the property list"
    CALL h5pset_fapl_mpio_f(file_acces_prop_list, MPI_COMM_WORLD, MPI_INFO_NULL, ierr)
    if (ierr.ne.0) error stop "createAndOpenH5EmptyFile: error setting the property list"
    CALL h5fcreate_f(file_name, H5F_ACC_TRUNC_F, file_id, ierr, access_prp = file_acces_prop_list)
    if (ierr.ne.0) error stop "createAndOpenH5EmptyFile: error creating the file"
  end subroutine

  subroutine createGroup(group_name, file_id)
    use mpi, only : MPI_COMM_WORLD, MPI_INFO_NULL
    implicit none
    character(len=*), intent(in) :: group_name
    INTEGER(HID_T), intent(in) :: file_id
    INTEGER(HID_T) :: group_id
    integer :: ierr

    call h5gcreate_f(file_id, group_name, group_id, ierr)
    if (ierr.ne.0) error stop "createGroup: error creating group"
    call h5gclose_f(group_id, ierr)
    if (ierr.ne.0) error stop "createGroup: error closing group"
  end subroutine
  
  subroutine openExistingH5ReadOnly(file_name, file_id, file_acces_prop_list)
    use mpi, only : MPI_COMM_WORLD, MPI_INFO_NULL
    implicit none
    character(len=*), intent(in) :: file_name
    INTEGER(HID_T), intent(out) :: file_id, file_acces_prop_list
    integer :: ierr

    CALL h5pcreate_f(H5P_FILE_ACCESS_F, file_acces_prop_list, ierr)
    if (ierr.ne.0) error stop "openExistingH5ReadOnly: error creating the property list"
    CALL h5pset_fapl_mpio_f(file_acces_prop_list, MPI_COMM_WORLD, MPI_INFO_NULL, ierr)
    if (ierr.ne.0) error stop "openExistingH5ReadOnly: error setting the property list"
    CALL h5fopen_f(file_name, H5F_ACC_RDONLY_F, file_id, ierr, access_prp = file_acces_prop_list)
    if (ierr.ne.0) error stop "openExistingH5ReadOnly: error opening the file"
  end subroutine

  subroutine openExistingH5ReadWrite(file_name, file_id, file_acces_prop_list)
    use mpi, only : MPI_COMM_WORLD, MPI_INFO_NULL
    implicit none
    character(len=*), intent(in) :: file_name
    INTEGER(HID_T), intent(out) :: file_id, file_acces_prop_list
    integer :: ierr

    CALL h5pcreate_f(H5P_FILE_ACCESS_F, file_acces_prop_list, ierr)
    if (ierr.ne.0) error stop "openExistingH5ReadWrite: error creating the property list"
    CALL h5pset_fapl_mpio_f(file_acces_prop_list, MPI_COMM_WORLD, MPI_INFO_NULL, ierr)
    if (ierr.ne.0) error stop "openExistingH5ReadWrite: error setting the property list"
    CALL h5fopen_f(file_name, H5F_ACC_RDWR_F, file_id, ierr, access_prp = file_acces_prop_list)
    if (ierr.ne.0) error stop "openExistingH5ReadWrite: error opening the file"
  end subroutine

  subroutine CloseH5File(file_id, file_acces_prop_list)
    implicit none
    INTEGER(HID_T), intent(in) :: file_id, file_acces_prop_list
    integer :: ierr

    CALL h5pclose_f(file_acces_prop_list, ierr)
    if (ierr.ne.0) error stop "CloseH5File: closing file access property list"
    CALL h5fclose_f(file_id, ierr)
    if (ierr.ne.0) error stop "CloseH5File: closing file"
  end subroutine

  subroutine findArraySizesStateful(nx_ny_nz, xyz_offset, data_shape, array_size, array_subsize, array_start)
    INTEGER, intent(in), dimension(3) ::  nx_ny_nz, xyz_offset
    INTEGER, intent(in), dimension(:) ::  data_shape
    INTEGER(HSIZE_T), dimension(:), intent(out) :: array_size, array_subsize, array_start
    integer :: z_dim
    z_dim = size(data_shape) -1
    array_subsize = data_shape
    array_size = array_subsize; array_size(z_dim-2:z_dim) = nx_ny_nz
    array_start = 0; array_start(z_dim-2:z_dim)= xyz_offset
  end subroutine

  subroutine getDataType(integer_real, dtype_id)
    character(len=*), intent(in) :: integer_real
    INTEGER(HID_T), intent(out) :: dtype_id

    ! select data type
    select case (integer_real)
    case ("integer")
      dtype_id= H5T_NATIVE_INTEGER
    case ("real")
      dtype_id= H5T_NATIVE_DOUBLE
    case default
      error stop "readOrWriteInit read_write can only be integer or real"
    end select
  endsubroutine

  subroutine readOrWrite3DStatefulReal(array_subsize, dset_id, filespace, memspace, p_access_list, data, read_write)
    INTEGER(HID_T), intent(in) :: dset_id, filespace, memspace, p_access_list
    real(k_real), dimension(:,:,:,:), intent(inout) :: data
    character(len=*), intent(in) :: read_write
    INTEGER(HSIZE_T), intent(in) :: array_subsize(4)
    integer :: ierr
    INTEGER(HID_T) :: dtype_id
    call getDataType("real", dtype_id)
    select case (read_write)
      case ("read")
        call h5dread_f(dset_id, dtype_id, data, array_subsize, ierr, &
                       mem_space_id=memspace, file_space_id=filespace, xfer_prp=p_access_list)
        if (ierr.ne.0) error stop "readOrWrite3DStatefulReal h5dread_f"
      case ("write")
        call h5dwrite_f(dset_id, dtype_id, data, array_subsize, ierr, memspace, filespace, p_access_list)
        if (ierr.ne.0) error stop "readOrWrite3DStatefulReal h5dwrite_f"
      case default
        error stop "readOrWrite3DStatefulReal wrong read_write option"
    end select
  end subroutine

  subroutine readOrWrite3DStatefulInteger(array_subsize, dset_id, filespace, memspace, p_access_list, data, read_write)
    INTEGER(HID_T), intent(in) :: dset_id, filespace, memspace, p_access_list
    integer, dimension(:,:,:,:), intent(inout) :: data
    character(len=*), intent(in) :: read_write
    INTEGER(HSIZE_T), intent(in) :: array_subsize(4)
    integer :: ierr
    INTEGER(HID_T) :: dtype_id
    call getDataType("integer", dtype_id)
    select case (read_write)
      case ("read")
        call h5dread_f(dset_id, dtype_id, data, array_subsize, ierr, &
                       mem_space_id=memspace, file_space_id=filespace, xfer_prp=p_access_list)
        if (ierr.ne.0) error stop "readOrWrite3DStatefulInteger h5dread_f"
      case ("write")
        call h5dwrite_f(dset_id, dtype_id, data, array_subsize, ierr, memspace, filespace, p_access_list)
        if (ierr.ne.0) error stop "readOrWrite3DStatefulInteger h5dwrite_f"
      case default
        error stop "readOrWrite3DStatefulInteger wrong read_write option"
    end select
  end subroutine

  subroutine readOrWrite4DStatefulReal(array_subsize, dset_id, filespace, memspace, p_access_list, data, read_write)
    INTEGER(HID_T), intent(in) :: dset_id, filespace, memspace, p_access_list
    real(k_real), dimension(:,:,:,:,:), intent(inout) :: data
    character(len=*), intent(in) :: read_write
    INTEGER(HSIZE_T), intent(in) :: array_subsize(5)
    integer :: ierr
    INTEGER(HID_T) :: dtype_id
    call getDataType("real", dtype_id)
    select case (read_write)
      case ("read")
        call h5dread_f(dset_id, dtype_id, data, array_subsize, ierr, &
                       mem_space_id=memspace, file_space_id=filespace, xfer_prp=p_access_list)
        if (ierr.ne.0) error stop "readOrWrite3DStatefulInteger h5dread_f"
      case ("write")
        call h5dwrite_f(dset_id, dtype_id, data, array_subsize, ierr, memspace, filespace, p_access_list)
        if (ierr.ne.0) error stop "readOrWrite4DStatefulReal h5dwrite_f"
      case default
        error stop "readOrWrite4DStatefulReal wrong read_write option"
    end select
  end subroutine

  subroutine readOrWrite4DStatefulInteger(array_subsize, dset_id, filespace, memspace, p_access_list, data, read_write)
    INTEGER(HID_T), intent(in) :: dset_id, filespace, memspace, p_access_list
    integer, dimension(:,:,:,:,:), intent(inout) :: data
    character(len=*), intent(in) :: read_write
    INTEGER(HSIZE_T), intent(in) :: array_subsize(5)
    integer :: ierr
    INTEGER(HID_T) :: dtype_id
    call getDataType("integer", dtype_id)
    select case (read_write)
      case ("read")
        call h5dread_f(dset_id, dtype_id, data, array_subsize, ierr, &
                       mem_space_id=memspace, file_space_id=filespace, xfer_prp=p_access_list)
        if (ierr.ne.0) error stop "readOrWrite4DStatefulInteger h5dread_f"
      case ("write")
        call h5dwrite_f(dset_id, dtype_id, data, array_subsize, ierr, memspace, filespace, p_access_list)
        if (ierr.ne.0) error stop "readOrWrite4DStatefulInteger h5dwrite_f"
      case default
        error stop "readOrWrite4DStatefulInteger wrong read_write option"
    end select
  end subroutine

  subroutine readOrWrite5DStatefulReal(array_subsize, dset_id, filespace, memspace, p_access_list, data, read_write)
    INTEGER(HID_T), intent(in) :: dset_id, filespace, memspace, p_access_list
    real(k_real), dimension(:,:,:,:,:,:), intent(inout) :: data
    character(len=*), intent(in) :: read_write
    INTEGER(HSIZE_T), intent(in) :: array_subsize(6)
    integer :: ierr
    INTEGER(HID_T) :: dtype_id
    call getDataType("real", dtype_id)
    select case (read_write)
      case ("read")
        call h5dread_f(dset_id, dtype_id, data, array_subsize, ierr, &
                       mem_space_id=memspace, file_space_id=filespace, xfer_prp=p_access_list)
        if (ierr.ne.0) error stop "readOrWrite5DStatefulReal h5dread_f"
      case ("write")
        call h5dwrite_f(dset_id, dtype_id, data, array_subsize, ierr, memspace, filespace, p_access_list)
        if (ierr.ne.0) error stop "readOrWrite5DStatefulReal h5dwrite_f"
      case default
        error stop "readOrWrite5DStatefulReal wrong read_write option"
    end select
  end subroutine

  subroutine readOrWrite5DStatefulInteger(array_subsize, dset_id, filespace, memspace, p_access_list, data, read_write)
    INTEGER(HID_T), intent(in) :: dset_id, filespace, memspace, p_access_list
    integer, dimension(:,:,:,:,:,:), intent(inout) :: data
    character(len=*), intent(in) :: read_write
    INTEGER(HSIZE_T), intent(in) :: array_subsize(6)
    integer :: ierr
    INTEGER(HID_T) :: dtype_id
    call getDataType("integer", dtype_id)
    select case (read_write)
      case ("read")
        call h5dread_f(dset_id, dtype_id, data, array_subsize, ierr, &
                       mem_space_id=memspace, file_space_id=filespace, xfer_prp=p_access_list)
        if (ierr.ne.0) error stop "readOrWrite5DStatefulInteger h5dread_f"
      case ("write")
        call h5dwrite_f(dset_id, dtype_id, data, array_subsize, ierr, &
        mem_space_id=memspace, file_space_id=filespace, xfer_prp=p_access_list)
        if (ierr.ne.0) error stop "readOrWrite5DStatefulInteger h5dwrite_f"
      case default
        error stop "readOrWrite5DStatefulInteger wrong read_write option"
    end select
  end subroutine

  subroutine readOrWrite6DStatefulReal(array_subsize, dset_id, filespace, memspace, p_access_list, data, read_write)
    INTEGER(HID_T), intent(in) :: dset_id, filespace, memspace, p_access_list
    real(k_real), dimension(:,:,:,:,:,:,:), intent(inout) :: data
    character(len=*), intent(in) :: read_write
    INTEGER(HSIZE_T), intent(in) :: array_subsize(7)
    integer :: ierr
    INTEGER(HID_T) :: dtype_id
    call getDataType("real", dtype_id)
    select case (read_write)
      case ("read")
        call h5dread_f(dset_id, dtype_id, data, array_subsize, ierr, &
                       mem_space_id=memspace, file_space_id=filespace, xfer_prp=p_access_list)
        if (ierr.ne.0) error stop "h5dread_f readOrWrite6DStatefulReal"
      case ("write")
        call h5dwrite_f(dset_id, dtype_id, data, array_subsize, ierr, memspace, filespace, p_access_list)
        if (ierr.ne.0) error stop "h5dwrite_f readOrWrite6DStatefulReal"
      case default
        error stop "readOrWrite6DStatefulReal wrong read_write option"
    end select
  end subroutine

  subroutine readOrWrite6DStatefulInteger(array_subsize, dset_id, filespace, memspace, p_access_list, data, read_write)
    INTEGER(HID_T), intent(in) :: dset_id, filespace, memspace, p_access_list
    integer, dimension(:,:,:,:,:,:,:), intent(inout) :: data
    character(len=*), intent(in) :: read_write
    INTEGER(HSIZE_T), intent(in) :: array_subsize(7)
    integer :: ierr
    INTEGER(HID_T) :: dtype_id
    call getDataType("integer", dtype_id)
    select case (read_write)
      case ("read")
        call h5dread_f(dset_id, dtype_id, data, array_subsize, ierr, &
                       mem_space_id=memspace, file_space_id=filespace, xfer_prp=p_access_list)
        if (ierr.ne.0) error stop "h5dread_f readOrWrite6DStatefulInteger"
      case ("write")
        call h5dwrite_f(dset_id, dtype_id, data, array_subsize, ierr, memspace, filespace, p_access_list)
        if (ierr.ne.0) error stop "h5dwrite_f readOrWrite6DStatefulInteger"
      case default
        error stop "readOrWrite6DStatefulInteger wrong read_write option"
    end select
  end subroutine

  subroutine readOrWriteDatasetPreamble(dataset_name, file_id, read_write, dtype_id, data_rank, array_subsize, array_size, array_start, p_access_list, memspace, filespace, dset_id)
    implicit none
    CHARACTER(len=*), intent(in) :: dataset_name
    integer, intent(in) :: data_rank
    INTEGER(HID_T), intent(in) :: file_id
    character(len=*), intent(in) :: read_write
    INTEGER(HSIZE_T), dimension(:) :: array_subsize(:), array_size(:), array_start(:)
    INTEGER(HSIZE_T) :: chunk_size(data_rank), max_chunk_size(data_rank)
    INTEGER(HID_T), intent(out) :: p_access_list, memspace, filespace, dset_id
    INTEGER(HID_T), intent(in) :: dtype_id
    INTEGER(HID_T) :: p_create_list
    integer :: ierr, i

    max_chunk_size(data_rank) = 1
    max_chunk_size(data_rank-3:data_rank-1) = 3
    if (data_rank>4) &
       max_chunk_size(1:data_rank-4) = array_size(1:data_rank-4)

    do i = 1, data_rank
      chunk_size(i) = min(array_size(i), max_chunk_size(i))
    enddo

    ! write(*,*) "chunk_size ", chunk_size
    select case (read_write)
      case ("read")
        CALL h5pcreate_f(H5P_DATASET_XFER_F, p_access_list, ierr)
        CALL h5pset_dxpl_mpio_f(p_access_list, H5FD_MPIO_collective_F, ierr)
        CALL h5screate_simple_f(data_rank, array_subsize, memspace, ierr)
        call h5dopen_f(file_id, dataset_name, dset_id, ierr)
        call h5dget_space_f(dset_id, filespace, ierr)
        call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, array_start, array_subsize, ierr)
      case ("write")
        CALL h5pcreate_f(H5P_DATASET_XFER_F, p_access_list, ierr)
        if (ierr.ne.0) error stop "readOrWriteDatasetPreamble h5pcreate_f p_access_list"
        CALL h5pset_dxpl_mpio_f(p_access_list, H5FD_MPIO_collective_F, ierr)
        if (ierr.ne.0) error stop "readOrWriteDatasetPreamble h5pset_dxpl_mpio_f H5FD_MPIO_COLLECTIVE"
        CALL h5screate_simple_f(data_rank, array_subsize, memspace, ierr)
        if (ierr.ne.0) error stop "readOrWriteDatasetPreamble h5screate_simple_f"
        CALL h5pcreate_f(H5P_DATASET_CREATE_F, p_create_list, ierr)
        if (ierr.ne.0) error stop "readOrWriteDatasetPreamble h5pcreate_f"
        CALL h5pset_chunk_f(p_create_list, data_rank, chunk_size, ierr)
        if (ierr.ne.0) error stop "readOrWriteDatasetPreamble h5pset_chunk_f"
        CALL h5pset_deflate_f(p_create_list, 6, ierr)
        if (ierr.ne.0) error stop "readOrWriteDatasetPreamble h5pset_deflate_f"
        CALL h5screate_simple_f(data_rank, array_size, filespace, ierr)
        if (ierr.ne.0) error stop "readOrWriteDatasetPreamble h5screate_simple_f"
        CALL h5dcreate_f(file_id, dataset_name, dtype_id, filespace, dset_id, ierr, dcpl_id=p_create_list)
        if (ierr.ne.0) error stop "readOrWriteDatasetPreamble h5dcreate_f"
        CALL h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, array_start, array_subsize, ierr)
        if (ierr.ne.0) error stop "readOrWriteDatasetPreamble h5sselect_hyperslab_f"
        CALL h5pclose_f(p_create_list, ierr)
        if (ierr.ne.0) error stop "readOrWriteDatasetPreamble h5pclose_f"
      case default
        error stop "readOrWriteDatasetPreamble wrong read_write option"
    end select

  end subroutine

  subroutine readOrWriteDatasetClosing(dset_id, filespace, memspace, p_access_list)
    INTEGER(HID_T), intent(inout) :: p_access_list, memspace, filespace, dset_id
    integer :: ierr
    CALL h5dclose_f(dset_id, ierr)
    if (ierr.ne.0) error stop "readOrWriteDatasetClosing h5dclose_f"
    call h5sclose_f(filespace, ierr)
    if (ierr.ne.0) error stop "readOrWriteDatasetClosing h5sclose_f"
    call h5sclose_f(memspace, ierr)
    if (ierr.ne.0) error stop "readOrWriteDatasetClosing h5sclose_f"
    call h5pclose_f(p_access_list, ierr)
    if (ierr.ne.0) error stop "readOrWriteDatasetClosing h5pclose_f"

  end subroutine

  subroutine readWrite3DRealDatasetStateful(data, dataset_name, nx_ny_nz, xyz_offset, file_id, read_write)
    real(k_real), intent(inout), dimension(:,:,:,:) :: data
    CHARACTER(len=*), intent(in) :: dataset_name
    INTEGER, intent(in), dimension(3) :: nx_ny_nz, xyz_offset
    INTEGER(HID_T), intent(in) :: file_id
    character(len=*), intent(in) :: read_write
    INTEGER(HID_T) :: filespace, dset_id, memspace, p_access_list
    INTEGER(HSIZE_T) :: array_size(rank(data)), array_subsize(rank(data)), array_start(rank(data))
    integer :: data_rank
    INTEGER(HID_T) :: dtype_id

    call getDataType("real", dtype_id)
    call findArraySizesStateful( nx_ny_nz, xyz_offset, shape(data),  array_size, array_subsize, array_start)
    data_rank = rank(data)
    call readOrWriteDatasetPreamble(dataset_name, file_id, read_write, dtype_id, data_rank, array_subsize, array_size, array_start, p_access_list, memspace, filespace, dset_id)
    call readOrWriteStateful(array_subsize, dset_id, filespace, memspace, p_access_list, data, read_write)
    call readOrWriteDatasetClosing(dset_id, filespace, memspace, p_access_list)
  end subroutine

  subroutine readWrite3DIntegerDatasetStateful(data, dataset_name,  nx_ny_nz, xyz_offset, file_id, read_write)
    integer, intent(inout), dimension(:,:,:,:) :: data
    CHARACTER(len=*), intent(in) :: dataset_name
    INTEGER, intent(in), dimension(3) :: nx_ny_nz, xyz_offset
    INTEGER(HID_T), intent(in) :: file_id
    character(len=*), intent(in) :: read_write
    INTEGER(HID_T) :: filespace, dset_id, memspace, p_access_list
    INTEGER(HSIZE_T) :: array_size(rank(data)), array_subsize(rank(data)), array_start(rank(data))
    integer :: data_rank
    INTEGER(HID_T) :: dtype_id

    call getDataType("integer", dtype_id)
    call findArraySizesStateful( nx_ny_nz, xyz_offset, shape(data),  array_size, array_subsize, array_start)
    data_rank = rank(data)
    call readOrWriteDatasetPreamble(dataset_name, file_id, read_write, dtype_id, data_rank, array_subsize, array_size, array_start, p_access_list, memspace, filespace, dset_id)
    call readOrWriteStateful(array_subsize, dset_id, filespace, memspace, p_access_list, data, read_write)
    call readOrWriteDatasetClosing(dset_id, filespace, memspace, p_access_list)
  end subroutine

  subroutine readWrite4DRealDatasetStateful(data, dataset_name, nx_ny_nz, xyz_offset, file_id, read_write)
    real(k_real), intent(inout), dimension(:,:,:,:,:) :: data
    CHARACTER(len=*), intent(in) :: dataset_name
    INTEGER, intent(in), dimension(3) :: nx_ny_nz, xyz_offset
    INTEGER(HID_T), intent(in) :: file_id
    character(len=*), intent(in) :: read_write
    INTEGER(HID_T) :: filespace, dset_id, memspace, p_access_list
    INTEGER(HSIZE_T) :: array_size(rank(data)), array_subsize(rank(data)), array_start(rank(data))
    integer :: data_rank
    INTEGER(HID_T) :: dtype_id

    call getDataType("real", dtype_id)
    call findArraySizesStateful(nx_ny_nz, xyz_offset, shape(data),  array_size, array_subsize, array_start)
    data_rank = rank(data)
    call readOrWriteDatasetPreamble(dataset_name, file_id, read_write, dtype_id, data_rank, array_subsize, array_size, array_start, p_access_list, memspace, filespace, dset_id)
    call readOrWriteStateful(array_subsize, dset_id, filespace, memspace, p_access_list, data, read_write)
    call readOrWriteDatasetClosing(dset_id, filespace, memspace, p_access_list)
  end subroutine

  subroutine readWrite4DIntegerDatasetStateful(data, dataset_name, nx_ny_nz, xyz_offset, file_id, read_write)
    integer, intent(inout), dimension(:,:,:,:,:) :: data
    CHARACTER(len=*), intent(in) :: dataset_name
    INTEGER, intent(in), dimension(3) :: nx_ny_nz, xyz_offset
    INTEGER(HID_T), intent(in) :: file_id
    character(len=*), intent(in) :: read_write
    INTEGER(HID_T) :: filespace, dset_id, memspace, p_access_list
    INTEGER(HSIZE_T) :: array_size(rank(data)), array_subsize(rank(data)), array_start(rank(data))
    integer :: data_rank
    INTEGER(HID_T) :: dtype_id

    call getDataType("integer", dtype_id)
    call findArraySizesStateful(nx_ny_nz, xyz_offset, shape(data),  array_size, array_subsize, array_start)
    data_rank = rank(data)
    call readOrWriteDatasetPreamble(dataset_name, file_id, read_write, dtype_id, data_rank, array_subsize, array_size, array_start, p_access_list, memspace, filespace, dset_id)
    call readOrWriteStateful(array_subsize, dset_id, filespace, memspace, p_access_list, data, read_write)
    call readOrWriteDatasetClosing(dset_id, filespace, memspace, p_access_list)
  end subroutine

  subroutine readWrite5DRealDatasetStateful(data, dataset_name, nx_ny_nz, xyz_offset, file_id, read_write)
    real(k_real), intent(inout), dimension(:,:,:,:,:,:) :: data
    CHARACTER(len=*), intent(in) :: dataset_name
    INTEGER, intent(in), dimension(3) :: nx_ny_nz, xyz_offset
    INTEGER(HID_T), intent(in) :: file_id
    character(len=*), intent(in) :: read_write
    INTEGER(HID_T) :: filespace, dset_id, memspace, p_access_list
    INTEGER(HSIZE_T) :: array_size(rank(data)), array_subsize(rank(data)), array_start(rank(data))
    integer :: data_rank
    INTEGER(HID_T) :: dtype_id

    call getDataType("real", dtype_id)
    call findArraySizesStateful(nx_ny_nz, xyz_offset, shape(data),  array_size, array_subsize, array_start)
    data_rank = rank(data)
    call readOrWriteDatasetPreamble(dataset_name, file_id, read_write, dtype_id, data_rank, array_subsize, array_size, array_start, p_access_list, memspace, filespace, dset_id)
    call readOrWriteStateful(array_subsize, dset_id, filespace, memspace, p_access_list, data, read_write)
    call readOrWriteDatasetClosing(dset_id, filespace, memspace, p_access_list)
  end subroutine

  subroutine readWrite5DIntegerDatasetStateful(data, dataset_name, nx_ny_nz, xyz_offset, file_id, read_write)
    integer, intent(inout), dimension(:,:,:,:,:,:) :: data
    CHARACTER(len=*), intent(in) :: dataset_name
    INTEGER, intent(in), dimension(3) :: nx_ny_nz, xyz_offset
    INTEGER(HID_T), intent(in) :: file_id
    character(len=*), intent(in) :: read_write
    INTEGER(HID_T) :: filespace, dset_id, memspace, p_access_list
    INTEGER(HSIZE_T) :: array_size(rank(data)), array_subsize(rank(data)), array_start(rank(data))
    integer :: data_rank
    INTEGER(HID_T) :: dtype_id

    call getDataType("integer", dtype_id)
    call findArraySizesStateful(nx_ny_nz, xyz_offset, shape(data),  array_size, array_subsize, array_start)
    data_rank = rank(data)
    call readOrWriteDatasetPreamble(dataset_name, file_id, read_write, dtype_id, data_rank, array_subsize, array_size, array_start, p_access_list, memspace, filespace, dset_id)
    call readOrWriteStateful(array_subsize, dset_id, filespace, memspace, p_access_list, data, read_write)
    call readOrWriteDatasetClosing(dset_id, filespace, memspace, p_access_list)
  end subroutine

  subroutine readWrite6DRealDatasetStateful(data, dataset_name, nx_ny_nz, xyz_offset, file_id, read_write)
    real(k_real), intent(inout), dimension(:,:,:,:,:,:,:) :: data
    CHARACTER(len=*), intent(in) :: dataset_name
    INTEGER, intent(in), dimension(3) :: nx_ny_nz, xyz_offset
    INTEGER(HID_T), intent(in) :: file_id
    character(len=*), intent(in) :: read_write
    INTEGER(HID_T) :: filespace, dset_id, memspace, p_access_list
    INTEGER(HSIZE_T) :: array_size(rank(data)), array_subsize(rank(data)), array_start(rank(data))
    integer :: data_rank
    INTEGER(HID_T) :: dtype_id

    call getDataType("real", dtype_id)
    call findArraySizesStateful(nx_ny_nz, xyz_offset, shape(data),  array_size, array_subsize, array_start)
    data_rank = rank(data)
    call readOrWriteDatasetPreamble(dataset_name, file_id, read_write, dtype_id, data_rank, array_subsize, array_size, array_start, p_access_list, memspace, filespace, dset_id)
    call readOrWriteStateful(array_subsize, dset_id, filespace, memspace, p_access_list, data, read_write)
    call readOrWriteDatasetClosing(dset_id, filespace, memspace, p_access_list)
  end subroutine

  subroutine readWrite6DIntegerDatasetStateful(data, dataset_name, nx_ny_nz, xyz_offset, file_id, read_write)
    integer, intent(inout), dimension(:,:,:,:,:,:,:) :: data
    CHARACTER(len=*), intent(in) :: dataset_name
    INTEGER, intent(in), dimension(3) :: nx_ny_nz, xyz_offset
    INTEGER(HID_T), intent(in) :: file_id
    character(len=*), intent(in) :: read_write
    INTEGER(HID_T) :: filespace, dset_id, memspace, p_access_list
    INTEGER(HSIZE_T) :: array_size(rank(data)), array_subsize(rank(data)), array_start(rank(data))
    integer :: data_rank
    INTEGER(HID_T) :: dtype_id

    call getDataType("integer", dtype_id)
    call findArraySizesStateful(nx_ny_nz, xyz_offset, shape(data),  array_size, array_subsize, array_start)
    data_rank = rank(data)
    call readOrWriteDatasetPreamble(dataset_name, file_id, read_write, dtype_id, data_rank, array_subsize, array_size, array_start, p_access_list, memspace, filespace, dset_id)
    call readOrWriteStateful(array_subsize, dset_id, filespace, memspace, p_access_list, data, read_write)
    call readOrWriteDatasetClosing(dset_id, filespace, memspace, p_access_list)
  end subroutine

  subroutine addRealAttribute(file_id, aname, attr_data)
    INTEGER(HID_T), intent(in) :: file_id
    CHARACTER(LEN=*), intent(in) :: aname
    real(k_real), intent(in) :: attr_data
    integer(HID_T) :: attr_id, &       ! Attribute identifier
                      aspace_id        ! Attribute Dataspace identifier
    integer(HSIZE_T), DIMENSION(1) :: data_dims !unused but needed
    integer :: ierr

    CALL H5Screate_f(H5S_SCALAR_F, aspace_id, ierr)
    if (ierr.ne.0) error stop "addRealAttribute H5Screate_f"
    CALL h5acreate_f(file_id, aname, H5T_NATIVE_DOUBLE, aspace_id, attr_id, ierr)
    if (ierr.ne.0) error stop "addRealAttribute h5acreate_f"
    data_dims(1) = 1
    CALL h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, attr_data, data_dims, ierr)
    if (ierr.ne.0) error stop "addRealAttribute h5awrite_f"
    CALL h5aclose_f(attr_id, ierr)
    if (ierr.ne.0) error stop "addRealAttribute h5aclose_f"
    CALL h5sclose_f(aspace_id, ierr)
    if (ierr.ne.0) error stop "addRealAttribute h5sclose_f"

  endsubroutine

  subroutine readRealAttribute(file_id, aname, attr_data)
    INTEGER(HID_T), intent(in) :: file_id
    CHARACTER(LEN=*), intent(in) :: aname
    real(k_real), intent(out) :: attr_data
    integer(HID_T) :: attr_id, &       ! Attribute identifier
                      atype_id
    integer(HSIZE_T), DIMENSION(1) :: data_dims !unused but needed
    integer :: ierr

    CALL h5aopen_name_f(file_id, aname, attr_id, ierr)
    if (ierr.ne.0) error stop "readRealAttribute h5aopen_name_f"
    CALL h5aget_type_f(attr_id, atype_id, ierr)
    if (ierr.ne.0) error stop "readRealAttribute h5aget_type_f"
    data_dims(1) = 1
    CALL h5aread_f(attr_id, atype_id,  attr_data, data_dims, ierr)
    if (ierr.ne.0) error stop "readRealAttribute h5aread_f"
    CALL h5aclose_f(attr_id, ierr)
    if (ierr.ne.0) error stop "readRealAttribute h5aclose_f"
    CALL h5tclose_f(atype_id, ierr)
    if (ierr.ne.0) error stop "readRealAttribute h5aclose_f"

  endsubroutine

  subroutine addRealMatrixAttribute(file_id, aname, attr_data)
    INTEGER(HID_T), intent(in) :: file_id
    CHARACTER(LEN=*), intent(in) :: aname
    real(k_real), intent(in), DIMENSION(:,:) :: attr_data
    integer(HID_T) :: attr_id, &       ! Attribute identifier
                      aspace_id        ! Attribute Dataspace identifier
    integer(HSIZE_T), DIMENSION(2) :: data_dims !unused but needed
    integer :: ierr

    data_dims = shape(attr_data)
    CALL h5screate_simple_f(2, data_dims, aspace_id, ierr)
    if (ierr.ne.0) error stop "addRealMatrixAttribute h5screate_simple_f"
    CALL h5acreate_f(file_id, aname, H5T_NATIVE_DOUBLE, aspace_id, attr_id, ierr)
    if (ierr.ne.0) error stop "addRealMatrixAttribute h5acreate_f"
    CALL h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, attr_data, data_dims, ierr)
    if (ierr.ne.0) error stop "addRealMatrixAttribute h5awrite_f"
    CALL h5aclose_f(attr_id, ierr)
    if (ierr.ne.0) error stop "addRealMatrixAttribute h5aclose_f"
    CALL h5sclose_f(aspace_id, ierr)
    if (ierr.ne.0) error stop "addRealMatrixAttribute h5sclose_f"

  endsubroutine

  subroutine readRealMatrixAttribute(file_id, aname, attr_data)
    INTEGER(HID_T), intent(in) :: file_id
    CHARACTER(LEN=*), intent(in) :: aname
    real(k_real), dimension(:,:), intent(out) :: attr_data
    integer(HID_T) :: attr_id, &       ! Attribute identifier
                      atype_id
    integer(HSIZE_T), DIMENSION(2) :: data_dims !unused but needed
    integer :: ierr

    CALL h5aopen_name_f(file_id, aname, attr_id, ierr)
    if (ierr.ne.0) error stop "readRealMatrixAttribute h5aopen_name_f"
    CALL h5aget_type_f(attr_id, atype_id, ierr)
    if (ierr.ne.0) error stop "readRealMatrixAttribute h5aget_type_f"
    data_dims = shape(attr_data)
    CALL h5aread_f(attr_id, atype_id,  attr_data, data_dims, ierr)
    if (ierr.ne.0) error stop "readRealMatrixAttribute h5aread_f"
    CALL h5aclose_f(attr_id, ierr)
    if (ierr.ne.0) error stop "readRealMatrixAttribute h5aclose_f"
    CALL h5tclose_f(atype_id, ierr)
    if (ierr.ne.0) error stop "readRealMatrixAttribute h5tclose_f"

  endsubroutine

  subroutine addIntegerAttribute(file_id, aname, attr_data)
    INTEGER(HID_T), intent(in) :: file_id
    CHARACTER(LEN=*), intent(in) :: aname
    integer, intent(in) :: attr_data
    integer(HID_T) :: attr_id, &       ! Attribute identifier
                      aspace_id        ! Attribute Dataspace identifier
    integer(HSIZE_T), DIMENSION(1) :: data_dims !unused but needed
    integer :: ierr


    CALL H5Screate_f(H5S_SCALAR_F, aspace_id, ierr)
    if (ierr.ne.0) error stop "addIntegerAttribute H5Screate_f"
    CALL h5acreate_f(file_id, aname, H5T_NATIVE_INTEGER, aspace_id, attr_id, ierr)
    if (ierr.ne.0) error stop "addIntegerAttribute h5acreate_f"
    data_dims(1) = 1
    CALL h5awrite_f(attr_id, H5T_NATIVE_INTEGER, attr_data, data_dims, ierr)
    if (ierr.ne.0) error stop "addIntegerAttribute h5awrite_f"
    CALL h5aclose_f(attr_id, ierr)
    if (ierr.ne.0) error stop "addIntegerAttribute h5aclose_f"
    CALL h5sclose_f(aspace_id, ierr)
    if (ierr.ne.0) error stop "addIntegerAttribute h5sclose_f"

  endsubroutine

  subroutine readIntegerAttribute(file_id, aname, attr_data)
    INTEGER(HID_T), intent(in) :: file_id
    CHARACTER(LEN=*), intent(in) :: aname
    integer, intent(out) :: attr_data
    integer(HID_T) :: attr_id, &       ! Attribute identifier
                      atype_id
    integer(HSIZE_T), DIMENSION(1) :: data_dims !unused but needed
    integer :: ierr

    CALL h5aopen_name_f(file_id, aname, attr_id, ierr)
    if (ierr.ne.0) error stop "readIntegerAttribute h5aopen_name_f"
    CALL h5aget_type_f(attr_id, atype_id, ierr)
    if (ierr.ne.0) error stop "readIntegerAttribute h5aget_type_f"
    data_dims(1) = 1
    CALL h5aread_f(attr_id, atype_id,  attr_data, data_dims, ierr)
    if (ierr.ne.0) error stop "readIntegerAttribute h5aread_f"
    CALL h5aclose_f(attr_id, ierr)
    if (ierr.ne.0) error stop "readIntegerAttribute h5aclose_f"
    CALL h5tclose_f(atype_id, ierr)
    if (ierr.ne.0) error stop "readIntegerAttribute h5tclose_f"

  endsubroutine

end module
