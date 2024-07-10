module mpi_useful_routines_mod
  use mpi
  use mpi_variables_mod, only : mpi_size, mpi_master_rank, my_mpi_err, i_am_mpi_master
  use kinds
  implicit none
contains

  subroutine writeMPIMaster(string)
    implicit none
    character(len=*), intent(in) :: string
    if (AmIMPIMaster()) write(*,*) string
    call MPI_BARRIER(MPI_COMM_WORLD, my_mpi_err)
  end subroutine

  logical function AmIMPIMaster()
    implicit none
    AmIMPIMaster = .FALSE.
    AmIMPIMaster = i_am_mpi_master
  end function

  subroutine MPIBarrier()
    implicit none
    call MPI_Barrier(MPI_COMM_WORLD, my_mpi_err)
  end subroutine

  subroutine MPIBroadcastScalar(scalar_to_broadcast)
    implicit none
    real(k_real), intent(inout) :: scalar_to_broadcast
    call MPI_Bcast(scalar_to_broadcast, 1, &
                   MPI_DOUBLE_PRECISION, mpi_master_rank, MPI_COMM_WORLD, my_mpi_err)
  end subroutine

  subroutine MPIBroadcastScalarInteger(scalar_to_broadcast)
    implicit none
    integer, intent(inout) :: scalar_to_broadcast
    call MPI_Bcast(scalar_to_broadcast, 1, &
                   MPI_INTEGER, mpi_master_rank, MPI_COMM_WORLD, my_mpi_err)
  end subroutine

  subroutine MPIBroadcastVector(vector_to_broadcast)
    implicit none
    real(k_real), dimension(:), intent(inout) :: vector_to_broadcast
    integer :: n
    n=product(shape(vector_to_broadcast))
    call MPI_Bcast(vector_to_broadcast, n, &
                   MPI_DOUBLE_PRECISION, mpi_master_rank, MPI_COMM_WORLD, my_mpi_err)
  end subroutine

  subroutine MPIBroadcastMatrix(matrix_to_broadcast)
    implicit none
    real(k_real), dimension(:,:), intent(inout) :: matrix_to_broadcast
    integer :: n

    n=product(shape(matrix_to_broadcast))
    call MPI_Bcast(matrix_to_broadcast, n, &
                   MPI_DOUBLE_PRECISION, mpi_master_rank, MPI_COMM_WORLD, my_mpi_err)
  end subroutine

  subroutine MPILogicalOR(logical_value, mpi_logical_or)
    implicit none
    logical, intent(in) :: logical_value
    logical, intent(out) :: mpi_logical_or
    integer :: logical_value_int, mpi_logical_or_int
    ! mpich does not recognizes MPI_LOGICAL, so we need to use integers to make it compatible
    
    logical_value_int = 0
    if (logical_value) logical_value_int =1

    call MPIMAXScalarInteger(logical_value_int, mpi_logical_or_int)

    mpi_logical_or = .FALSE.
    if (mpi_logical_or_int.eq.1) mpi_logical_or = .TRUE.
  end subroutine

  subroutine MPISumScalar(scalar, my_mpi_sum)
    implicit none
    real(k_real), intent(in) :: scalar
    real(k_real), intent(out) :: my_mpi_sum

    call MPI_Reduce(scalar, my_mpi_sum, 1, MPI_DOUBLE_PRECISION, &
                    MPI_SUM, mpi_master_rank, MPI_COMM_WORLD, my_mpi_err)

    call MPI_Bcast(my_mpi_sum, 1, &
                   MPI_DOUBLE_PRECISION, mpi_master_rank, MPI_COMM_WORLD, my_mpi_err)

  end subroutine

  subroutine MPIMAXScalar(scalar, my_mpi_max)
    implicit none
    real(k_real), intent(in) :: scalar
    real(k_real), intent(out) :: my_mpi_max

    call MPI_Reduce(scalar, my_mpi_max, 1, MPI_DOUBLE_PRECISION, &
                    MPI_MAX, mpi_master_rank, MPI_COMM_WORLD, my_mpi_err)

    call MPI_Bcast(my_mpi_max, 1, &
                   MPI_DOUBLE_PRECISION, mpi_master_rank, MPI_COMM_WORLD, my_mpi_err)

  end subroutine

  subroutine MPIMAXScalarInteger(scalar, my_mpi_max)
    implicit none
    integer, intent(in) :: scalar
    integer, intent(out) :: my_mpi_max

    call MPI_Reduce(scalar, my_mpi_max, 1, MPI_INTEGER, &
                    MPI_MAX, mpi_master_rank, MPI_COMM_WORLD, my_mpi_err)

    call MPI_Bcast(my_mpi_max, 1, &
                  MPI_INTEGER, mpi_master_rank, MPI_COMM_WORLD, my_mpi_err)

  end subroutine

  subroutine MPIMaxIncrementGridScalar(grid_scalar, grid_scalar_old, my_mpi_max_increment)
    implicit none
    real(k_real), intent(in), dimension(:,:,:) :: grid_scalar, grid_scalar_old
    real(k_real), intent(out) :: my_mpi_max_increment
    real(k_real) :: max_increment_rank

    max_increment_rank = maxval(abs(grid_scalar-grid_scalar_old))

    call MPI_Reduce(max_increment_rank, my_mpi_max_increment, 1, MPI_DOUBLE_PRECISION, &
                    MPI_MAX, mpi_master_rank, MPI_COMM_WORLD, my_mpi_err)

    call MPI_Bcast(my_mpi_max_increment, 1, &
                   MPI_DOUBLE_PRECISION, mpi_master_rank, MPI_COMM_WORLD, my_mpi_err)

  end subroutine

  subroutine MPIMaxIncrementGridTensor2(grid_tensor2, grid_tensor2_old, my_mpi_max_increment)
    implicit none
    real(k_real), intent(in), dimension(:,:,:,:,:) :: grid_tensor2, grid_tensor2_old
    real(k_real), intent(out) :: my_mpi_max_increment
    real(k_real) :: max_increment_rank

    max_increment_rank = maxval(abs(grid_tensor2-grid_tensor2_old))

    call MPI_Reduce(max_increment_rank, my_mpi_max_increment, 1, MPI_DOUBLE_PRECISION, &
                    MPI_MAX, mpi_master_rank, MPI_COMM_WORLD, my_mpi_err)

    call MPI_Bcast(my_mpi_max_increment, 1, &
                   MPI_DOUBLE_PRECISION, mpi_master_rank, MPI_COMM_WORLD, my_mpi_err)

  end subroutine

  subroutine MPIMaxIncrementGridSSScalar(grid_ss_scalar, grid_ss_scalar_old, my_mpi_max_increment)
    implicit none
    real(k_real), intent(in), dimension(:,:,:,:) :: grid_ss_scalar, grid_ss_scalar_old
    real(k_real), intent(out) :: my_mpi_max_increment
    real(k_real) :: max_increment_rank

    max_increment_rank = maxval(abs(grid_ss_scalar-grid_ss_scalar_old))

    call MPI_Reduce(max_increment_rank, my_mpi_max_increment, 1, MPI_DOUBLE_PRECISION, &
                    MPI_MAX, mpi_master_rank, MPI_COMM_WORLD, my_mpi_err)

    call MPI_Bcast(my_mpi_max_increment, 1, &
                   MPI_DOUBLE_PRECISION, mpi_master_rank, MPI_COMM_WORLD, my_mpi_err)

  end subroutine

  subroutine MPIComputeQualitySSScalar(grid_ss_scalar, grid_ss_scalar_old, grid_ss_scalar_rate, grid_ss_scalar_rate_old, n_points, dt, L, T, Q, my_mpi_max_increment)
    use tensor_math_mod, only : WeightedvectorNorm
    implicit none
    real(k_real), intent(in), dimension(:,:,:,:) :: grid_ss_scalar, grid_ss_scalar_old, grid_ss_scalar_rate, grid_ss_scalar_rate_old
    real(k_real), intent(in) :: dt, L, T
    integer, intent(in) :: n_points
    integer :: n_vals
    real(k_real), intent(out) :: Q, my_mpi_max_increment
    real(k_real) :: q_point_square
    n_vals = size(grid_ss_scalar_old, 1)

    call MPIMaxIncrementGridSSScalar(grid_ss_scalar, grid_ss_scalar_old, my_mpi_max_increment)
    q_point_square = sum( ((grid_ss_scalar_rate-grid_ss_scalar_rate_old)*dt/( L*(grid_ss_scalar - grid_ss_scalar_old) + T))**2 )

    call MPI_Reduce(q_point_square, Q, 1, MPI_DOUBLE_PRECISION, &
                    MPI_SUM, mpi_master_rank, MPI_COMM_WORLD, my_mpi_err)
    Q = sqrt(Q/int2real(n_points))

    call MPI_Bcast(Q, 1, &
                   MPI_DOUBLE_PRECISION, mpi_master_rank, MPI_COMM_WORLD, my_mpi_err)
  end subroutine

  subroutine MPIComputeQualityTensor2(grid_tensor2, grid_tensor2_old, grid_tensor2_rate, grid_tensor2_rate_old, n_points, dt, L, T, Q, my_mpi_max_increment)
    use tensor_math_mod, only : Weightedtensor2Norm
    implicit none
    real(k_real), intent(in), dimension(:,:,:,:,:) :: grid_tensor2, grid_tensor2_old, grid_tensor2_rate, grid_tensor2_rate_old
    real(k_real), intent(in) :: dt, L, T
    integer, intent(in) :: n_points
    real(k_real), intent(out) :: Q, my_mpi_max_increment
    real(k_real) :: q_point_square

    call MPIMaxIncrementGridTensor2(grid_tensor2, grid_tensor2_old, my_mpi_max_increment)
    q_point_square = sum( ((grid_tensor2_rate-grid_tensor2_rate_old)*dt/( L*(grid_tensor2 - grid_tensor2_old) + T))**2 )

    call MPI_Reduce(q_point_square, Q, 1, MPI_DOUBLE_PRECISION, &
                    MPI_SUM, mpi_master_rank, MPI_COMM_WORLD, my_mpi_err)
    Q = sqrt(Q/int2real(n_points))
    call MPI_Bcast(Q, 1, &
                   MPI_DOUBLE_PRECISION, mpi_master_rank, MPI_COMM_WORLD, my_mpi_err)

  end subroutine

  subroutine MPIComputeQualityScalar(grid_scalar, grid_scalar_old, grid_scalar_rate, grid_scalar_rate_old, n_points, dt, L, T, Q, my_mpi_max_increment)
    implicit none
    real(k_real), intent(in), dimension(:,:,:) :: grid_scalar, grid_scalar_old, grid_scalar_rate, grid_scalar_rate_old
    real(k_real), intent(in) :: dt, L, T
    integer, intent(in) :: n_points
    real(k_real), intent(out) :: Q, my_mpi_max_increment
    real(k_real) :: q_point_square

    call MPIMaxIncrementGridScalar(grid_scalar, grid_scalar_old, my_mpi_max_increment)
    q_point_square = sum( ((grid_scalar_rate-grid_scalar_rate_old)*dt/( L*(grid_scalar - grid_scalar_old) + T))**2 )

    call MPI_Reduce(q_point_square, Q, 1, MPI_DOUBLE_PRECISION, &
                    MPI_SUM, mpi_master_rank, MPI_COMM_WORLD, my_mpi_err)
    Q = sqrt(Q/int2real(n_points))
    
    call MPI_Bcast(Q, 1, &
                   MPI_DOUBLE_PRECISION, mpi_master_rank, MPI_COMM_WORLD, my_mpi_err)

  end subroutine

  subroutine MPIMaxRelativeIncrementGridSSScalar(grid_ss_scalar, grid_ss_scalar_old, my_mpi_max_relative_increment)
    implicit none
    real(k_real), intent(in), dimension(:,:,:,:) :: grid_ss_scalar, grid_ss_scalar_old
    real(k_real), intent(out) :: my_mpi_max_relative_increment
    real(k_real) :: max_relative_increment_rank

    max_relative_increment_rank = maxval(abs(grid_ss_scalar-grid_ss_scalar_old)/abs(grid_ss_scalar_old))

    call MPI_Reduce(max_relative_increment_rank, my_mpi_max_relative_increment, 1, MPI_DOUBLE_PRECISION, &
                    MPI_MAX, mpi_master_rank, MPI_COMM_WORLD, my_mpi_err)

    call MPI_Bcast(my_mpi_max_relative_increment, 1, &
                   MPI_DOUBLE_PRECISION, mpi_master_rank, MPI_COMM_WORLD, my_mpi_err)

  end subroutine

  subroutine MPIMaxGridScalar(grid_scalar, my_mpi_max)
    implicit none
    real(k_real), intent(in), dimension(:,:,:) :: grid_scalar
    real(k_real), intent(out) :: my_mpi_max
    real(k_real) :: max_val_rank

    max_val_rank = maxval(grid_scalar)

    call MPI_Reduce(max_val_rank, my_mpi_max, 1, MPI_DOUBLE_PRECISION, &
                    MPI_MAX, mpi_master_rank, MPI_COMM_WORLD, my_mpi_err)

    call MPI_Bcast(my_mpi_max, 1, &
                   MPI_DOUBLE_PRECISION, mpi_master_rank, MPI_COMM_WORLD, my_mpi_err)

  end subroutine

  subroutine MPICountTrueInMask(grid_mask, mpi_num_true)
    implicit none
    logical, intent(in), dimension(:,:,:) :: grid_mask
    integer, intent(out) :: mpi_num_true
    integer :: mpi_num_true_rank

    
    mpi_num_true_rank = count(grid_mask)

    call MPI_Reduce(mpi_num_true_rank, mpi_num_true, 1, MPI_INTEGER, &
                    MPI_SUM, mpi_master_rank, MPI_COMM_WORLD, my_mpi_err)

    call MPI_Bcast(mpi_num_true, 1, &
                   MPI_INTEGER, mpi_master_rank, MPI_COMM_WORLD, my_mpi_err)

  end subroutine

  subroutine MPIMaxGridSSScalar(grid_ss_scalar, my_mpi_max)
    implicit none
    real(k_real), intent(in), dimension(:,:,:,:) :: grid_ss_scalar
    real(k_real), intent(out) :: my_mpi_max
    real(k_real) :: max_val_rank

    max_val_rank = maxval(grid_ss_scalar)

    call MPI_Reduce(max_val_rank, my_mpi_max, 1, MPI_DOUBLE_PRECISION, &
                    MPI_MAX, mpi_master_rank, MPI_COMM_WORLD, my_mpi_err)

    call MPI_Bcast(my_mpi_max, 1, &
                   MPI_DOUBLE_PRECISION, mpi_master_rank, MPI_COMM_WORLD, my_mpi_err)

  end subroutine

  subroutine MPISumScalarInteger(scalar, my_mpi_sum)
    implicit none
    integer, intent(in) :: scalar
    integer, intent(out) :: my_mpi_sum

    call MPI_Reduce(scalar, my_mpi_sum, 1, MPI_INTEGER, &
                    MPI_SUM, mpi_master_rank, MPI_COMM_WORLD, my_mpi_err)

    call MPI_Bcast(my_mpi_sum, 1, &
                   MPI_INTEGER, mpi_master_rank, MPI_COMM_WORLD, my_mpi_err)

  end subroutine

  subroutine MPISumMatrix(matrix, my_mpi_sum)
    implicit none
    real(k_real), intent(in) :: matrix(:,:)
    integer :: n_vals
    real(k_real), dimension(:,:), intent(out) :: my_mpi_sum

    n_vals = product(SHAPE(matrix))

    call MPI_Reduce(matrix, my_mpi_sum, n_vals, MPI_DOUBLE_PRECISION, &
                    MPI_SUM, mpi_master_rank, MPI_COMM_WORLD, my_mpi_err)

    call MPI_Bcast(my_mpi_sum, n_vals, &
                   MPI_DOUBLE_PRECISION, mpi_master_rank, MPI_COMM_WORLD, my_mpi_err)

  end subroutine

  subroutine MPIAverageGridScalar(scalar, mpi_average)
    use mpi_variables_mod, only : mpi_local_z_all_proc
    implicit none
    real(k_real), intent(in) :: scalar(:,:,:)
    integer :: array_size(3), n_vals, n_points, ngrid_points
    real(k_real), intent(out) :: mpi_average
    real(k_real) :: proc_total
    integer :: x, y, z

    array_size = SHAPE(scalar)
    n_vals = 1
    n_points = product(array_size(1:3))
    ngrid_points = product(array_size(1:2))*sum(mpi_local_z_all_proc)

    proc_total = 0._k_real
    do z=1,array_size(3)
    do y=1,array_size(2)
    do x=1,array_size(1)
      proc_total = proc_total + scalar(x,y,z)/int2real(ngrid_points)
    enddo
    enddo
    enddo


    call MPI_Reduce(proc_total, mpi_average, n_vals, MPI_DOUBLE_PRECISION, &
                    MPI_SUM, mpi_master_rank, MPI_COMM_WORLD, my_mpi_err)

    call MPI_Bcast(mpi_average, n_vals, &
                   MPI_DOUBLE_PRECISION, mpi_master_rank, MPI_COMM_WORLD, my_mpi_err)

  end subroutine

  subroutine MPIAverageGridScalarMasked(scalar, grid_mask, mpi_average)
    implicit none
    real(k_real), intent(in) :: scalar(:,:,:)
    logical, intent(in) :: grid_mask(:,:,:)
    integer :: array_size(3), n_vals
    real(k_real), intent(out) :: mpi_average
    real(k_real) :: proc_total
    integer :: x, y, z, num_true

    call MPICountTrueInMask(grid_mask, num_true)
    array_size = SHAPE(scalar)
    n_vals = 1

    proc_total = 0._k_real
    do z=1,array_size(3)
    do y=1,array_size(2)
    do x=1,array_size(1)
      if (grid_mask(x,y,z)) &
      proc_total = proc_total + scalar(x,y,z)/int2real(num_true)
    enddo
    enddo
    enddo


    call MPI_Reduce(proc_total, mpi_average, n_vals, MPI_DOUBLE_PRECISION, &
                    MPI_SUM, mpi_master_rank, MPI_COMM_WORLD, my_mpi_err)

    call MPI_Bcast(mpi_average, n_vals, &
                   MPI_DOUBLE_PRECISION, mpi_master_rank, MPI_COMM_WORLD, my_mpi_err)

  end subroutine

  subroutine MPIMaximumGridScalar(scalar, global_max)
    implicit none
    real(k_real), intent(in) :: scalar(:,:,:)
    integer :: n_vals
    real(k_real), intent(out) :: global_max
    real(k_real) :: proc_max

    n_vals = 1
    proc_max = maxval(scalar)

    call MPI_Reduce(proc_max, global_max, n_vals, MPI_DOUBLE_PRECISION, &
                    MPI_MAX, mpi_master_rank, MPI_COMM_WORLD, my_mpi_err)

    call MPI_Bcast(global_max, n_vals, &
                   MPI_DOUBLE_PRECISION, mpi_master_rank, MPI_COMM_WORLD, my_mpi_err)

  end subroutine


  subroutine MPIMaximumGridScalarMasked(scalar, grid_mask, global_max)
    implicit none
    real(k_real), intent(in) :: scalar(:,:,:)
    logical, intent(in) :: grid_mask(:,:,:)
    integer :: n_vals
    real(k_real), intent(out) :: global_max
    real(k_real) :: proc_max

    n_vals = 1
    proc_max = maxval(scalar, grid_mask)

    call MPI_Reduce(proc_max, global_max, n_vals, MPI_DOUBLE_PRECISION, &
                    MPI_MAX, mpi_master_rank, MPI_COMM_WORLD, my_mpi_err)

    call MPI_Bcast(global_max, n_vals, &
                   MPI_DOUBLE_PRECISION, mpi_master_rank, MPI_COMM_WORLD, my_mpi_err)

  end subroutine

  subroutine MPIAverageGridVector(grid_vector, mpi_average)
    use mpi_variables_mod, only : mpi_local_z_all_proc
    implicit none
    real(k_real), intent(in) :: grid_vector(:,:,:,:)
    integer :: array_size(4), n_vals, ngrid_points
    real(k_real), dimension(:), intent(out) :: mpi_average
    real(k_real) :: proc_total(size(grid_vector,1))
    integer :: x, y, z

    array_size =  SHAPE(grid_vector)
    n_vals = array_size(1)
    ngrid_points = product(array_size(2:3))*sum(mpi_local_z_all_proc)

    proc_total = 0._k_real
    do z=1,array_size(4)
    do y=1,array_size(3)
    do x=1,array_size(2)
      proc_total = proc_total + grid_vector(:,x,y,z)/int2real(ngrid_points)
    enddo
    enddo
    enddo


    call MPI_Reduce(proc_total, mpi_average, n_vals, MPI_DOUBLE_PRECISION, &
                    MPI_SUM, mpi_master_rank, MPI_COMM_WORLD, my_mpi_err)

    call MPI_Bcast(mpi_average, n_vals, &
                   MPI_DOUBLE_PRECISION, mpi_master_rank, MPI_COMM_WORLD, my_mpi_err)

  end subroutine

  subroutine MPIAverageVoxelWeightGridVector(grid_vector, voxel_weight, mpi_average)
    use mpi_variables_mod, only : mpi_local_z_all_proc
    implicit none
    real(k_real), intent(in) :: grid_vector(:,:,:,:), voxel_weight(:,:,:)
    integer :: array_size(4), n_vals, ngrid_points
    real(k_real), dimension(:), intent(out) :: mpi_average
    real(k_real) :: proc_total(size(grid_vector,1))
    integer :: x, y, z

    array_size =  SHAPE(grid_vector)
    n_vals = array_size(1)
    ngrid_points = product(array_size(2:3))*sum(mpi_local_z_all_proc)

    proc_total = 0._k_real
    do z=1,array_size(4)
    do y=1,array_size(3)
    do x=1,array_size(2)
      proc_total = proc_total + grid_vector(:,x,y,z)/int2real(ngrid_points) * voxel_weight(x,y,z)
    enddo
    enddo
    enddo


    call MPI_Reduce(proc_total, mpi_average, n_vals, MPI_DOUBLE_PRECISION, &
                    MPI_SUM, mpi_master_rank, MPI_COMM_WORLD, my_mpi_err)

    call MPI_Bcast(mpi_average, n_vals, &
                   MPI_DOUBLE_PRECISION, mpi_master_rank, MPI_COMM_WORLD, my_mpi_err)

  end subroutine

  subroutine MPIAverageGridVectorMasked(grid_vector, grid_mask, mpi_average)
    implicit none
    real(k_real), intent(in) :: grid_vector(:,:,:,:)
    logical, intent(in) :: grid_mask(:,:,:)
    integer :: array_size(4), n_vals, num_true
    real(k_real), dimension(:), intent(out) :: mpi_average
    real(k_real) :: proc_total(size(grid_vector,1))
    integer :: x, y, z

    array_size =  SHAPE(grid_vector)
    n_vals = array_size(1)
    call MPICountTrueInMask(grid_mask, num_true)

    proc_total = 0._k_real
    do z=1,array_size(4)
    do y=1,array_size(3)
    do x=1,array_size(2)
      if (grid_mask(x,y,z)) &
      proc_total = proc_total + grid_vector(:,x,y,z) / int2real(num_true)
    enddo
    enddo
    enddo


    call MPI_Reduce(proc_total, mpi_average, n_vals, MPI_DOUBLE_PRECISION, &
                    MPI_SUM, mpi_master_rank, MPI_COMM_WORLD, my_mpi_err)

    call MPI_Bcast(mpi_average, n_vals, &
                   MPI_DOUBLE_PRECISION, mpi_master_rank, MPI_COMM_WORLD, my_mpi_err)

  end subroutine

  subroutine MPIAverageGridMatrix(matrix, mpi_average)
    use mpi_variables_mod, only : mpi_local_z_all_proc
    implicit none
    real(k_real), intent(in) :: matrix(:,:,:,:,:)
    integer :: array_size(5), n_vals, ngrid_points
    real(k_real), dimension(:,:), intent(out) :: mpi_average
    real(k_real) :: proc_total(size(matrix,1), size(matrix,2))
    integer :: x, y, z

    array_size = SHAPE(matrix)
    n_vals = product(array_size(1:2))
    ngrid_points = product(array_size(3:4))*sum(mpi_local_z_all_proc)

    proc_total = 0._k_real
    do z=1,array_size(5)
    do y=1,array_size(4)
    do x=1,array_size(3)
      proc_total(:,:) = proc_total(:,:) + matrix(:,:,x,y,z)/int2real(ngrid_points)
    enddo
    enddo
    enddo


    call MPI_Reduce(proc_total, mpi_average, n_vals, MPI_DOUBLE_PRECISION, &
                    MPI_SUM, mpi_master_rank, MPI_COMM_WORLD, my_mpi_err)

    call MPI_Bcast(mpi_average, n_vals, &
                   MPI_DOUBLE_PRECISION, mpi_master_rank, MPI_COMM_WORLD, my_mpi_err)

  end subroutine

  subroutine MPIAverageVoxelWeightGridMatrix(matrix, voxel_weight, mpi_average)
    use mpi_variables_mod, only : mpi_local_z_all_proc
    implicit none
    real(k_real), intent(in) :: matrix(:,:,:,:,:), voxel_weight(:,:,:)
    integer :: array_size(5), n_vals, ngrid_points
    real(k_real), dimension(:,:), intent(out) :: mpi_average
    real(k_real) :: proc_total(size(matrix,1), size(matrix,2))
    integer :: x, y, z

    array_size = SHAPE(matrix)
    n_vals = product(array_size(1:2))
    ngrid_points = product(array_size(3:4))*sum(mpi_local_z_all_proc)

    proc_total = 0._k_real
    do z=1,array_size(5)
    do y=1,array_size(4)
    do x=1,array_size(3)
      proc_total(:,:) = proc_total(:,:) + matrix(:,:,x,y,z)/int2real(ngrid_points)*voxel_weight(x,y,z)
    enddo
    enddo
    enddo


    call MPI_Reduce(proc_total, mpi_average, n_vals, MPI_DOUBLE_PRECISION, &
                    MPI_SUM, mpi_master_rank, MPI_COMM_WORLD, my_mpi_err)

    call MPI_Bcast(mpi_average, n_vals, &
                   MPI_DOUBLE_PRECISION, mpi_master_rank, MPI_COMM_WORLD, my_mpi_err)

  end subroutine

  subroutine MPIAverageGridMatrixMasked(matrix, grid_mask, mpi_average)
    implicit none
    real(k_real), intent(in) :: matrix(:,:,:,:,:)
    logical, intent(in) :: grid_mask(:,:,:)
    integer :: array_size(5), n_vals, num_true
    real(k_real), dimension(:,:), intent(out) :: mpi_average
    real(k_real) :: proc_total(size(matrix,1), size(matrix,2))
    integer :: x, y, z

    array_size = SHAPE(matrix)
    n_vals = product(array_size(1:2))
    call MPICountTrueInMask(grid_mask, num_true)

    proc_total = 0._k_real
    do z=1,array_size(5)
    do y=1,array_size(4)
    do x=1,array_size(3)
      if (grid_mask(x,y,z)) &
      proc_total(:,:) = proc_total(:,:) + matrix(:,:,x,y,z)/int2real(num_true)
    enddo
    enddo
    enddo
    

    call MPI_Reduce(proc_total, mpi_average, n_vals, MPI_DOUBLE_PRECISION, &
                    MPI_SUM, mpi_master_rank, MPI_COMM_WORLD, my_mpi_err)

    call MPI_Bcast(mpi_average, n_vals, &
                   MPI_DOUBLE_PRECISION, mpi_master_rank, MPI_COMM_WORLD, my_mpi_err)

  end subroutine

  subroutine MPIAverageGridSSMatrix(ssmatrix, mpi_average)
    use mpi_variables_mod, only : mpi_local_z_all_proc
    implicit none
    real(k_real), intent(in) :: ssmatrix(:,:,:,:,:,:)
    integer :: array_size(6), n_vals, n_points, ngrid_points
    real(k_real), dimension(:,:,:), intent(out) :: mpi_average
    real(k_real) :: proc_total(size(ssmatrix,1), size(ssmatrix,2), size(ssmatrix,3))
    integer :: x, y, z

    array_size = SHAPE(ssmatrix)
    n_vals = product(array_size(1:3))
    n_points = product(array_size(4:6))
    ngrid_points = product(array_size(4:5))*sum(mpi_local_z_all_proc)

    proc_total = 0._k_real
    do z=1,array_size(6)
    do y=1,array_size(5)
    do x=1,array_size(4)
      proc_total = proc_total + ssmatrix(:,:,:,x,y,z)/int2real(ngrid_points)
    enddo
    enddo
    enddo


    call MPI_Reduce(proc_total, mpi_average, n_vals, MPI_DOUBLE_PRECISION, &
                    MPI_SUM, mpi_master_rank, MPI_COMM_WORLD, my_mpi_err)

    call MPI_Bcast(mpi_average, n_vals, &
                   MPI_DOUBLE_PRECISION, mpi_master_rank, MPI_COMM_WORLD, my_mpi_err)

  end subroutine

  subroutine MPIGatherVGridScalar(gridscalar, gridscalar_full)
    use mpi_variables_mod, only : mpi_local_z_all_proc
    implicit none
    real(k_real), intent(in) :: gridscalar(:,:,:)
    real(k_real), intent(out) :: gridscalar_full(:,:,:)
    integer :: npoints_proc, npoints_all(mpi_size), displacement_all(mpi_size), &
               i_rank, &
               nx, ny, nz, nvals

    nvals = 1
    nx=size(gridscalar,rank(gridscalar)-2)
    ny=size(gridscalar,rank(gridscalar)-1)
    nz=size(gridscalar,rank(gridscalar))
    npoints_proc = nvals*nx*ny*nz
    npoints_all = nvals*nx*ny*mpi_local_z_all_proc

    displacement_all(1) = 0
    do i_rank = 2,mpi_size
      displacement_all(i_rank) = displacement_all(i_rank-1) + npoints_all(i_rank-1)
    enddo

    call MPI_Gatherv(gridscalar, npoints_proc, MPI_DOUBLE_PRECISION, &
                     gridscalar_full, npoints_all, displacement_all, &
                     MPI_DOUBLE_PRECISION, mpi_master_rank, MPI_COMM_WORLD, my_mpi_err)

  end subroutine

  subroutine MPIGatherVGridScalarInteger(gridscalar, gridscalar_full)
    use mpi_variables_mod, only : mpi_local_z_all_proc
    implicit none
    integer, intent(in) :: gridscalar(:,:,:)
    integer, intent(out) :: gridscalar_full(:,:,:)
    integer :: npoints_proc, npoints_all(mpi_size), displacement_all(mpi_size), &
               i_rank, &
               nx, ny, nz, nvals

    nvals = 1
    nx=size(gridscalar,rank(gridscalar)-2)
    ny=size(gridscalar,rank(gridscalar)-1)
    nz=size(gridscalar,rank(gridscalar))
    npoints_proc = nvals*nx*ny*nz
    npoints_all = nvals*nx*ny*mpi_local_z_all_proc

    displacement_all(1) = 0
    do i_rank = 2,mpi_size
      displacement_all(i_rank) = displacement_all(i_rank-1) + npoints_all(i_rank-1)
    enddo

    call MPI_Gatherv(gridscalar, npoints_proc, MPI_INTEGER, &
                     gridscalar_full, npoints_all, displacement_all, &
                     MPI_INTEGER, mpi_master_rank, MPI_COMM_WORLD, my_mpi_err)

  end subroutine

  subroutine MPIGatherVGridTensor2(gridtensor, gridtensor_full)
    use mpi_variables_mod, only : mpi_local_z_all_proc
    implicit none
    real(k_real), intent(in) :: gridtensor(:,:,:,:,:)
    real(k_real), intent(out) :: gridtensor_full(:,:,:,:,:)
    integer :: npoints_proc, npoints_all(mpi_size), displacement_all(mpi_size), &
               i_rank, &
               nx, ny, nz, nvals

    nvals = 9
    nx=size(gridtensor,rank(gridtensor)-2)
    ny=size(gridtensor,rank(gridtensor)-1)
    nz=size(gridtensor,rank(gridtensor))
    npoints_proc = nvals*nx*ny*nz
    npoints_all = nvals*nx*ny*mpi_local_z_all_proc

    displacement_all(1) = 0
    do i_rank = 2,mpi_size
      displacement_all(i_rank) = displacement_all(i_rank-1) + npoints_all(i_rank-1)
    enddo

    call MPI_Gatherv(gridtensor, npoints_proc, MPI_DOUBLE_PRECISION, &
                     gridtensor_full, npoints_all, displacement_all, &
                     MPI_DOUBLE_PRECISION, mpi_master_rank, MPI_COMM_WORLD, my_mpi_err)

  end subroutine

  subroutine MPIGatherVGridRealVector(gridvector, gridvector_full)
    use mpi_variables_mod, only : mpi_local_z_all_proc
    implicit none
    real(k_real), intent(in) :: gridvector(:,:,:,:)
    real(k_real), intent(out) :: gridvector_full(:,:,:,:)
    integer :: npoints_proc, npoints_all(mpi_size), displacement_all(mpi_size), &
               i_rank, &
               nx, ny, nz, nvals

    nvals = 3
    nx=size(gridvector,rank(gridvector)-2)
    ny=size(gridvector,rank(gridvector)-1)
    nz=size(gridvector,rank(gridvector))
    npoints_proc = nvals*nx*ny*nz
    npoints_all = nvals*nx*ny*mpi_local_z_all_proc

    displacement_all(1) = 0
    do i_rank = 2,mpi_size
      displacement_all(i_rank) = displacement_all(i_rank-1) + npoints_all(i_rank-1)
    enddo

    call MPI_Gatherv(gridvector, npoints_proc, MPI_DOUBLE_PRECISION, &
                     gridvector_full, npoints_all, displacement_all, &
                     MPI_DOUBLE_PRECISION, mpi_master_rank, MPI_COMM_WORLD, my_mpi_err)

  end subroutine
end module
