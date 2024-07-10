module mpi_variables_mod
  use mpi
  implicit none
  ! MPI VARIABLES
  integer mpi_rank, mpi_size, my_mpi_err, my_mpi_tag, my_mpi_status(MPI_STATUS_SIZE)
  integer, parameter :: mpi_master_rank = 0
  logical :: i_am_mpi_master
  integer, allocatable, dimension(:) :: mpi_local_z_all_proc
  integer :: mpi_integer_buffer
  contains
    
subroutine init_mpi()
  implicit none
  call MPI_INIT(my_mpi_err)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, mpi_size, my_mpi_err)
  call MPI_COMM_RANK(MPI_COMM_WORLD, mpi_rank, my_mpi_err)
  allocate(mpi_local_z_all_proc(mpi_size))
  print *, 'Hello World from process: ', mpi_rank, 'of ', mpi_size
  i_am_mpi_master = .false.
  if (mpi_rank.eq.mpi_master_rank) i_am_mpi_master=.true.
end subroutine

subroutine finalize_mpi()
  implicit none
  call MPI_FINALIZE(my_mpi_err)
end subroutine
end module
