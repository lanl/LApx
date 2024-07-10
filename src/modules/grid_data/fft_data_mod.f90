module fft_grid_data_mod
USE, intrinsic :: iso_c_binding
use kinds
use mpi, only : MPI_COMM_WORLD, MPI_INTEGER
use mpi_variables_mod, only : mpi_rank, my_mpi_err
use all_grid_data_mod, only : all_grid_data
implicit none
include 'fftw3-mpi.f03'
! integer*8 :: fftw_fwd, fftw_bkwd, fftw_est, fftw_norm_order
integer, parameter :: fftw_fwd=-1 , fftw_bkwd=1 , fftw_est=0 ,fftw_norm_order=0
! integer*8 :: FFTW_ESTIMATE , FFTW_MEASURE
! integer , parameter :: FFTW_ESTIMATE = 0 , FFTW_MEASURE = 1


! integer*8 :: FFTW_OUT_OF_PLACE , FFTW_IN_PLACE , FFTW_USE_WISDOM
integer , parameter :: FFTW_OUT_OF_PLACE = 0
integer, parameter :: FFTW_IN_PLACE = 8 , FFTW_USE_WISDOM = 16

real(k_real), pointer :: fft_frequency_vector_k(:,:,:,:) => null(), &
                         fft_frequency_vector_norm(:,:,:) => null(), &  
                         fft_frequency_tensor_K(:,:,:,:,:) => null()     

type fft_data_container
  complex(C_DOUBLE_COMPLEX), pointer :: data(:,:,:) => null()
  type(C_PTR) :: cdata
end type

type fft_grid_data_type
  type(fft_data_container), dimension(:), pointer :: data_container => null()
  complex(C_DOUBLE_COMPLEX), pointer :: data(:,:,:) => null()
  logical :: initialized = .false.
  integer(8) :: nx=-1, ny=-1, nz=-1, nx_local=-1, x_offset=-1, ny_local=-1, y_offset=-1, nz_local=-1, z_offset=-1
  type(C_PTR) :: plan(9), iplan(9)
  integer :: scalar_index = 1
  class(all_grid_data), pointer :: grid_data => null()

  real(k_real), pointer, dimension(:,:,:,:) :: frequency_vector => null()
  real(k_real), pointer, dimension(:,:,:,:,:) :: frequency_tensor => null()
  real(k_real), pointer, dimension(:,:,:) :: frequency_vector_norm => null()

contains
  procedure initFFTData
  procedure :: getFFTDataAtPointRealTensor2
  procedure :: getFFTDataAtPointComplexTensor2
  procedure :: getFFTDataAtPointRealScalar
  procedure :: getFFTDataAtPointComplexScalar
  procedure, private :: setFFTDataAtPointFromComplexTensor2
  procedure, private :: setFFTDataAtPointFromRealTensor2
  procedure, private :: setFFTDataAtPointFromComplexScalar
  procedure, private :: setFFTDataAtPointFromRealScalar
  generic :: setFFTDataAtPointTensor2 => setFFTDataAtPointFromRealTensor2, setFFTDataAtPointFromComplexTensor2
  generic :: setFFTDataAtPointScalar => setFFTDataAtPointFromRealScalar, setFFTDataAtPointFromComplexScalar
  procedure :: executeDFTTensor2
  procedure :: executeDFTScalar
  procedure :: executeInverseDFTTensor2
  procedure :: executeInverseDFTScalar
  procedure :: destroyPlans
  procedure :: addFieldVariablesToGrid
  procedure :: initGridPointers
  procedure :: initFrequencyVectorAndTensor

end type fft_grid_data_type
  
contains

  subroutine addFieldVariablesToGrid(this, grid_data)
    use grid_data_var_type_mod
    implicit none
    class(fft_grid_data_type), intent(inout) :: this
    class(all_grid_data), intent(in), target :: grid_data
    this%grid_data => grid_data 

    associate (all_grid_data_vars => this%grid_data)

      call all_grid_data_vars%addVar("frequency_vector_k", real_vector)
      call all_grid_data_vars%addVar("frequency_vector_norm", scalar)
      call all_grid_data_vars%addVar("frequency_tensor_K", tensor2)
      
    end associate

  end subroutine

  subroutine initGridPointers(this)
    implicit none
    class(fft_grid_data_type), intent(inout) :: this

    call this%grid_data%getRealVectorDataPointerByName("frequency_vector_k", this%frequency_vector)
    call this%grid_data%getScalarDataPointerByName("frequency_vector_norm", this%frequency_vector_norm)
    call this%grid_data%getTensor2DataPointerByName("frequency_tensor_K", this%frequency_tensor)

    call this%grid_data%getRealVectorDataPointerByName("frequency_vector_k", fft_frequency_vector_k)
    call this%grid_data%getScalarDataPointerByName("frequency_vector_norm", fft_frequency_vector_norm)
    call this%grid_data%getTensor2DataPointerByName("frequency_tensor_K", fft_frequency_tensor_K)

  end subroutine

  subroutine initFrequencyVectorAndTensor(this, voxel_size)
    ! use green_function_mod
    use math_constants, only : PI
    use tensor_math_mod, only : vector3OuterProduct, vectorNorm
    implicit none
    class(fft_grid_data_type), intent(inout) :: this
    real(k_real), intent(in) :: voxel_size(3)
    integer(8) :: ix, iy, iz, kx, ky, kz, &
               ix_global, iy_global, iz_global

    do iz=1,this%nz_local
      iz_global = iz+this%z_offset
      do iy=1,this%ny_local
        iy_global = iy+this%y_offset
        do ix=1,this%nx_local
          ix_global = ix+this%x_offset

          associate( nx => this%nx,  ny => this%ny, nz => this%nz, &
                     fv=>this%frequency_vector(:,ix,iy,iz), fv_norm=>this%frequency_vector_norm(ix,iy,iz), &
                     fT2=>this%frequency_tensor(:,:,ix,iy,iz))
            if(ix_global.le.nx/2) kx=ix_global-1
            if(ix_global.gt.nx/2) kx=ix_global-nx-1
        
            if(iy_global.le.ny/2) ky=iy_global-1
            if(iy_global.gt.ny/2) ky=iy_global-ny-1
        
            if(iz_global.le.nz/2) kz=iz_global-1
            if(iz_global.gt.nz/2) kz=iz_global-nz-1
            
            fv(1)=2._k_real*PI*kx/(voxel_size(1)*nx)
            fv(2)=2._k_real*PI*ky/(voxel_size(2)*ny)
            fv(3)=2._k_real*PI*kz/(voxel_size(3)*nz)
            fv_norm=vectorNorm(fv)
            if (fv_norm.ne.0._k_real) fv=fv/fv_norm

            fT2 = vector3OuterProduct(fv, fv)
          end associate 

        enddo
      enddo
    enddo

  end subroutine



  subroutine initFFTData(this, nx,ny,nz, nz_local, z_offset, nz_start, nz_end)
    use mpi_variables_mod, only : mpi_size, finalize_mpi, mpi_local_z_all_proc, mpi_master_rank
    use mpi_useful_routines_mod, only : MPILogicalOR
    implicit none
    class(fft_grid_data_type), intent(inout) :: this
    integer, intent(in) :: nx, ny, nz
    integer, intent(out):: nz_local, z_offset, nz_start, nz_end
    integer(C_INTPTR_T) :: nz_local_cptr, z_offset_cptr
    integer(C_INTPTR_T) :: alloc_local
    logical :: error_on_slice_distribution_proc = .false., &
               error_on_slice_distribution
    integer :: i
    if (nx<1) error stop " initFFTData nx <1"
    if (ny<1) error stop " initFFTData ny <1"
    if (nz<1) error stop " initFFTData nz <1"

    this%nx = nx
    this%ny = ny
    this%nz = nz
    call fftw_mpi_init()

    allocate(this%data_container(9))
    alloc_local = fftw_mpi_local_size_3d(this%nz, this%ny, this%nx, MPI_COMM_WORLD, nz_local_cptr, z_offset_cptr)
    z_offset = int(z_offset_cptr)
    nz_local = int(nz_local_cptr)
    this%nx_local = nx
    this%x_offset = 0
    this%ny_local = ny
    this%y_offset = 0
    this%nz_local = nz_local
    this%z_offset = z_offset
    nz_start = 1 + z_offset
    nz_end = nz_start + nz_local - 1

    error_on_slice_distribution_proc = .false.
    IF (nz_local.lt.1) then
      write(*,*) "with ", mpi_size, "processes and ", nz, " voxels in the Z direction, process ", mpi_rank, " will have no work to do. Abort"
      write(*,*) "number of slices for process ", mpi_rank, " is ", nz_local
      error_on_slice_distribution_proc = .true.
    endif

    call MPILogicalOR(error_on_slice_distribution_proc, error_on_slice_distribution)
    if (error_on_slice_distribution) then
      call finalize_mpi()
    endif

    do i=1,9
      this%data_container(i)%cdata = fftw_alloc_complex(alloc_local)
      call c_f_pointer(this%data_container(i)%cdata, this%data_container(i)%data, [this%nx,this%ny,this%nz_local])
    enddo

    do i=1,9
          this%plan(i) = fftw_mpi_plan_dft_3d(this%nz, this%ny, this%nx , &
               this%data_container(i)%data(:,:,:),this%data_container(i)%data(:,:,:), MPI_COMM_WORLD, fftw_fwd , FFTW_MEASURE ) !FFTW_IN_PLACE
          this%iplan(i)= fftw_mpi_plan_dft_3d(this%nz, this%ny , this%nx , &
               this%data_container(i)%data(:,:,:),this%data_container(i)%data(:,:,:), MPI_COMM_WORLD, fftw_bkwd , FFTW_MEASURE )
    end do

    this%initialized = .true.

    CALL MPI_BARRIER(MPI_COMM_WORLD, my_mpi_err)
    if (mpi_rank.eq.mpi_master_rank) write(*,*) "all processes succesfully initialized fftw"

    call MPI_AllGather(nz_local, 1, MPI_INTEGER, &
                     mpi_local_z_all_proc, 1, MPI_INTEGER, &
                     MPI_COMM_WORLD, my_mpi_err)


    write(*,*) mpi_rank, "fftw allocation details"
    write(*,*) mpi_rank, "alloc_local ", alloc_local
    write(*,*) mpi_rank, "this%nx ", this%nx
    write(*,*) mpi_rank, "this%ny ", this%ny
    write(*,*) mpi_rank, "this%nz ", this%nz
    write(*,*) mpi_rank, "this%nz_local ", this%nz_local
    write(*,*) mpi_rank, "this%z_offset ", this%z_offset
    write(*,*) mpi_rank, "nz_start ", nz_start
    write(*,*) mpi_rank, "nz_end ", nz_end
    write(*,*) mpi_rank, "mpi_local_z_all_proc", mpi_local_z_all_proc

  end subroutine

  function getFFTDataAtPointComplexTensor2(this, ix, iy, iz) result(cmplx_tensor2)
    class(fft_grid_data_type), intent(in) :: this
    integer, intent(in) :: ix,iy,iz
    complex(k_real) :: cmplx_tensor2(3,3)
    integer :: i,j, c
    c = 0
    do j=1,3
      do i=1,3
        c=c+1
        cmplx_tensor2(i,j) = this%data_container(c)%data(ix,iy,iz)
      enddo
    enddo
  end function

  function getFFTDataAtPointComplexScalar(this, ix, iy, iz) result(cmplx_scalar)
    class(fft_grid_data_type), intent(in) :: this
    integer, intent(in) :: ix,iy,iz
    complex(k_real) :: cmplx_scalar

    cmplx_scalar = this%data_container(this%scalar_index)%data(ix,iy,iz)

  end function

  function getFFTDataAtPointRealTensor2(this, ix, iy, iz) result(tensor2)
    class(fft_grid_data_type), intent(in) :: this
    integer, intent(in) :: ix,iy,iz
    real(k_real) :: tensor2(3,3)
    integer :: i,j, c
    c = 0
    do j=1,3
      do i=1,3
        c=c+1
        tensor2(i,j) = this%data_container(c)%data(ix,iy,iz)%re
      enddo
    enddo
  end function

  function getFFTDataAtPointRealScalar(this, ix, iy, iz) result(scalar)
    class(fft_grid_data_type), intent(in) :: this
    integer, intent(in) :: ix,iy,iz
    real(k_real) :: scalar

    scalar = this%data_container(this%scalar_index)%data(ix,iy,iz)%re

  end function

  subroutine setFFTDataAtPointFromComplexTensor2(this, ix, iy, iz, cmplx_tensor2)
    class(fft_grid_data_type), intent(inout) :: this
    integer, intent(in) :: ix,iy,iz
    complex(k_real), intent(in) :: cmplx_tensor2(3,3)
    integer :: i, j, c

    c = 0
    do j=1,3
      do i=1,3
        c=c+1
        this%data_container(c)%data(ix,iy,iz) = cmplx_tensor2(i,j)
      enddo
    enddo

  end subroutine

  subroutine setFFTDataAtPointFromComplexScalar(this, ix, iy, iz, cmplx_scalar)
    class(fft_grid_data_type), intent(inout) :: this
    integer, intent(in) :: ix,iy,iz
    complex(k_real), intent(in) :: cmplx_scalar

    this%data_container(this%scalar_index)%data(ix,iy,iz) = cmplx_scalar

  end subroutine

  subroutine setFFTDataAtPointFromRealTensor2(this, ix, iy, iz, tensor2)
    class(fft_grid_data_type), intent(inout) :: this
    integer, intent(in) :: ix,iy,iz
    real(k_real), intent(in) :: tensor2(3,3)
    integer :: i, j, c

    c = 0
    do j=1,3
      do i=1,3
        c=c+1
        this%data_container(c)%data(ix,iy,iz)= cmplx( tensor2(i,j), 0., k_real)
      enddo
    enddo

  end subroutine

  subroutine setFFTDataAtPointFromRealScalar(this, ix, iy, iz, scalar)
    class(fft_grid_data_type), intent(inout) :: this
    integer, intent(in) :: ix,iy,iz
    real(k_real), intent(in) :: scalar

    this%data_container(this%scalar_index)%data(ix,iy,iz)= cmplx( scalar, 0., k_real)

  end subroutine

  subroutine executeDFTTensor2(this)
    class(fft_grid_data_type), intent(inout) :: this
    integer :: i
    ! compute the FFT transform of the polarization field
    do i=1,9
      call fftw_mpi_execute_dft(this%plan(i), this%data_container(i)%data(:,:,:), this%data_container(i)%data(:,:,:))
    end do
  end subroutine

  subroutine executeDFTScalar(this)
    class(fft_grid_data_type), intent(inout) :: this
    call fftw_mpi_execute_dft(this%plan(this%scalar_index), this%data_container(this%scalar_index)%data(:,:,:), this%data_container(this%scalar_index)%data(:,:,:))
  end subroutine

  subroutine executeInverseDFTTensor2(this)
    class(fft_grid_data_type), intent(inout) :: this
    integer :: i
    ! compute the FFT transform of the polarization field
    do i=1,9
      call fftw_mpi_execute_dft(this%iplan(i), this%data_container(i)%data(:,:,:), this%data_container(i)%data(:,:,:))
    end do
  end subroutine

  subroutine executeInverseDFTScalar(this)
    class(fft_grid_data_type), intent(inout) :: this
    ! compute the FFT transform of the polarization field
    call fftw_mpi_execute_dft(this%iplan(this%scalar_index), this%data_container(this%scalar_index)%data(:,:,:), this%data_container(this%scalar_index)%data(:,:,:))
  end subroutine

  subroutine destroyPlans(this)
    class(fft_grid_data_type), intent(inout) :: this
    integer :: i
    do i=1,9
      call fftw_destroy_plan(this%plan(i))
      call fftw_destroy_plan(this%iplan(i))
    end do
  end subroutine
end module
