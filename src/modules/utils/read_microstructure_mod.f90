module read_microstructure_file_mod
use kinds

contains

subroutine read_microstructure_from_file(microstructure_fname, all_mighty_grid_in)
  use all_grid_data_mod, only :all_grid_data
  use all_mighty_grid_mod, only : all_mighty_grid_type, rank2MultiPhase
  use bunge_mod, only : bungeRotationMatrixCrystal2Sample
  use units_conversion_mod, only : deg2rad
  use gb_normals_mod, only : initGBproperties
  use read_from_file_utils, only : file_reader
  use string_module, only : string_array
  implicit none
  character(len=*), intent(in) :: microstructure_fname
  type(all_mighty_grid_type), intent(inout) :: all_mighty_grid_in
  class(all_grid_data), pointer :: common_grid_data => null()
  type(rank2MultiPhase) :: R_crystal2sample
  integer, pointer, dimension(:,:,:) :: grain_id=>null()
  real(k_real), pointer, dimension(:,:,:,:) :: phase_fraction=>null()
  integer :: ix, iy, iz, iix, iiy, iiz, z_local, &
             nx,ny,nz, n_phases, iph, counter, num_columns

  real(k_real) :: ph, th, om
  real(k_real), allocatable, dimension(:) :: read_buffer
  ! the following arrays are need for gb normals calcualtion
  integer, allocatable :: grain_id_full(:,:,:)
  type(file_reader) :: microstructure_reader
  type(string_array) :: dummy_string_array

  nullify(grain_id, phase_fraction, common_grid_data)
  call all_mighty_grid_in%AMGGetCommonGridPtr(common_grid_data)
  call common_grid_data%getGenericVectorDataPointerByName("phase_fraction", phase_fraction)
  call common_grid_data%getScalarIntegerDataPointerByName("grain_id", grain_id)
  call all_mighty_grid_in%AMGGetAllPhaseTensor2("R_crystal2sample", R_crystal2sample)

  call all_mighty_grid_in%AMGGetGlobalGridDimension(nx,ny,nz)
    allocate(grain_id_full(nx,ny,nz))
    n_phases = all_mighty_grid_in%n_phases
    !x,y,z, gid, (fraction, alpha, beta, gamma)*n_phases
    num_columns = 3 + 1 + (1 + 3)*n_phases
    allocate(read_buffer(num_columns))

    call microstructure_reader%openReadTextFile(microstructure_fname)
    call microstructure_reader%readLineAndCheckStringsAreEqual("LApx-wPF", smart_string_array=dummy_string_array, &
              extra_error_string="convert old input to new fromat using convert_old_microstructure_to_LApx-PF.py script")

    do iz=1,nz; do iy=1,ny; do ix=1,nx

      call microstructure_reader%readVector(num_columns, read_buffer)
      iix = real2int(read_buffer(1))
      iiy = real2int(read_buffer(2))
      iiz = real2int(read_buffer(3))  
      
      if ((iix.ne.ix).or.(iiy.ne.iy).or.(iiz.ne.iz)) then
        write(*,*) "the indeces in the mictrostrucure file are not properlyordered"
        write(*,*) "I'm expecting the point ", ix, iy, iz, "instead i have", iix, iiy, iiz
        error stop "check your microstructure file"
      endif

      grain_id_full(ix,iy,iz) = real2int(read_buffer(4))  
      z_local = common_grid_data%getZRankFromZGlobal(iz)
      if (z_local.gt.0) then
        grain_id(ix,iy,z_local) = real2int(read_buffer(4))
        counter = 4
        do iph = 1,n_phases
          counter = counter +1; phase_fraction(iph, ix,iy,z_local) = read_buffer(counter);
          counter = counter +1; ph = read_buffer(counter);
          counter = counter +1; th = read_buffer(counter);
          counter = counter +1; om = read_buffer(counter);
          call bungeRotationMatrixCrystal2Sample(ph*deg2rad,th*deg2rad,om*deg2rad,  R_crystal2sample%ph(iph)%data(:,:,ix,iy,z_local))
        enddo

        if (any(phase_fraction(:, ix,iy,z_local).lt.0._k_real)) then
          write(*,*) "Error Reading the microstructure file!"
          write(*,*) "one of the phase fraction for point ", ix, iy, iz, " is negative:", phase_fraction(:, ix,iy,z_local)
          error stop "Abort"
        endif
        if (abs(sum(phase_fraction(:, ix,iy,z_local) ) -1._k_real) .gt. 1e-10_k_real ) then
          write(*,*) "Error Reading the microstructure file!"
          write(*,*) "the total phase fraciton for point ", ix, iy, iz, " is not equal to 1, instead it is ", sum(phase_fraction(:, ix,iy,z_local))
          error stop "Abort"
        endif
      endif

    enddo; enddo; enddo
    call microstructure_reader%closeTextFile()
    call initGBproperties(common_grid_data, grain_id_full)
    deallocate(grain_id_full)

  call R_crystal2sample%FreeMemory()
  nullify(grain_id, phase_fraction, common_grid_data)

end subroutine

end module
