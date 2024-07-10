module all_mighty_grid_mod
use kinds
use all_grid_data_mod, only : all_grid_data, all_grid_data_multi_phase
use fft_grid_data_mod, only : fft_grid_data_type
use, intrinsic :: iso_c_binding
use grid_data_types, only : var_dimension_type
use string_module
implicit none

type all_mighty_grid_type
    integer :: n_phases = -1
    type(all_grid_data), pointer :: common_grid => null()
    type(all_grid_data_multi_phase), pointer :: per_phase_grid => null()
    type(fft_grid_data_type), pointer :: fft_grid_data => null()
    integer :: nx, ny, nz, num_voxel
    integer ::  nx_rank, ny_rank, nz_rank, &
                x_offset_rank, x_start_rank, x_end_rank, &
                y_offset_rank, y_start_rank, y_end_rank, &
                z_offset_rank, z_start_rank, z_end_rank
    logical :: dimensions_initialized= .false.
    logical :: num_slip_systems_initialized = .false.
    class(var_dimension_type), pointer :: dimension_obj => null()



    ! write field option variables
    integer :: n_file_to_keep = 0
    integer :: dump_every_n_steps = 100 
    type(string_array) :: restart_file_names, xmf_restart_file_names
    type(string_type) :: restart_file_base_name
    integer :: write_fields_every_n_steps = 0
    type(string_type) :: field_file_base_name

    contains 
    procedure :: AMGInit
    procedure :: AMGSetDimension
    procedure :: AMGSetNumSlipSystems
    procedure :: AMGgetCommonGridPtr
    procedure :: AMGgetFFTGridPtr
    procedure :: AMGgetPhaseGridPtrByIdx
    procedure :: AMGallocateGridVariables
    procedure :: AMGgetLoopLimitsRank
    procedure :: AMGGetLoopLimitsGlobal
    procedure :: AMGGetRankBoxLimit
    procedure :: AMGgetGlobalGridDimension
    procedure :: AMGgetNumVoxel
    procedure :: AMGgetAllPhaseTensor2
    procedure :: AMGGetAllPhaseMatrix66
    procedure :: AMGsetDumpForRestartOptions
    procedure :: AMGsetWriteFieldOptions
    procedure :: AMGReloadFromDump
    procedure :: AMGDumpForRestart
    procedure :: AMGUpdateAllStatefulVariables
    procedure :: AMGResetAllStatefulVariables
    procedure :: AMGXIndexdxLocal2Global
    procedure :: AMGYIndexdxLocal2Global
    procedure :: AMGZIndexdxLocal2Global
    


end type


type rank2SinglePhase
    real(k_real), dimension(:,:,:,:,:), pointer :: data => null()
end type

type rank2MultiPhase
    type(rank2SinglePhase), dimension(:), pointer :: ph => null()
    contains
    procedure :: FreeMemory => FreeMemoryRank2MultiPhase
end type

type Matrix66SinglePhase
    real(k_real), dimension(:,:,:,:,:), pointer :: data => null()
end type

type Matrix66MultiPhase
    type(Matrix66SinglePhase), dimension(:), pointer :: ph => null()
    contains
    procedure :: FreeMemory => FreeMemoryMatrix66MultiPhase
end type

contains 

subroutine FreeMemoryRank2MultiPhase(this)
    implicit none
    class(rank2MultiPhase), intent(inout) :: this
    deallocate(this%ph)
end subroutine

subroutine FreeMemoryMatrix66MultiPhase(this)
    implicit none
    class(Matrix66MultiPhase), intent(inout) :: this
    deallocate(this%ph)
end subroutine


subroutine AMGinit(this, n_phases)
    implicit none
    class(all_mighty_grid_type), intent(inout) :: this
    integer, intent(in) :: n_phases
    integer :: idx
    class(all_grid_data), pointer :: temp_grid => null()
    this%n_phases = n_phases
    allocate(this%fft_grid_data)
    allocate(this%common_grid)
    allocate(this%per_phase_grid)
    do idx = 1, n_phases
        allocate(temp_grid)
        call this%per_phase_grid%addElement(temp_grid)
        nullify(temp_grid)
    enddo
end subroutine

function AMGXIndexdxLocal2Global(this, ix_local) result(ix_global)
    implicit none
    class(all_mighty_grid_type), intent(in) :: this
    integer, intent(in) :: ix_local
    integer ::  ix_global
    ix_global = ix_local + int(this%x_offset_rank)
end function

function AMGYIndexdxLocal2Global(this, iy_local) result(iy_global)
    implicit none
    class(all_mighty_grid_type), intent(in) :: this
    integer, intent(in) :: iy_local
    integer ::  iy_global
    iy_global = iy_local + int(this%y_offset_rank)
end function

function AMGZIndexdxLocal2Global(this, iz_local) result(iz_global)
    implicit none
    class(all_mighty_grid_type), intent(in) :: this
    integer, intent(in) :: iz_local
    integer ::  iz_global
    iz_global = iz_local + int(this%z_offset_rank)
end function

subroutine AMGGetCommonGridPtr(this, c_grid_ptr)
    implicit none
    class(all_mighty_grid_type), intent(in) :: this
    class(all_grid_data), pointer, intent(inout) :: c_grid_ptr
    if (associated(c_grid_ptr)) error stop "getCommonGridPtr: c_grid_ptr already associated"
    if (.not.associated(this%common_grid)) error stop "getCommonGridPtr: this%common_grid not associated"
    c_grid_ptr => this%common_grid
end subroutine

subroutine AMGGetFFTGridPtr(this, fft_grid_ptr)
    implicit none
    class(all_mighty_grid_type), intent(in) :: this
    type(fft_grid_data_type), pointer, intent(inout) :: fft_grid_ptr
    if (associated(fft_grid_ptr)) error stop "getFFTGridPtr: fft_grid_ptr already associated"
    if (.not.associated(this%fft_grid_data)) error stop "getFFTGridPtr: this%fft_grid_data not associated"
    fft_grid_ptr => this%fft_grid_data
end subroutine

subroutine AMGGetPhaseGridPtrByIdx(this, ph_idx, ph_grid_ptr)
    implicit none
    class(all_mighty_grid_type), intent(in) :: this
    integer, intent(in) :: ph_idx
    class(all_grid_data), pointer, intent(inout) :: ph_grid_ptr
    if (associated(ph_grid_ptr)) error stop "getPhaseGridPtrByIdx: ph_grid_ptr already associated"
    call this%per_phase_grid%getElementPtr(ph_idx, ph_grid_ptr)
end subroutine

subroutine AMGSetDimension(this, n_voxel_xyz )
    implicit none
    class(all_mighty_grid_type), intent(inout) :: this
    integer, intent(in) :: n_voxel_xyz(3)
    class(all_grid_data), pointer :: temp_grid => null()
    integer :: idx
    integer :: xyz_start_rank(3), &
               xyz_end_rank(3), &
               xyz_offset_rank(3), &
               nx_ny_nz_rank(3)
    associate(nx => n_voxel_xyz(1), ny => n_voxel_xyz(2), nz => n_voxel_xyz(3))
    this%nx = nx; this%ny = ny; this%nz = nz;
    this%num_voxel = product(n_voxel_xyz)

    call this%fft_grid_data%initFFTData(nx,ny,nz, this%nz_rank, this%z_offset_rank, this%z_start_rank, this%z_end_rank) 


    this%nx_rank = this%nx
    this%x_offset_rank=0;
    this%x_start_rank=this%x_offset_rank+1; 
    this%x_end_rank=this%x_offset_rank + this%nx_rank

    this%ny_rank = this%ny
    this%y_offset_rank=0;
    this%y_start_rank=this%y_offset_rank+1; 
    this%y_end_rank=this%y_offset_rank + this%ny_rank
    
    xyz_start_rank(1) = this%x_start_rank
    xyz_start_rank(2) = this%y_start_rank
    xyz_start_rank(3) = this%z_start_rank

    xyz_end_rank(1) = this%x_end_rank
    xyz_end_rank(2) = this%y_end_rank
    xyz_end_rank(3) = this%z_end_rank

    xyz_offset_rank(1) = this%x_offset_rank
    xyz_offset_rank(2) = this%y_offset_rank
    xyz_offset_rank(3) = this%z_offset_rank
    
    nx_ny_nz_rank(1) = this%nx_rank
    nx_ny_nz_rank(2) = this%ny_rank
    nx_ny_nz_rank(3) = this%nz_rank

    allocate(this%dimension_obj)
    call this%dimension_obj%init(n_voxel_xyz, nx_ny_nz_rank, &
                            xyz_start_rank, xyz_end_rank, xyz_offset_rank)

    ! call this%common_grid%setDimensions( n_voxel_xyz, xyz_start_rank, xyz_end_rank, xyz_offset_rank)
    call this%common_grid%setDimensions(this%dimension_obj)

    do idx = 1, this%n_phases
        call this%per_phase_grid%getElementPtr(idx, temp_grid)
        call temp_grid%setDimensions( this%dimension_obj)
        nullify(temp_grid)
    enddo
    end associate

    this%dimensions_initialized =.true.
end subroutine


subroutine AMGSetNumSlipSystems(this, n_ss )
    implicit none
    class(all_mighty_grid_type), intent(inout) :: this
    integer, dimension(:) :: n_ss
    class(all_grid_data), pointer :: temp_grid => null()
    integer :: idx

    call this%common_grid%setNumSlipSystems(1)

    do idx = 1, this%n_phases
        call this%per_phase_grid%getElementPtr(idx, temp_grid)
        call temp_grid%setNumSlipSystems(n_ss(idx))
        nullify(temp_grid)
    enddo

    this%num_slip_systems_initialized =.true.
end subroutine

subroutine AMGAllocateGridVariables(this)
    implicit none
    class(all_mighty_grid_type), intent(inout) :: this
    class(all_grid_data), pointer :: temp_grid => null()
    integer :: idx

    if (.not.(this%dimensions_initialized)) error stop "AMGAllocateGridVariables: cannot allocate field variables, global dimensions not initialized"
    if (.not.(this%num_slip_systems_initialized)) error stop "AMGAllocateGridVariables: cannot allocate field variables, number of slip systems not initialized"

    
    call this%common_grid%init()
    do idx = 1, this%n_phases
        call this%per_phase_grid%getElementPtr(idx, temp_grid)
        call temp_grid%init(phase_id=idx)
        nullify(temp_grid)
    enddo
    
end subroutine


subroutine AMGGetRankBoxLimit(this, x_box_start_rank, x_box_end_rank, y_box_start_rank, y_box_end_rank, z_box_start_rank, z_box_end_rank)
    implicit none
    class(all_mighty_grid_type), intent(inout) :: this
    integer, intent(out) :: x_box_start_rank, x_box_end_rank, y_box_start_rank, y_box_end_rank, z_box_start_rank, z_box_end_rank

    if (.not.(this%dimensions_initialized)) error stop "AMGGetRankBoxLimit: dimensions have not been initailized. Abort!"

    x_box_start_rank = this%x_start_rank
    y_box_start_rank = this%y_start_rank
    z_box_start_rank = this%z_start_rank
    x_box_end_rank = this%x_end_rank
    y_box_end_rank = this%y_end_rank
    z_box_end_rank = this%z_end_rank
    
end subroutine

subroutine AMGGetLoopLimitsRank(this, x_start_rank, x_end_rank, y_start_rank, y_end_rank, z_start_rank, z_end_rank)
    implicit none
    class(all_mighty_grid_type), intent(inout) :: this
    integer, intent(out) ::   x_start_rank, x_end_rank, y_start_rank, y_end_rank, z_start_rank, z_end_rank

    if (.not.(this%dimensions_initialized)) error stop "AMGGetLoopLimitsRank: dimensions have not been initailized. Abort!"

    x_start_rank = 1
    y_start_rank = 1
    z_start_rank = 1
    x_end_rank = this%nx_rank
    y_end_rank = this%ny_rank
    z_end_rank = this%nz_rank
    
end subroutine

subroutine AMGGetLoopLimitsGlobal(this, x_start, x_end, y_start, y_end, z_start, z_end)
    implicit none
    class(all_mighty_grid_type), intent(inout) :: this
    integer, intent(out) ::  x_start, x_end, y_start, y_end, z_start, z_end

    x_start = 1
    y_start = 1
    z_start = 1
    x_end = this%nx
    y_end = this%ny
    z_end = this%nz
    
end subroutine

subroutine AMGGetGlobalGridDimension(this, nx, ny, nz)
    implicit none
    class(all_mighty_grid_type), intent(inout) :: this
    integer, intent(out) ::  nx, ny, nz

    nx = this%nx 
    ny = this%ny 
    nz = this%nz
    
end subroutine

function AMGGetNumVoxel(this) result(num_voxel)
    implicit none
    class(all_mighty_grid_type), intent(inout) :: this
    integer ::  num_voxel
    num_voxel = this%num_voxel
    
end function

subroutine AMGGetAllPhaseTensor2(this, var_name, all_phase_rank_2)
    implicit none
    class(all_mighty_grid_type), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    class(rank2MultiPhase), intent(inout) :: all_phase_rank_2 
    class(all_grid_data), pointer :: temp_grid
    integer :: idx

    if (associated(all_phase_rank_2%ph)) error stop "AMGGetAllPhaseTensor2:, all_phase_rank_2%data_all already associated"
    allocate(all_phase_rank_2%ph(this%n_phases))

    do idx = 1, this%n_phases
        call this%per_phase_grid%getElementPtr(idx, temp_grid)
        call temp_grid%getTensor2DataPointerByName(var_name, all_phase_rank_2%ph(idx)%data)
        nullify(temp_grid)
    enddo
    
end subroutine


subroutine AMGGetAllPhaseMatrix66(this, var_name, all_phase_matrix66)
    implicit none
    class(all_mighty_grid_type), intent(inout) :: this
    character(len=*), intent(in) :: var_name
    class(Matrix66MultiPhase), intent(inout) :: all_phase_matrix66 
    class(all_grid_data), pointer :: temp_grid
    integer :: idx

    if (associated(all_phase_matrix66%ph)) error stop "AMGGetAllPhaseMatrix66:, all_phase_matrix66%data_all already associated"
    allocate(all_phase_matrix66%ph(this%n_phases))

    do idx = 1, this%n_phases
        call this%per_phase_grid%getElementPtr(idx, temp_grid)
        call temp_grid%getMatrix66DataPointerByName(var_name, all_phase_matrix66%ph(idx)%data)
        nullify(temp_grid)
    enddo
    
end subroutine

subroutine AMGSetDumpForRestartOptions(this, n_file_to_keep, dump_every_n_steps, restart_file_base_name )
    implicit none
    class(all_mighty_grid_type), intent(inout) :: this
    integer, intent(in) :: n_file_to_keep, dump_every_n_steps
    character(len=*) :: restart_file_base_name
    class(all_grid_data), pointer :: temp_grid
    integer :: idx

    this%n_file_to_keep=n_file_to_keep
    this%dump_every_n_steps= dump_every_n_steps
    call this%restart_file_base_name%setString(restart_file_base_name)

    call this%common_grid%setDumpForRestartOptions(n_file_to_keep, dump_every_n_steps, restart_file_base_name)

    do idx = 1, this%n_phases
        call this%per_phase_grid%getElementPtr(idx, temp_grid)
        call temp_grid%setDumpForRestartOptions(n_file_to_keep, dump_every_n_steps, restart_file_base_name)
        nullify(temp_grid)
    enddo

end subroutine


module subroutine AMGDumpForRestart(this, time, dt, step, XLSEC, XLSEC_old, force_write, write_only_current_value)
    use mpi_variables_mod, only : i_am_mpi_master
    use number_to_string_mod, only : int2string
    use mpi, only : MPI_COMM_WORLD, MPI_INFO_NULL, MPI_Barrier
    use read_write_parallel_h5_mod, only : createAndOpenH5EmptyFile, CloseH5File, addRealAttribute, addIntegerAttribute, addRealMatrixAttribute
    use hdf5
    use xmf_writer_mod
    implicit none
    class(all_mighty_grid_type), intent(inout) :: this
    real(k_real), intent(in) :: time, dt
    integer, intent(in) :: step
    real(k_real), intent(in) :: XLSEC(6,6), XLSEC_old(6,6)
    logical, intent(in) :: force_write, write_only_current_value
    class(all_grid_data), pointer :: temp_grid
    integer :: idx, i, ierr
    character(len=100) :: file_name, xmf_filename
    INTEGER(HID_T) :: plist_id1, file_id
    type(xmf_writer_type) :: xmf_writer

    if (.not.write_only_current_value) then
        if (this%dump_every_n_steps.eq.0) return ! 0 means not write in any case
        if (.not.(force_write)) then
            if (this%n_file_to_keep==0) return
            if (step==0) return
            if (mod(step,this%dump_every_n_steps).ne.0) return
        endif
        ! if we reached this place we need to write a dump file
        file_name = this%restart_file_base_name%getString()//"_"//trim(adjustl(int2string(step)))//".h5"
        xmf_filename = this%restart_file_base_name%getString()//"_"//trim(adjustl(int2string(step)))//".xmf"
        
        ! check if we need to delete a dump file before writing the new one
        if (this%restart_file_names%getNumStrings().ge.this%n_file_to_keep) then
            ! delete the oldest file
            call MPI_Barrier(MPI_COMM_WORLD, ierr)
            if (i_am_mpi_master) call system('rm '//this%restart_file_names%getStringByIndex(1))
            call MPI_Barrier(MPI_COMM_WORLD, ierr)
            if (i_am_mpi_master) call system('rm '//this%xmf_restart_file_names%getStringByIndex(1))
            call MPI_Barrier(MPI_COMM_WORLD, ierr)
            if (this%n_file_to_keep.gt.1) then
            do i=2, this%n_file_to_keep
                call this%restart_file_names%setStringByIndex(i-1, this%restart_file_names%getStringByIndex(i))
                call this%xmf_restart_file_names%setStringByIndex(i-1, this%xmf_restart_file_names%getStringByIndex(i))
            enddo
            endif
            call this%restart_file_names%setStringByIndex(this%n_file_to_keep, file_name)
            call this%xmf_restart_file_names%setStringByIndex(this%n_file_to_keep, xmf_filename)
        else
            ! set the first entry to the actual name
            call this%restart_file_names%addString(file_name)
            call this%xmf_restart_file_names%addString(xmf_filename)
        endif
  
    else
        ! we are not writing a dump file but a field file
        if (this%write_fields_every_n_steps.eq.0) return ! 0 means not write in any case
        if (.not.(force_write)) then
            if (mod(step,this%write_fields_every_n_steps).ne.0) return
        endif
        ! if we reached this place we need to write a field file
        file_name = this%field_file_base_name%getString()//"_"//trim(adjustl(int2string(step)))//".h5"
        xmf_filename = this%field_file_base_name%getString()//"_"//trim(adjustl(int2string(step)))//".xmf"
    endif
    
    !if we reached this point we want to write the file

    ! first create the new file and close it
    call createAndOpenH5EmptyFile(trim(adjustl(file_name)), file_id, plist_id1)
    call addRealAttribute(file_id, "time", time)
    call addRealAttribute(file_id, "dt", dt)
    call addIntegerAttribute(file_id, "step", step)
    call addRealMatrixAttribute(file_id, "XLSEC", XLSEC)
    call addRealMatrixAttribute(file_id, "XLSEC_old", XLSEC_old)
    call CloseH5File(file_id, plist_id1)

    !write xmf
    call xmf_writer%createNewXMF(trim(adjustl(xmf_filename)), this%dimension_obj%nx_ny_nz)
    call xmf_writer%writeXMFHeader(time)
    call this%common_grid%AGDDumpForRestart(trim(adjustl(file_name)), xmf_writer, &
                                    write_only_current_value)

    do idx = 1, this%n_phases
        call this%per_phase_grid%getElementPtr(idx, temp_grid)
        call temp_grid%AGDDumpForRestart(trim(adjustl(file_name)), xmf_writer, &
                                        write_only_current_value)
        nullify(temp_grid)
    enddo
    call xmf_writer%writeXMFEndOfFile()
    call xmf_writer%closeXMFFile()

  end subroutine

subroutine AMGSetWriteFieldOptions(this, write_fields_every_n_steps, field_file_base_name )
    use number_to_string_mod, only : int2string
    implicit none
    class(all_mighty_grid_type), intent(inout) :: this
    integer, intent(in) :: write_fields_every_n_steps
    character(len=*) :: field_file_base_name
    class(all_grid_data), pointer :: temp_grid
    integer :: idx

    this%write_fields_every_n_steps= write_fields_every_n_steps
    call this%field_file_base_name%setString(field_file_base_name)

    call this%common_grid%setWriteFieldOptions(write_fields_every_n_steps, field_file_base_name)

    do idx = 1, this%n_phases
        call this%per_phase_grid%getElementPtr(idx, temp_grid)
        call temp_grid%setWriteFieldOptions(write_fields_every_n_steps, field_file_base_name)
        nullify(temp_grid)
    enddo

end subroutine

subroutine AMGReloadFromDump(this, file_name, time, dt, step, XLSEC, XLSEC_old )
    implicit none
    class(all_mighty_grid_type), intent(inout) :: this
    character(len=*), intent(in) :: file_name
    real(k_real), intent(out) :: time, dt
    integer, intent(out) :: step
    real(k_real), intent(out) :: XLSEC(6,6), XLSEC_old(6,6)
    class(all_grid_data), pointer :: temp_grid
    integer :: idx
    call this%common_grid%AGDReloadFromDump(file_name, time, dt, step, XLSEC, XLSEC_old)

    do idx = 1, this%n_phases
        call this%per_phase_grid%getElementPtr(idx, temp_grid)
        call temp_grid%AGDReloadFromDump(file_name, time, dt, step, XLSEC, XLSEC_old)
        nullify(temp_grid)
    enddo

end subroutine

subroutine AMGUpdateAllStatefulVariables(this )
    implicit none
    class(all_mighty_grid_type), intent(inout) :: this
    class(all_grid_data), pointer :: temp_grid
    integer :: idx
    call this%common_grid%AGDupdateAllStatefulVariables()

    do idx = 1, this%n_phases
        call this%per_phase_grid%getElementPtr(idx, temp_grid)
        call temp_grid%AGDupdateAllStatefulVariables()
        nullify(temp_grid)
    enddo

end subroutine

subroutine AMGResetAllStatefulVariables(this )
    implicit none
    class(all_mighty_grid_type), intent(inout) :: this
    class(all_grid_data), pointer :: temp_grid
    integer :: idx
    call this%common_grid%AGDResetAllStatefulVariables()

    do idx = 1, this%n_phases
        call this%per_phase_grid%getElementPtr(idx, temp_grid)
        call temp_grid%AGDResetAllStatefulVariables()
        nullify(temp_grid)
    enddo

end subroutine

end module