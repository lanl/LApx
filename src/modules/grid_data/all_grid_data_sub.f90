submodule(all_grid_data_mod) all_grid_data_sub

contains

  module function getGlobalGridNPoints(this) result(npoints)
    implicit none
    class(all_grid_data), intent(in) :: this
    integer :: npoints
    npoints = product(this%nx_ny_nz)
  end function

  module function getGlobalGridDimension(this) result(nx_ny_nz)
    implicit none
    class(all_grid_data), intent(in) :: this
    integer, dimension(3) :: nx_ny_nz
    nx_ny_nz = this%nx_ny_nz
  end function

  module function getRankGridDimension(this) result(nx_ny_nz)
    implicit none
    class(all_grid_data), intent(in) :: this
    integer, dimension(3) :: nx_ny_nz
    nx_ny_nz = this%nx_ny_nz_rank
  end function

  module function getRankZStart(this) result(z_start_rank)
    implicit none
    class(all_grid_data), intent(in) :: this
    integer :: z_start_rank
    z_start_rank = this%z_start_rank
  end function

  module function getRankZEnd(this) result(z_end_rank)
    implicit none
    class(all_grid_data), intent(in) :: this
    integer :: z_end_rank
    z_end_rank = this%z_end_rank
  end function

  module function getZRankFromZGlobal(this, z_global) result(z_local)
    implicit none
    class(all_grid_data), intent(in) :: this
    integer, intent(in):: z_global
    integer :: z_local
    z_local = -1
    if (z_global>=this%z_start_rank.and.z_global<=this%z_end_rank) z_local = z_global - this%z_offset_rank
  end function

  module function getZGlobalFromZRank(this, z_local) result(z_global)
    implicit none
    class(all_grid_data), intent(in) :: this
    integer, intent(in):: z_local
    integer :: z_global
    z_global = z_local + this%z_offset_rank
  end function

  module subroutine  setDimensionsAGD(this, dimension_obj) !n_points_global, xyz_start_rank, xyz_end_rank, xyz_offset_rank)
    use grid_data_types, only : var_dimension_type
    implicit none
    class(all_grid_data), intent(inout) :: this
    class(var_dimension_type), intent(in), pointer :: dimension_obj
    ! integer, intent(in) :: n_points_global(3), &
    !                         xyz_start_rank(3), &
    !                         xyz_end_rank(3), &
    !                         xyz_offset_rank(3)
    if (.not.(associated(dimension_obj))) error stop "setDimensionsAGD: dimension_obj not associated"
    
    this%dimension_obj => dimension_obj
    this%nx_ny_nz = dimension_obj%nx_ny_nz
    this%nx_ny_nz_rank = dimension_obj%nx_ny_nz_rank

    this%xyz_start_rank = dimension_obj%xyz_start_rank
    this%xyz_end_rank = dimension_obj%xyz_end_rank
    this%xyz_offset_rank = dimension_obj%xyz_offset_rank

    this%x_offset_rank = this%xyz_offset_rank(1)
    this%x_start_rank = this%xyz_start_rank(1)
    this%x_end_rank = this%xyz_end_rank(1)

    this%y_offset_rank = this%xyz_offset_rank(2)
    this%y_start_rank = this%xyz_start_rank(2)
    this%y_end_rank = this%xyz_end_rank(2)

    this%z_offset_rank = this%xyz_offset_rank(3)
    this%z_start_rank = this%xyz_start_rank(3)
    this%z_end_rank = this%xyz_end_rank(3)

    this%n_points = product(this%nx_ny_nz)
    this%n_points_rank = product(this%nx_ny_nz_rank)
    this%voxel_weight = 1._k_real/int2real(this%n_points)
    this%dimensions_were_set =.true.

  end subroutine

  module subroutine  setNumSlipSystemsAGD(this, nss)
    class(all_grid_data), intent(inout) :: this
    integer, intent(in) :: nss
    this%n_ss = nss
    this%num_ss_was_set =.true.

  end subroutine

  module subroutine initAGD(this, phase_id)
    implicit none
    class(all_grid_data), intent(inout) :: this
    integer, intent(in), optional :: phase_id
    integer :: i, j, num_var_by_type
    integer(kind(grid_data_var_type_enum)) :: var_type
    type(grid_data_var_type), pointer :: var_ptr
    type(string_array) :: var_name_list
    integer, dimension(:), allocatable :: current_var_size
    integer :: stateful_level
    character(len=:), allocatable :: var_name

    if (.not.(this%dimensions_were_set)) error stop "initAGD cannot allocate grid because global diemnsion were not set "
    if (.not.(this%num_ss_was_set)) error stop "initAGD cannot allocate grid because the number of slip system was not set "
    this%phase_id = 0
    if (present(phase_id)) this%phase_id=phase_id

    this%all_scalar_has_data=.FALSE.
    this%all_scalar_integer_has_data=.FALSE.
    this%all_real_vector_has_data=.FALSE.
    this%all_vector5_has_data=.FALSE.
    this%all_vector6_has_data=.FALSE.
    this%all_generic_vector_has_data=.FALSE.
    this%all_ss_scalar_has_data=.FALSE.
    this%all_tensor2_has_data=.FALSE.
    this%all_matrix66_has_data=.FALSE.
    this%all_ss_generic_vector_has_data=.FALSE.
    this%all_ss_vector5_has_data=.FALSE.
    this%all_ss_vector6_has_data=.FALSE.
    this%all_ss_matrix_has_data=.FALSE.
    this%all_generic_matrix_has_data=.FALSE.

    associate(  nx => this%nx_ny_nz(1), &
                ny => this%nx_ny_nz(2), &
                nz_total => this%nx_ny_nz(3), &
                nz_rank => this%nx_ny_nz_rank(3), &
                z_offset_rank => this%z_offset_rank, &
                z_start_rank => this%z_offset_rank, &
                z_end_rank => this%z_offset_rank, &
                n_ss => this%n_ss )

    do i=1,size(all_grid_data_var_type_enum)
      var_type = all_grid_data_var_type_enum(i)
      num_var_by_type = this%all_grid_vars_list %countVarByType(var_type, var_name_list)

      if (num_var_by_type>0) then
        do j=1, num_var_by_type

          call this%all_grid_vars_list%checkVarExists(var_name_list%getStringByIndex(j), var_type, var_ptr)
          if (.not.(associated(var_ptr))) then
            write(*,*) "I can't find a variable I should be able to find!!"
            stop
          end if
          current_var_size = var_ptr%getSize()
          stateful_level = var_ptr%getStatefulLevel()
          var_name = var_ptr%getName()

          select case (var_type)
          case (scalar)
            if (j==1) allocate(this%all_scalar(num_var_by_type))
            call this%all_scalar(j)%init(var_name, this%dimension_obj, stateful_level)
            this%all_scalar_has_data=.true.

          case (scalar_integer)
            if (j==1) allocate(this%all_scalar_integer(num_var_by_type))
            call this%all_scalar_integer(j)%init(var_name, this%dimension_obj, stateful_level)
            this%all_scalar_integer_has_data=.true.

          case (real_vector)
            if (j==1) allocate(this%all_real_vector(num_var_by_type))
            call this%all_real_vector(j)%init(var_name, this%dimension_obj, stateful_level)
            this%all_real_vector_has_data=.true.

          case (vector5)
            if (j==1) allocate(this%all_vector5(num_var_by_type))
            call this%all_vector5(j)%init(var_name, this%dimension_obj, stateful_level)
            this%all_vector5_has_data=.true.

          case (vector6)
            if (j==1) allocate(this%all_vector6(num_var_by_type))
            call this%all_vector6(j)%init(var_name, this%dimension_obj, stateful_level)
            this%all_vector6_has_data=.true.

          case (generic_vector)
            if (j==1) allocate(this%all_generic_vector(num_var_by_type))
            call this%all_generic_vector(j)%init(var_name, this%dimension_obj, current_var_size(1), stateful_level)
            this%all_generic_vector_has_data=.true.

          case (ss_scalar)
            if (j==1) allocate(this%all_ss_scalar(num_var_by_type))
            call this%all_ss_scalar(j)%init(var_name, this%dimension_obj, n_ss, stateful_level)
            this%all_ss_scalar_has_data=.true.

          case (tensor2)
            if (j==1) allocate(this%all_tensor2(num_var_by_type))
            call this%all_tensor2(j)%init(var_name, this%dimension_obj, stateful_level)
            this%all_tensor2_has_data=.true.

          case (matrix66)
            if (j==1) allocate(this%all_matrix66(num_var_by_type))
            call this%all_matrix66(j)%init(var_name, this%dimension_obj, stateful_level)
            this%all_matrix66_has_data=.true.

          case (ss_generic_vector)
            if (j==1) allocate(this%all_ss_generic_vector(num_var_by_type))
            call this%all_ss_generic_vector(j)%init(var_name, this%dimension_obj, n_ss, current_var_size(1), stateful_level)
            this%all_ss_generic_vector_has_data=.true.

          case (ss_vector5)
            if (j==1) allocate(this%all_ss_vector5(num_var_by_type))
            call this%all_ss_vector5(j)%init(var_name, this%dimension_obj, n_ss, stateful_level)
            this%all_ss_vector5_has_data=.true.

          case (ss_vector6)
            if (j==1) allocate(this%all_ss_vector6(num_var_by_type))
            call this%all_ss_vector6(j)%init(var_name, this%dimension_obj, n_ss, stateful_level)
            this%all_ss_vector6_has_data=.true.

          case (ss_generic_matrix)
            if (j==1) allocate(this%all_ss_generic_matrix(num_var_by_type))
            call this%all_ss_generic_matrix(j)%init(var_name, this%dimension_obj, n_ss, current_var_size(1), current_var_size(2), stateful_level)
            this%all_ss_matrix_has_data=.true.

          case (generic_matrix)
            if (j==1) allocate(this%all_generic_matrix(num_var_by_type))
            call this%all_generic_matrix(j)%init(var_name, this%dimension_obj, current_var_size(1), current_var_size(2), stateful_level)
            this%all_generic_matrix_has_data=.true.

          case default
            write(*,*) "the selected variable type was not coded into initAGD!"
            stop
          end select

        end do
      end if
    end do
    end associate
  this%initialized = .TRUE.

  end subroutine

  module subroutine AGDDumpForRestart(this, file_name, xmf_writer, write_only_current_value)
    use hdf5
    use read_write_parallel_h5_mod, only : openExistingH5ReadWrite, createGroup, CloseH5File
    use mpi, only : MPI_COMM_WORLD, MPI_INFO_NULL, MPI_Barrier
    use number_to_string_mod, only : int2string
    use xmf_writer_mod, only :xmf_writer_type
    use string_module, only : string_type
    implicit none
    class(all_grid_data), intent(inout) :: this
    character(len=*), intent(in) :: file_name
    type(xmf_writer_type), intent(inout) :: xmf_writer
    INTEGER(HID_T) :: plist_id1, file_id
    integer :: i
    logical, intent(in) :: write_only_current_value
    type(string_type) :: suffix

    ! if (.not.write_only_current_value) then
    !   if (this%dump_every_n_steps.eq.0) return ! 0 means not write in any case
    !   if (.not.(force_write)) then
    !     if (this%n_file_to_keep==0) return
    !     if (step==0) return
    !     if (mod(step,this%dump_every_n_steps).ne.0) return
    !   endif
    !   ! if we reached this place we need to write a dump file
    !   file_name = this%restart_file_base_name%getString()//"_"//trim(adjustl(int2string(step)))//".h5"
    !   xmf_filename = this%restart_file_base_name%getString()//"_"//trim(adjustl(int2string(step)))//".xmf"
      
    !   ! check if we need to delete a dump file before writing the new one
    !   if (this%restart_file_names%getNumStrings().ge.this%n_file_to_keep) then
    !     ! delete the oldest file
    !     call MPI_Barrier(MPI_COMM_WORLD, ierr)
    !     if (i_am_mpi_master) call system('rm '//this%restart_file_names%getStringByIndex(1))
    !     call MPI_Barrier(MPI_COMM_WORLD, ierr)
    !     if (i_am_mpi_master) call system('rm '//this%xmf_restart_file_names%getStringByIndex(1))
    !     call MPI_Barrier(MPI_COMM_WORLD, ierr)
    !     if (this%n_file_to_keep.gt.1) then
    !       do i=2, this%n_file_to_keep
    !         call this%restart_file_names%setStringByIndex(i-1, this%restart_file_names%getStringByIndex(i))
    !         call this%xmf_restart_file_names%setStringByIndex(i-1, this%xmf_restart_file_names%getStringByIndex(i))
    !       enddo
    !     endif
    !     call this%restart_file_names%setStringByIndex(this%n_file_to_keep, file_name)
    !     call this%xmf_restart_file_names%setStringByIndex(this%n_file_to_keep, xmf_filename)
    !   else
    !     ! set the first entry to the actual name
    !     call this%restart_file_names%addString(file_name)
    !     call this%xmf_restart_file_names%addString(xmf_filename)
    !   endif

    ! else
    !   ! we are not writing a dump file but a field file
    !   if (this%write_fields_every_n_steps.eq.0) return ! 0 means not write in any case
    !   if (.not.(force_write)) then
    !     if (mod(step,this%write_fields_every_n_steps).ne.0) return
    !   endif
    !   ! if we reached this place we need to write a field file
    !   file_name = this%field_file_base_name%getString()//"_"//trim(adjustl(int2string(step)))//".h5"
    !   xmf_filename = this%field_file_base_name%getString()//"_"//trim(adjustl(int2string(step)))//".xmf"
    ! endif
    
    ! CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id1, ierr)
    ! CALL h5pset_fapl_mpio_f(plist_id1, MPI_COMM_WORLD, MPI_INFO_NULL, ierr)
    ! CALL h5fcreate_f(file_name, H5F_ACC_TRUNC_F, file_id, ierr, access_prp = plist_id1)
    ! call openExistingH5ReadWrite(file_name, file_id, plist_id1)

    if (this%phase_id>0) then
      call openExistingH5ReadWrite(file_name, file_id, plist_id1)
      call createGroup(this%phase_group_name%getString(), file_id)
      call CloseH5File(file_id, plist_id1)
    endif


    if (this%all_scalar_has_data) then
      do i=1,size(this%all_scalar)
        call this%all_scalar(i)%dumpForRestart(file_name, trim(adjustl(this%field_absolute_path%getString())), &
             only_current_value=write_only_current_value)
      end do
    end if

    if (this%all_scalar_integer_has_data) then
      do i=1,size(this%all_scalar_integer)
        call this%all_scalar_integer(i)%dumpForRestart(file_name, trim(adjustl(this%field_absolute_path%getString())), &
             only_current_value=write_only_current_value)
      end do
    end if

    if (this%all_real_vector_has_data) then
      do i=1,size(this%all_real_vector)
        call this%all_real_vector(i)%dumpForRestart(file_name, trim(adjustl(this%field_absolute_path%getString())), &
             only_current_value=write_only_current_value)
      end do
    end if

    if (this%all_vector5_has_data) then
      do i=1,size(this%all_vector5)
        call this%all_vector5(i)%dumpForRestart(file_name, trim(adjustl(this%field_absolute_path%getString())), &
             only_current_value=write_only_current_value)
      end do
    end if

    if (this%all_vector6_has_data) then
      do i=1,size(this%all_vector6)
        call this%all_vector6(i)%dumpForRestart(file_name, trim(adjustl(this%field_absolute_path%getString())), &
             only_current_value=write_only_current_value)
      end do
    end if

    if (this%all_generic_vector_has_data) then
      do i=1,size(this%all_generic_vector)
        call this%all_generic_vector(i)%dumpForRestart(file_name, trim(adjustl(this%field_absolute_path%getString())), &
             only_current_value=write_only_current_value)
      end do
    end if

    if (this%all_ss_scalar_has_data) then
      do i=1,size(this%all_ss_scalar)
        call this%all_ss_scalar(i)%dumpForRestart(file_name, trim(adjustl(this%field_absolute_path%getString())), &
             only_current_value=write_only_current_value)
      end do
    end if

    if (this%all_tensor2_has_data) then
      do i=1,size(this%all_tensor2)
        call this%all_tensor2(i)%dumpForRestart(file_name, trim(adjustl(this%field_absolute_path%getString())), &
             only_current_value=write_only_current_value)
      end do
    end if

    if (this%all_matrix66_has_data) then
      do i=1,size(this%all_matrix66)
        call this%all_matrix66(i)%dumpForRestart(file_name, trim(adjustl(this%field_absolute_path%getString())), &
             only_current_value=write_only_current_value)
      end do
    end if

    if (this%all_ss_generic_vector_has_data) then
      do i=1,size(this%all_ss_generic_vector)
        call this%all_ss_generic_vector(i)%dumpForRestart(file_name, trim(adjustl(this%field_absolute_path%getString())), &
             only_current_value=write_only_current_value)
      end do
    end if

    if (this%all_ss_vector5_has_data) then
      do i=1,size(this%all_ss_vector5)
        call this%all_ss_vector5(i)%dumpForRestart(file_name, trim(adjustl(this%field_absolute_path%getString())), &
             only_current_value=write_only_current_value)
      end do
    end if

    if (this%all_ss_vector6_has_data) then
      do i=1,size(this%all_ss_vector6)
        call this%all_ss_vector6(i)%dumpForRestart(file_name, trim(adjustl(this%field_absolute_path%getString())), &
             only_current_value=write_only_current_value)
      end do
    end if

    if (this%all_ss_matrix_has_data) then
      do i=1,size(this%all_ss_generic_matrix)
        call this%all_ss_generic_matrix(i)%dumpForRestart(file_name, trim(adjustl(this%field_absolute_path%getString())), &
             only_current_value=write_only_current_value)
      end do
    end if
    
    if (this%all_generic_matrix_has_data) then
      do i=1,size(this%all_generic_matrix)
        call this%all_generic_matrix(i)%dumpForRestart(file_name, trim(adjustl(this%field_absolute_path%getString())), &
             only_current_value=write_only_current_value)
      end do
    end if


    ! !write xmf
    if (this%phase_id.gt.0) then
      call suffix%setString("_ph"//trim(adjustl(int2string(this%phase_id))))
    else 
      call suffix%resetString()
    endif

    if (this%all_scalar_has_data) then
      do i=1,size(this%all_scalar)
        call this%all_scalar(i)%writeXDMFAttribute(xmf_writer, file_name, trim(adjustl(this%field_absolute_path%getString())), &
        write_only_current=write_only_current_value, suffix = trim(adjustl(suffix%getString())))
      end do
    end if

    if (this%all_scalar_integer_has_data) then
      do i=1,size(this%all_scalar_integer)
        call this%all_scalar_integer(i)%writeXDMFAttribute(xmf_writer, file_name, trim(adjustl(this%field_absolute_path%getString())), &
        write_only_current=write_only_current_value, suffix = trim(adjustl(suffix%getString())))
      end do
    end if

    if (this%all_real_vector_has_data) then
      do i=1,size(this%all_real_vector)
        call this%all_real_vector(i)%writeXDMFAttribute(xmf_writer, file_name, trim(adjustl(this%field_absolute_path%getString())), &
        write_only_current=write_only_current_value, suffix = trim(adjustl(suffix%getString())))
      end do
    end if

    if (this%all_vector5_has_data) then
      do i=1,size(this%all_vector5)
        call this%all_vector5(i)%writeXDMFAttribute(xmf_writer, file_name, trim(adjustl(this%field_absolute_path%getString())), &
        write_only_current=write_only_current_value, suffix = trim(adjustl(suffix%getString())))
      end do
    end if

    if (this%all_vector6_has_data) then
      do i=1,size(this%all_vector6)
        call this%all_vector6(i)%writeXDMFAttribute(xmf_writer, file_name, trim(adjustl(this%field_absolute_path%getString())), &
        write_only_current=write_only_current_value, suffix = trim(adjustl(suffix%getString())))
      end do
    end if

    if (this%all_generic_vector_has_data) then
      do i=1,size(this%all_generic_vector)
        call this%all_generic_vector(i)%writeXDMFAttribute(xmf_writer, file_name, trim(adjustl(this%field_absolute_path%getString())), &
        write_only_current=write_only_current_value, suffix = trim(adjustl(suffix%getString())))
      end do
    end if

    if (this%all_ss_scalar_has_data) then
      do i=1,size(this%all_ss_scalar)
        call this%all_ss_scalar(i)%writeXDMFAttribute(xmf_writer, file_name, trim(adjustl(this%field_absolute_path%getString())), &
        write_only_current=write_only_current_value, suffix = trim(adjustl(suffix%getString())))
      end do
    end if

    if (this%all_tensor2_has_data) then
      do i=1,size(this%all_tensor2)
        call this%all_tensor2(i)%writeXDMFAttribute(xmf_writer, file_name, trim(adjustl(this%field_absolute_path%getString())), &
        write_only_current=write_only_current_value, suffix = trim(adjustl(suffix%getString())))
      end do
    end if

    if (this%all_matrix66_has_data) then
      do i=1,size(this%all_matrix66)
        call this%all_matrix66(i)%writeXDMFAttribute(xmf_writer, file_name, trim(adjustl(this%field_absolute_path%getString())), &
        write_only_current=write_only_current_value, suffix = trim(adjustl(suffix%getString())))
      end do
    end if

    if (this%all_ss_generic_vector_has_data) then
      do i=1,size(this%all_ss_generic_vector)
        call this%all_ss_generic_vector(i)%writeXDMFAttribute(xmf_writer, file_name, trim(adjustl(this%field_absolute_path%getString())), &
        write_only_current=write_only_current_value, suffix = trim(adjustl(suffix%getString())))
      end do
    end if

    if (this%all_ss_vector5_has_data) then
      do i=1,size(this%all_ss_vector5)
        call this%all_ss_vector5(i)%writeXDMFAttribute(xmf_writer, file_name, trim(adjustl(this%field_absolute_path%getString())), &
        write_only_current=write_only_current_value, suffix = trim(adjustl(suffix%getString())))
      end do
    end if

    if (this%all_ss_vector6_has_data) then
      do i=1,size(this%all_ss_vector6)
        call this%all_ss_vector6(i)%writeXDMFAttribute(xmf_writer, file_name, trim(adjustl(this%field_absolute_path%getString())), &
        write_only_current=write_only_current_value, suffix = trim(adjustl(suffix%getString())))
      end do
    end if

    if (this%all_ss_matrix_has_data) then
      do i=1,size(this%all_ss_generic_matrix)
        call this%all_ss_generic_matrix(i)%writeXDMFAttribute(xmf_writer, file_name, trim(adjustl(this%field_absolute_path%getString())), &
        write_only_current=write_only_current_value, suffix = trim(adjustl(suffix%getString())))
      end do
    end if

    if (this%all_generic_matrix_has_data) then
      do i=1,size(this%all_generic_matrix)
        call this%all_generic_matrix(i)%writeXDMFAttribute(xmf_writer, file_name, trim(adjustl(this%field_absolute_path%getString())), &
        write_only_current=write_only_current_value, suffix = trim(adjustl(suffix%getString())))
      end do
    end if


  end subroutine


    module subroutine DumpMaterialPointValuesToTextFile(this, ix,iy,iz, file_id)
      implicit none
      class(all_grid_data), intent(inout) :: this
      integer, intent(in) :: ix,iy,iz, file_id
      integer :: i

      if (this%all_scalar_has_data) then
        do i=1,size(this%all_scalar)
          call this%all_scalar(i)%writePointToFile(ix,iy,iz, file_id)
        end do
      end if

      if (this%all_scalar_integer_has_data) then
        do i=1,size(this%all_scalar_integer)
          call this%all_scalar_integer(i)%writePointToFile(ix,iy,iz, file_id)
        end do
      end if

      if (this%all_real_vector_has_data) then
        do i=1,size(this%all_real_vector)
          call this%all_real_vector(i)%writePointToFile(ix,iy,iz, file_id)
        end do
      end if

      if (this%all_vector5_has_data) then
        do i=1,size(this%all_vector5)
          call this%all_vector5(i)%writePointToFile(ix,iy,iz, file_id)
        end do
      end if

      if (this%all_vector6_has_data) then
        do i=1,size(this%all_vector6)
          call this%all_vector6(i)%writePointToFile(ix,iy,iz, file_id)
        end do
      end if

      if (this%all_generic_vector_has_data) then
        do i=1,size(this%all_generic_vector)
          call this%all_generic_vector(i)%writePointToFile(ix,iy,iz, file_id)
        end do
      end if

      if (this%all_ss_scalar_has_data) then
        do i=1,size(this%all_ss_scalar)
          call this%all_ss_scalar(i)%writePointToFile(ix,iy,iz, file_id)
        end do
      end if

      if (this%all_tensor2_has_data) then
        do i=1,size(this%all_tensor2)
          call this%all_tensor2(i)%writePointToFile(ix,iy,iz, file_id)
        end do
      end if

      if (this%all_matrix66_has_data) then
        do i=1,size(this%all_matrix66)
          call this%all_matrix66(i)%writePointToFile(ix,iy,iz, file_id)
        end do
      end if

      if (this%all_ss_generic_vector_has_data) then
        do i=1,size(this%all_ss_generic_vector)
          call this%all_ss_generic_vector(i)%writePointToFile(ix,iy,iz, file_id)
        end do
      end if

      if (this%all_ss_vector5_has_data) then
        do i=1,size(this%all_ss_vector5)
          call this%all_ss_vector5(i)%writePointToFile(ix,iy,iz, file_id)
        end do
      end if

      if (this%all_ss_vector6_has_data) then
        do i=1,size(this%all_ss_vector6)
          call this%all_ss_vector6(i)%writePointToFile(ix,iy,iz, file_id)
        end do
      end if

      if (this%all_ss_matrix_has_data) then
        do i=1,size(this%all_ss_generic_matrix)
          call this%all_ss_generic_matrix(i)%writePointToFile(ix,iy,iz, file_id)
        end do
      end if

      if (this%all_generic_matrix_has_data) then
        do i=1,size(this%all_generic_matrix)
          call this%all_generic_matrix(i)%writePointToFile(ix,iy,iz, file_id)
        end do
      end if

    end subroutine


  module subroutine AGDReloadFromDump(this, file_name, time, dt, step, XLSEC, XLSEC_old)
    use hdf5
    use read_write_parallel_h5_mod, only : readRealAttribute, readIntegerAttribute, readRealMatrixAttribute
    use mpi, only : MPI_COMM_WORLD, MPI_INFO_NULL
    use mpi_useful_routines_mod, only : MPIBroadcastScalar, MPIBroadcastScalarInteger, &
                                        MPIBroadcastMatrix, MPIBarrier, AmIMPIMaster
    implicit none
    class(all_grid_data), intent(inout) :: this
    character(len=*), intent(in) :: file_name
    real(k_real), intent(out) :: time, dt
    integer, intent(out) :: step
    real(k_real), intent(out) ::  XLSEC(6,6), XLSEC_old(6,6)
    integer :: ierr, i
    INTEGER(HID_T) :: plist_id1, file_id

    CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id1, ierr)
    CALL h5pset_fapl_mpio_f(plist_id1, MPI_COMM_WORLD, MPI_INFO_NULL, ierr)
    CALL h5fopen_f(file_name, H5F_ACC_RDONLY_F, file_id, ierr, access_prp = plist_id1)

    if (this%all_scalar_has_data) then
      do i=1,size(this%all_scalar)
        call this%all_scalar(i)%reloadFromDump(file_id)
      end do
    end if

    if (this%all_scalar_integer_has_data) then
      do i=1,size(this%all_scalar_integer)
        call this%all_scalar_integer(i)%reloadFromDump(file_id)
      end do
    end if

    if (this%all_real_vector_has_data) then
      do i=1,size(this%all_real_vector)
        call this%all_real_vector(i)%reloadFromDump(file_id)
      end do
    end if

    if (this%all_vector5_has_data) then
      do i=1,size(this%all_vector5)
        call this%all_vector5(i)%reloadFromDump(file_id)
      end do
    end if

    if (this%all_vector6_has_data) then
      do i=1,size(this%all_vector6)
        call this%all_vector6(i)%reloadFromDump(file_id)
      end do
    end if

    if (this%all_generic_vector_has_data) then
      do i=1,size(this%all_generic_vector)
        call this%all_generic_vector(i)%reloadFromDump(file_id)
      end do
    end if

    if (this%all_ss_scalar_has_data) then
      do i=1,size(this%all_ss_scalar)
        call this%all_ss_scalar(i)%reloadFromDump(file_id)
      end do
    end if

    if (this%all_tensor2_has_data) then
      do i=1,size(this%all_tensor2)
        call this%all_tensor2(i)%reloadFromDump(file_id)
      end do
    end if

    if (this%all_matrix66_has_data) then
      do i=1,size(this%all_matrix66)
        call this%all_matrix66(i)%reloadFromDump(file_id)
      end do
    end if

    if (this%all_ss_generic_vector_has_data) then
      do i=1,size(this%all_ss_generic_vector)
        call this%all_ss_generic_vector(i)%reloadFromDump(file_id)
      end do
    end if

    if (this%all_ss_vector5_has_data) then
      do i=1,size(this%all_ss_vector5)
        call this%all_ss_vector5(i)%reloadFromDump(file_id)
      end do
    end if

    if (this%all_ss_vector6_has_data) then
      do i=1,size(this%all_ss_vector6)
        call this%all_ss_vector6(i)%reloadFromDump(file_id)
      end do
    end if

    if (this%all_ss_matrix_has_data) then
      do i=1,size(this%all_ss_generic_matrix)
        call this%all_ss_generic_matrix(i)%reloadFromDump(file_id)
      end do
    end if

    if (this%all_generic_matrix_has_data) then
      do i=1,size(this%all_generic_matrix)
        call this%all_generic_matrix(i)%reloadFromDump(file_id)
      end do
    end if

    CALL h5pclose_f(plist_id1, ierr)
    CALL h5fclose_f(file_id, ierr)
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    if (AmIMPIMaster()) then
      call H5fopen_f (file_name, H5F_ACC_RDWR_F, file_id, ierr)
      call readRealAttribute(file_id, "time", time)
      call readRealAttribute(file_id, "dt", dt)
      call readIntegerAttribute(file_id, "step", step)
      call readRealMatrixAttribute(file_id, "XLSEC", XLSEC)
      call readRealMatrixAttribute(file_id, "XLSEC_old", XLSEC_old)
      CALL h5fclose_f(file_id, ierr)
    endif

    call MPIBarrier()
    call MPIBroadcastScalar(time)
    call MPIBroadcastScalar(dt)
    call MPIBroadcastScalarInteger(step)
    call MPIBroadcastMatrix(XLSEC)
    call MPIBroadcastMatrix(XLSEC_old)
  end subroutine

  module subroutine AGDupdateAllStatefulVariables(this)
    implicit none
    class(all_grid_data), intent(inout) :: this
    integer :: i

    if (this%all_scalar_has_data) then
      do i=1,size(this%all_scalar)
        call this%all_scalar(i)%updateHistory()
      end do
    end if

    if (this%all_scalar_integer_has_data) then
      do i=1,size(this%all_scalar_integer)
        call this%all_scalar_integer(i)%updateHistory()
      end do
    end if

    if (this%all_real_vector_has_data) then
      do i=1,size(this%all_real_vector)
        call this%all_real_vector(i)%updateHistory()
      end do
    end if

    if (this%all_vector5_has_data) then
      do i=1,size(this%all_vector5)
        call this%all_vector5(i)%updateHistory()
      end do
    end if

    if (this%all_vector6_has_data) then
      do i=1,size(this%all_vector6)
        call this%all_vector6(i)%updateHistory()
      end do
    end if

    if (this%all_generic_vector_has_data) then
      do i=1,size(this%all_generic_vector)
        call this%all_generic_vector(i)%updateHistory()
      end do
    end if

    if (this%all_ss_scalar_has_data) then
      do i=1,size(this%all_ss_scalar)
        call this%all_ss_scalar(i)%updateHistory()
      end do
    end if

    if (this%all_tensor2_has_data) then
      do i=1,size(this%all_tensor2)
        call this%all_tensor2(i)%updateHistory()
      end do
    end if

    if (this%all_matrix66_has_data) then
      do i=1,size(this%all_matrix66)
        call this%all_matrix66(i)%updateHistory()
      end do
    end if

    if (this%all_ss_generic_vector_has_data) then
      do i=1,size(this%all_ss_generic_vector)
        call this%all_ss_generic_vector(i)%updateHistory()
      end do
    end if

    if (this%all_ss_vector5_has_data) then
      do i=1,size(this%all_ss_vector5)
        call this%all_ss_vector5(i)%updateHistory()
      end do
    end if

    if (this%all_ss_vector6_has_data) then
      do i=1,size(this%all_ss_vector6)
        call this%all_ss_vector6(i)%updateHistory()
      end do
    end if

    if (this%all_ss_matrix_has_data) then
      do i=1,size(this%all_ss_generic_matrix)
        call this%all_ss_generic_matrix(i)%updateHistory()
      end do
    end if

    if (this%all_generic_matrix_has_data) then
      do i=1,size(this%all_generic_matrix)
        call this%all_generic_matrix(i)%updateHistory()
      end do
    end if

  end subroutine

  module subroutine AGDResetAllStatefulVariables(this)
    implicit none
    class(all_grid_data), intent(inout) :: this
    integer :: i

    if (this%all_scalar_has_data) then
      do i=1,size(this%all_scalar)
        call this%all_scalar(i)%ResetHistory()
      end do
    end if

    if (this%all_scalar_integer_has_data) then
      do i=1,size(this%all_scalar_integer)
        call this%all_scalar_integer(i)%ResetHistory()
      end do
    end if

    if (this%all_real_vector_has_data) then
      do i=1,size(this%all_real_vector)
        call this%all_real_vector(i)%ResetHistory()
      end do
    end if

    if (this%all_vector5_has_data) then
      do i=1,size(this%all_vector5)
        call this%all_vector5(i)%ResetHistory()
      end do
    end if

    if (this%all_vector6_has_data) then
      do i=1,size(this%all_vector6)
        call this%all_vector6(i)%ResetHistory()
      end do
    end if

    if (this%all_generic_vector_has_data) then
      do i=1,size(this%all_generic_vector)
        call this%all_generic_vector(i)%ResetHistory()
      end do
    end if

    if (this%all_ss_scalar_has_data) then
      do i=1,size(this%all_ss_scalar)
        call this%all_ss_scalar(i)%ResetHistory()
      end do
    end if

    if (this%all_tensor2_has_data) then
      do i=1,size(this%all_tensor2)
        call this%all_tensor2(i)%ResetHistory()
      end do
    end if

    if (this%all_matrix66_has_data) then
      do i=1,size(this%all_matrix66)
        call this%all_matrix66(i)%ResetHistory()
      end do
    end if

    if (this%all_ss_generic_vector_has_data) then
      do i=1,size(this%all_ss_generic_vector)
        call this%all_ss_generic_vector(i)%ResetHistory()
      end do
    end if

    if (this%all_ss_vector5_has_data) then
      do i=1,size(this%all_ss_vector5)
        call this%all_ss_vector5(i)%ResetHistory()
      end do
    end if

    if (this%all_ss_vector6_has_data) then
      do i=1,size(this%all_ss_vector6)
        call this%all_ss_vector6(i)%ResetHistory()
      end do
    end if

    if (this%all_ss_matrix_has_data) then
      do i=1,size(this%all_ss_generic_matrix)
        call this%all_ss_generic_matrix(i)%ResetHistory()
      end do
    end if

    if (this%all_generic_matrix_has_data) then
      do i=1,size(this%all_generic_matrix)
        call this%all_generic_matrix(i)%ResetHistory()
      end do
    end if

  end subroutine

  module subroutine AGDGetGridDataScalarPointerByType(this, scalar_type, generic_griddata_scalar_ptr, has_data)
      implicit none
      class(all_grid_data), intent(in) :: this
      integer(kind(grid_data_var_type_enum)), intent(in) :: scalar_type
      class(griddata_scalar), dimension(:), pointer, intent(out):: generic_griddata_scalar_ptr
      logical, intent(out) :: has_data

      select case (scalar_type)
      case (scalar)
        generic_griddata_scalar_ptr => this%all_scalar
        has_data = this%all_scalar_has_data
      case default
        write(*,*) "AGDGetGridDataScalarPointerByType: Type not found!"
        stop
      end select

  end subroutine

  module subroutine AGDGetScalarDataPointerByNameTypeAndStatefulIndex(this, scalar_name, scalar_type, stateful_index, pointer_2_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: scalar_name
    integer(kind(grid_data_var_type_enum)), intent(in) :: scalar_type
    integer, intent(in) :: stateful_index
    real(k_real), dimension(:,:,:), pointer, intent(inout) :: pointer_2_data
    class(griddata_scalar), dimension(:), pointer :: generic_scalar_ptr => null()
    integer :: i
    logical :: has_data

    if (associated(pointer_2_data)) then
      write(*,*) "the Scalar pointer you want to link with the variable named ", scalar_name, " is already associated "
      write(*,*) "This is forbidden to prevent unwanted mistakes!"
      stop
    endif

    call this%AGDGetGridDataScalarPointerByType(scalar_type, generic_scalar_ptr, has_data)

    if (has_data) then
      do i=1,size(generic_scalar_ptr)
        if (generic_scalar_ptr(i)%name == scalar_name) then
          call generic_scalar_ptr(i)%getDataPointerByStatefulLevelIndex(stateful_index, pointer_2_data)
          return
        end if
      enddo
      write(*,*) "AGDGetScalarDataPointerByNameTypeAndStatefulIndex: can't find the variable named ", scalar_name
      stop
    else
      write(*,*) "AGDGetScalarDataPointerByNameTypeAndStatefulIndex: can't find the variable named ", scalar_name
      write(*,*) "AGDGetScalarDataPointerByNameTypeAndStatefulIndex: the scalar container of type ", scalar_type, "doesn't have data"
      stop
    endif
  end subroutine

  module subroutine AGDGetAvgScalarDataPointerByNameTypeAndStatefulIndex(this, scalar_name, scalar_type, stateful_index, pointer_2_avg_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: scalar_name
    integer(kind(grid_data_var_type_enum)), intent(in) :: scalar_type
    integer, intent(in) :: stateful_index
    real(k_real), pointer, intent(inout) :: pointer_2_avg_data
    class(griddata_scalar), dimension(:), pointer :: generic_scalar_ptr => null()
    integer :: i
    logical :: has_data

    if (associated(pointer_2_avg_data)) then
      write(*,*) "the AvgScalar pointer you want to link with the variable named ", scalar_name, " is already associated "
      write(*,*) "This is forbidden to prevent unwanted mistakes!"
      stop
    endif

    call this%AGDGetGridDataScalarPointerByType(scalar_type, generic_scalar_ptr, has_data)

    if (has_data) then
      do i=1,size(generic_scalar_ptr)
        if (generic_scalar_ptr(i)%name == scalar_name) then
          call generic_scalar_ptr(i)%getAvgDataPointerByStatefulLevelIndex(stateful_index, pointer_2_avg_data)
          return
        end if
      enddo
      write(*,*) "AGDGetAvgScalarDataPointerByNameTypeAndStatefulIndex: can't find the variable named ", scalar_name
      stop
    else
      write(*,*) "AGDGetAvgScalarDataPointerByNameTypeAndStatefulIndex: can't find the variable named ", scalar_name
      write(*,*) "AGDGetAvgScalarDataPointerByNameTypeAndStatefulIndex: the scalar container of type ", scalar_type, "doesn't have data"
      stop
    endif
  end subroutine

  module subroutine AGDGetGridDataScalarVariablePointerByName(this, scalar_name, scalar_type, pointer_2_variable)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: scalar_name
    integer(kind(grid_data_var_type_enum)), intent(in) :: scalar_type
    class(griddata_scalar), pointer, intent(inout) :: pointer_2_variable
    class(griddata_scalar), dimension(:), pointer :: generic_griddata_scalar_ptr  => null()
    integer :: i
    logical :: has_data

    nullify(generic_griddata_scalar_ptr)

    if (associated(pointer_2_variable)) then
      write(*,*) "the Scalar pointer you want to link with the variable named ", scalar_name, " is already associated "
      write(*,*) "This is forbidden to prevent unwanted mistakes!"
      stop
    endif

    call this%AGDGetGridDataScalarPointerByType(scalar_type, generic_griddata_scalar_ptr, has_data)

    if (has_data) then
      do i=1,size(generic_griddata_scalar_ptr)
        if (generic_griddata_scalar_ptr(i)%name == scalar_name) then
          pointer_2_variable => generic_griddata_scalar_ptr(i)
          return
        end if
      enddo
      write(*,*) "AGDGetGridDataScalarPointerByName: can't find the variable named ", scalar_name
      stop
    else
      write(*,*) "AGDGetGridDataScalarPointerByName: can't find the variable named ", scalar_name
      write(*,*) "AGDGetGridDataScalarPointerByName: the scalar container of type ", scalar_type, "doesn't have data"
      stop
    endif

  end subroutine

  module subroutine AGDGetScalarIntegerDataPointerByNameTypeAndStatefulIndex(this, scalar_integer_name, scalar_integer_type, stateful_index, pointer_2_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: scalar_integer_name
    integer(kind(grid_data_var_type_enum)), intent(in) :: scalar_integer_type
    integer, intent(in) :: stateful_index
    integer(k_int), dimension(:,:,:), pointer, intent(inout) :: pointer_2_data
    class(griddata_scalar_integer), dimension(:), pointer :: generic_scalar_integer_ptr
    integer :: i
    logical :: has_data

    if (associated(pointer_2_data)) then
      write(*,*) "the ScalarInteger pointer you want to link with the variable named ", scalar_integer_name, " is already associated "
      write(*,*) "This is forbidden to prevent unwanted mistakes!"
      stop
    endif

    select case (scalar_integer_type)
    case (scalar_integer)
      generic_scalar_integer_ptr => this%all_scalar_integer
      has_data = this%all_scalar_integer_has_data
    case default
      write(*,*) "Using the wrong function! "
      stop
    end select

    if (has_data) then
      do i=1,size(generic_scalar_integer_ptr)
        if (generic_scalar_integer_ptr(i)%name == scalar_integer_name) then
          call generic_scalar_integer_ptr(i)%getDataPointerByStatefulLevelIndex(stateful_index, pointer_2_data)
          return
        end if
      enddo
      write(*,*) "AGDGetScalarIntegerDataPointerByNameTypeAndStatefulIndex: can't find the variable named ", scalar_integer_name
      stop
    else
      write(*,*) "associated(this%all_scalar_integer) ", associated(this%all_scalar_integer)
      write(*,*) "size(this%all_scalar_integer) ", size(this%all_scalar_integer)
      write(*,*) "scalar name ", scalar_integer_name
      write(*,*) "AGDGetScalarIntegerDataPointerByNameTypeAndStatefulIndex: the scalar container of type ", scalar_integer_type, "doesn't have data"
      stop
    endif
  end subroutine

  module subroutine AGDGetGridDataVectorPointerByType(this, vector_type, generic_griddata_vector_ptr, has_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    integer(kind(grid_data_var_type_enum)), intent(in) :: vector_type
    class(griddata_vector), dimension(:), pointer, intent(out):: generic_griddata_vector_ptr
    logical, intent(out) :: has_data

    select case (vector_type)
    case (real_vector)
      generic_griddata_vector_ptr => this%all_real_vector
      has_data = this%all_real_vector_has_data
    case (vector5)
      generic_griddata_vector_ptr => this%all_vector5
      has_data = this%all_vector5_has_data
    case (vector6)
      generic_griddata_vector_ptr => this%all_vector6
      has_data = this%all_vector6_has_data
    case (generic_vector)
      generic_griddata_vector_ptr => this%all_generic_vector
      has_data = this%all_generic_vector_has_data
    case (ss_scalar)
      generic_griddata_vector_ptr => this%all_ss_scalar
      has_data = this%all_ss_scalar_has_data
    case default
      write(*,*) "AGDGetGridDataVectorPointerByType: Type not found!"
      stop
    end select

  end subroutine

  module subroutine AGDGetVectorDataPointerByNameTypeAndStatefulIndex(this, vector_name, vector_type, stateful_index, pointer_2_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: vector_name
    integer(kind(grid_data_var_type_enum)), intent(in) :: vector_type
    integer, intent(in) :: stateful_index
    real(k_real), dimension(:,:,:,:), pointer, intent(inout) :: pointer_2_data
    class(griddata_vector), dimension(:), pointer :: generic_griddata_vector_ptr => null()
    integer :: i
    logical :: has_data

    if (associated(pointer_2_data)) then
      write(*,*) "the Vector pointer you want to link with the variable named ", vector_name, " is already associated "
      write(*,*) "This is forbidden to prevent unwanted mistakes!"
      stop
    endif

    call this%AGDGetGridDataVectorPointerByType(vector_type, generic_griddata_vector_ptr, has_data)

    if (has_data) then
      do i=1,size(generic_griddata_vector_ptr)
        if (generic_griddata_vector_ptr(i)%name == vector_name) then
          call generic_griddata_vector_ptr(i)%getDataPointerByStatefulLevelIndex(stateful_index, pointer_2_data)
          return
        end if
      enddo
      write(*,*) "AGDGetVectorDataPointerByNameTypeAndStatefulIndex: can't find the variable named ", vector_name
      stop
    else
      write(*,*) "AGDGetVectorDataPointerByNameTypeAndStatefulIndex: can't find the variable named ", vector_name
      write(*,*) "AGDGetVectorDataPointerByNameTypeAndStatefulIndex: the vector container of type ", vector_type, "doesn't have data"
      stop
    endif
  end subroutine

  module subroutine AGDGetAvgVectorDataPointerByNameTypeAndStatefulIndex(this, vector_name, vector_type, stateful_index, pointer_2_avg_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: vector_name
    integer(kind(grid_data_var_type_enum)), intent(in) :: vector_type
    integer, intent(in) :: stateful_index
    real(k_real), dimension(:), pointer, intent(inout) :: pointer_2_avg_data
    class(griddata_vector), dimension(:), pointer :: generic_griddata_vector_ptr => null()
    integer :: i
    logical :: has_data

    if (associated(pointer_2_avg_data)) then
      write(*,*) "the AvgVector pointer you want to link with the variable named ", vector_name, " is already associated "
      write(*,*) "This is forbidden to prevent unwanted mistakes!"
      stop
    endif

    call this%AGDGetGridDataVectorPointerByType(vector_type, generic_griddata_vector_ptr, has_data)

    if (has_data) then
      do i=1,size(generic_griddata_vector_ptr)
        if (generic_griddata_vector_ptr(i)%name == vector_name) then
          call generic_griddata_vector_ptr(i)%getAvgDataPointerByStatefulLevelIndex(stateful_index, pointer_2_avg_data)
          return
        end if
      enddo
      write(*,*) "AGDGetAvgVectorDataPointerByNameTypeAndStatefulIndex: can't find the variable named ", vector_name
      stop
    else
      write(*,*) "AGDGetAvgVectorDataPointerByNameTypeAndStatefulIndex: the vector container of type ", vector_type, "doesn't have data"
      stop
    endif
  end subroutine

  module subroutine AGDGetGridDataVectorVariablePointerByName(this, vector_name, vector_type, pointer_2_variable)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: vector_name
    integer(kind(grid_data_var_type_enum)), intent(in) :: vector_type
    class(griddata_vector), pointer, intent(inout) :: pointer_2_variable
    class(griddata_vector), dimension(:), pointer :: generic_griddata_vector_ptr  => null()
    integer :: i
    logical :: has_data

    nullify(generic_griddata_vector_ptr)

    if (associated(pointer_2_variable)) then
      write(*,*) "the Vector pointer you want to link with the variable named ", vector_name, " is already associated "
      write(*,*) "This is forbidden to prevent unwanted mistakes!"
      stop
    endif

    call this%AGDGetGridDataVectorPointerByType(vector_type, generic_griddata_vector_ptr, has_data)

    if (has_data) then
      do i=1,size(generic_griddata_vector_ptr)
        if (generic_griddata_vector_ptr(i)%name == vector_name) then
          pointer_2_variable => generic_griddata_vector_ptr(i)
          return
        end if
      enddo
      write(*,*) "AGDGetGridDataVectorVariablePointerByName: can't find the variable named ", vector_name
      stop
    else
      write(*,*) "AGDGetGridDataVectorVariablePointerByName: the vector container of type ", vector_type, "doesn't have data"
      stop
    endif

  end subroutine

  module subroutine AGDGetGridDataMatrixPointerByType(this, matrix_type, generic_griddata_matrix_ptr, has_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    integer(kind(grid_data_var_type_enum)), intent(in) :: matrix_type
    class(griddata_matrix), dimension(:), pointer, intent(out):: generic_griddata_matrix_ptr
    logical, intent(out) :: has_data

    select case (matrix_type)
    case (tensor2)
      generic_griddata_matrix_ptr => this%all_tensor2
      has_data = this%all_tensor2_has_data
    case (matrix66)
      generic_griddata_matrix_ptr => this%all_matrix66
      has_data = this%all_matrix66_has_data
    case (ss_generic_vector)
      generic_griddata_matrix_ptr => this%all_ss_generic_vector
      has_data = this%all_ss_generic_vector_has_data
    case (ss_vector5)
      generic_griddata_matrix_ptr => this%all_ss_vector5
      has_data = this%all_ss_vector5_has_data
    case (ss_vector6)
      generic_griddata_matrix_ptr => this%all_ss_vector6
      has_data = this%all_ss_vector6_has_data
    case (generic_matrix)
      generic_griddata_matrix_ptr => this%all_generic_matrix
      has_data = this%all_generic_matrix_has_data
    case default
      write(*,*) "getGridDataMatrixPointerByType: Type not found!"
      stop
    end select

  end subroutine

  module subroutine AGDGetMatrixDataPointerByNameTypeAndStatefulIndex(this, matrix_name, matrix_type, stateful_index, pointer_2_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: matrix_name
    integer(kind(grid_data_var_type_enum)), intent(in) :: matrix_type
    integer, intent(in) :: stateful_index
    real(k_real), dimension(:,:,:,:,:), pointer, intent(inout) :: pointer_2_data
    class(griddata_matrix), dimension(:), pointer :: generic_griddata_matrix_ptr => null()
    integer :: i
    logical :: has_data

    if (associated(pointer_2_data)) then
      write(*,*) "the Matrix pointer you want to link with the variable named ", matrix_name, " is already associated "
      write(*,*) "This is forbidden to prevent unwanted mistakes!"
      stop
    endif

    call this%AGDGetGridDataMatrixPointerByType(matrix_type, generic_griddata_matrix_ptr, has_data)

    if (has_data) then
      do i=1,size(generic_griddata_matrix_ptr)
        if (generic_griddata_matrix_ptr(i)%name == matrix_name) then
          call generic_griddata_matrix_ptr(i)%getDataPointerByStatefulLevelIndex(stateful_index, pointer_2_data)
          return
        end if
      enddo
      write(*,*) "AGDGetMatrixDataPointerByNameTypeAndStatefulIndex: can't find the variable named ", matrix_name
      stop
    else
      write(*,*) "AGDGetMatrixDataPointerByNameTypeAndStatefulIndex: can't find the variable named ", matrix_name
      write(*,*) "AGDGetMatrixDataPointerByNameTypeAndStatefulIndex: the matrix container of type ", matrix_type, "doesn't have data"
      stop
    endif
  end subroutine

  module subroutine AGDGetAvgMatrixDataPointerByNameTypeAndStatefulIndex(this, matrix_name, matrix_type, stateful_index, pointer_2_avg_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: matrix_name
    integer(kind(grid_data_var_type_enum)), intent(in) :: matrix_type
    integer, intent(in) :: stateful_index
    real(k_real), dimension(:,:), pointer, intent(inout) :: pointer_2_avg_data
    class(griddata_matrix), dimension(:), pointer :: generic_griddata_matrix_ptr  => null()
    integer :: i
    logical :: has_data

    if (associated(pointer_2_avg_data)) then
      write(*,*) "the AvgMatrix pointer you want to link with the variable named ", matrix_name, " is already associated "
      write(*,*) "This is forbidden to prevent unwanted mistakes!"
      stop
    endif

    call this%AGDGetGridDataMatrixPointerByType(matrix_type, generic_griddata_matrix_ptr, has_data)

    if (has_data) then
      do i=1,size(generic_griddata_matrix_ptr)
        if (generic_griddata_matrix_ptr(i)%name == matrix_name) then
          call generic_griddata_matrix_ptr(i)%getAvgDataPointerByStatefulLevelIndex(stateful_index, pointer_2_avg_data)
          return
        end if
      enddo
      write(*,*) "AGDGetAvgMatrixDataPointerByNameTypeAndStatefulIndex: can't find the variable named ", matrix_name
      stop
    else
      write(*,*) "AGDGetAvgMatrixDataPointerByNameTypeAndStatefulIndex: can't find the variable named ", matrix_name
      write(*,*) "AGDGetAvgMatrixDataPointerByNameTypeAndStatefulIndex: the matrix container of type ", matrix_type, "doesn't have data"
      stop
    endif
  end subroutine

  module subroutine AGDGetGridDataMatrixVariablePointerByName(this, matrix_name, matrix_type, pointer_2_variable)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: matrix_name
    integer(kind(grid_data_var_type_enum)), intent(in) :: matrix_type
    class(griddata_matrix), pointer, intent(inout) :: pointer_2_variable
    class(griddata_matrix), dimension(:), pointer :: generic_griddata_matrix_ptr  => null()
    integer :: i
    logical :: has_data

    nullify(generic_griddata_matrix_ptr)

    if (associated(pointer_2_variable)) then
      write(*,*) "the Matrix pointer you want to link with the variable named ", matrix_name, " is already associated "
      write(*,*) "This is forbidden to prevent unwanted mistakes!"
      stop
    endif

    call this%AGDGetGridDataMatrixPointerByType(matrix_type, generic_griddata_matrix_ptr, has_data)

    if (has_data) then
      do i=1,size(generic_griddata_matrix_ptr)
        if (generic_griddata_matrix_ptr(i)%name == matrix_name) then
          pointer_2_variable => generic_griddata_matrix_ptr(i)
          return
        end if
      enddo
      write(*,*) "AGDGetGridDataMatrixVariablePointerByName: can't find the variable named ", matrix_name
      stop
    else
      write(*,*) "AGDGetGridDataMatrixVariablePointerByName: the matrix container of type ", matrix_type, "doesn't have data"
      stop
    endif

  end subroutine

  module subroutine AGDGetGridDataSSMatrixPointerByType(this, ss_matrix_type, generic_griddata_ss_matrix_ptr, has_data)
      implicit none
      class(all_grid_data), intent(in) :: this
      integer(kind(grid_data_var_type_enum)), intent(in) :: ss_matrix_type
      class(griddata_ss_generic_matrix), dimension(:), pointer, intent(out):: generic_griddata_ss_matrix_ptr
      logical, intent(out) :: has_data

      select case (ss_matrix_type)
      case (ss_generic_matrix)
        generic_griddata_ss_matrix_ptr => this%all_ss_generic_matrix
        has_data = this%all_ss_matrix_has_data
      case default
        write(*,*) "GetGridDataSSMatrixPointerByType: Type not found!"
        stop
      end select

  end subroutine

  module subroutine AGDGetSSMatrixDataPointerByNameTypeAndStatefulIndex(this, ss_matrix_name, ss_matrix_type, stateful_index, pointer_2_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: ss_matrix_name
    integer(kind(grid_data_var_type_enum)), intent(in) :: ss_matrix_type
    integer, intent(in) :: stateful_index
    real(k_real), dimension(:,:,:,:,:,:), pointer, intent(inout) :: pointer_2_data
    class(griddata_ss_generic_matrix), dimension(:), pointer :: generic_griddata_ss_matrix_ptr => null()

    integer :: i
    logical :: has_data

    if (associated(pointer_2_data)) then
      write(*,*) "the SSMatrix pointer you want to link with the variable named ", ss_matrix_name, " is already associated "
      write(*,*) "This is forbidden to prevent unwanted mistakes!"
      stop
    endif

    call this%AGDGetGridDataSSMatrixPointerByType(ss_matrix_type, generic_griddata_ss_matrix_ptr, has_data)

    if (has_data) then
      do i=1,size(generic_griddata_ss_matrix_ptr)
        if (generic_griddata_ss_matrix_ptr(i)%name == ss_matrix_name) then
          call generic_griddata_ss_matrix_ptr(i)%getDataPointerByStatefulLevelIndex(stateful_index, pointer_2_data)
          return
        end if
      enddo
      write(*,*) "AGDGetSSMatrixDataPointerByNameTypeAndStatefulIndex: can't find the variable named ", ss_matrix_name
      stop
    else
      write(*,*) "AGDGetSSMatrixDataPointerByNameTypeAndStatefulIndex: the matrix container of type ", ss_matrix_type, "doesn't have data"
      stop
    endif
  end subroutine

  module subroutine AGDGetAvgSSMatrixDataPointerByNameTypeAndStatefulIndex(this, ss_matrix_name, ss_matrix_type, stateful_index, pointer_2_avg_data)
    implicit none
    class(all_grid_data), intent(in) :: this
    character(len=*), intent(in) :: ss_matrix_name
    integer(kind(grid_data_var_type_enum)), intent(in) :: ss_matrix_type
    integer, intent(in) :: stateful_index
    real(k_real), dimension(:,:,:), pointer, intent(inout) :: pointer_2_avg_data
    class(griddata_ss_generic_matrix), dimension(:), pointer :: generic_griddata_ss_matrix_ptr => null()
    integer :: i
    logical :: has_data

    if (associated(pointer_2_avg_data)) then
      write(*,*) "the AvgSSMatrix pointer you want to link with the variable named ", ss_matrix_name, " is already associated "
      write(*,*) "This is forbidden to prevent unwanted mistakes!"
      stop
    endif

    call this%AGDGetGridDataSSMatrixPointerByType(ss_matrix_type, generic_griddata_ss_matrix_ptr, has_data)

    if (has_data) then
      do i=1,size(generic_griddata_ss_matrix_ptr)
        if (generic_griddata_ss_matrix_ptr(i)%name == ss_matrix_name) then
          call generic_griddata_ss_matrix_ptr(i)%getAvgDataPointerByStatefulLevelIndex(stateful_index, pointer_2_avg_data)
          return
        end if
      enddo
      write(*,*) "AGDGetAvgSSMatrixDataPointerByNameTypeAndStatefulIndex: can't find the variable named ", ss_matrix_name
      stop
    else
      write(*,*) "AGDGetAvgSSMatrixDataPointerByNameTypeAndStatefulIndex: the matrix container of type ", ss_matrix_type, "doesn't have data hence we can't find the ,atrix named ", ss_matrix_name
      stop
    endif
  end subroutine

end submodule
