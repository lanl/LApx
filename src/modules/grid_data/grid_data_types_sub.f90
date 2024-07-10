submodule(grid_data_types) grid_data_types_sub

contains
  module subroutine checkDimensionBase(this)
    class(abstract_grid_data_type), intent(in) :: this

    if (LEN_TRIM(this%name) == 0) then
      write(*,*) "a griddata variable has an empty name"
      stop
    end if
    if (this%nx <= 0) then
      write(*,*) "nx < 0 for griddata named ", this%name
      stop
    end if
    if (this%ny <= 0) then
      write(*,*) "ny < 0 for griddata named ", this%name
      stop
    end if
    if (this%nz <= 0) then
      write(*,*) "nz < 0 for griddata named ", this%name
      stop
    end if
  end subroutine

end submodule
