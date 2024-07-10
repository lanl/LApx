module log_file_mod
  use kinds
  implicit none
  logical :: write_detailed_log_to_screen = .TRUE.
  character(len=10000) :: log_message

contains
subroutine writeToLogFile(message, create_file)
  use mpi_variables_mod, only : i_am_mpi_master
  implicit none
  logical, intent(in), optional :: create_file
  logical:: create_file_
  integer ::  file_id
  character(len=*), intent(in) :: message
  character(len=31) :: f_name="simualtion_log_file.log"

  create_file_ = .FALSE.
  if (present(create_file)) create_file_ = create_file

  file_id = 37
  f_name="simulation_log_file.log"

  if (i_am_mpi_master) then
    if (create_file_) then
      open(unit = file_id, file=f_name, status='unknown' )
      write(file_id,fmt='(*(A))') "LOG FILE"
      close(file_id)
    else
      open(unit = file_id, file=f_name, position='append', status='old')
      write(file_id,fmt='(*(A))') adjustl(trim(message))
      write(*,*) adjustl(trim(message))
      close(file_id)
    endif
  endif

end subroutine

subroutine writeToScreen(message)
    use mpi_variables_mod, only : i_am_mpi_master
    implicit none
    character(len=*), intent(in) :: message

    if (i_am_mpi_master.and.write_detailed_log_to_screen) &
      write(*,*) message
end subroutine
end module
