module print_utils_mod
use kinds
implicit none

public :: printToScreen
private :: printScalar, printMatrix, printMatrixInteger, printVector, &
           printVectorInteger, printTensor4

interface printToScreen
  procedure printMatrix
  procedure printVector
  procedure printTensor4
  procedure printScalar
  procedure printVectorInteger
  procedure printMatrixInteger
end interface

interface writeToString
  procedure writeToStringReal
  procedure writeToStringMatrixReal
  procedure writeToStringMatrixInteger
end interface

contains

  subroutine printScalar(scalar, name)
    use log_file_mod, only : write_detailed_log_to_screen
    use mpi_variables_mod, only : i_am_mpi_master
    implicit none
    real(kind=k_real), intent(in) :: scalar
    character(len=*), intent(in) :: name
    if (i_am_mpi_master.and.write_detailed_log_to_screen) then
      write(*,*) ""
      write(*,*) name, ": "
      write(*,*) scalar
    endif
  end subroutine

subroutine printMatrix(matrix, name)
  use log_file_mod, only : write_detailed_log_to_screen
  use mpi_variables_mod, only : i_am_mpi_master
  implicit none
  real(kind=k_real), intent(in), dimension(:,:) :: matrix
  character(len=*), intent(in) :: name
  integer :: dims(2), i ,j

  if (i_am_mpi_master.and.write_detailed_log_to_screen) then
    dims = shape(matrix)
    associate(ni=>dims(1), nj=>dims(2))
      write(*,*) ""
      write(*,*) name, ": "
      write(*,*) ( (matrix(i,j),j=1,nj), new_line("A"), i=1,ni )
    end associate
  endif
end subroutine

subroutine writeToStringReal(real_value, name, mystring)
  use number_to_string_mod
  implicit none
  real(kind=k_real), intent(in) :: real_value
  character(len=*), intent(in) :: name
  character(len=*), intent(out) :: mystring

    write(mystring,*) trim(adjustl(real2string(real_value))) // new_line("A")
    mystring = trim(adjustl(name)) // ": "//new_line("A") // trim(mystring)
end subroutine

subroutine writeToStringMatrixReal(matrix, name, mystring)
  use number_to_string_mod
  implicit none
  real(kind=k_real), intent(in), dimension(:,:) :: matrix
  character(len=*), intent(in) :: name
  character(len=*), intent(out) :: mystring
  integer :: dims(2), i ,j

  dims = shape(matrix)
  associate(ni=>dims(1), nj=>dims(2))
    write(mystring,*) ( (matrix(i,j),j=1,nj), new_line("A"), i=1,ni )
  end associate
  mystring = trim(adjustl(name)) // ": "//new_line("A") // trim(mystring)
end subroutine

subroutine writeToStringMatrixInteger(matrix, name, mystring)
  use number_to_string_mod
  implicit none
  integer, intent(in), dimension(:,:) :: matrix
  character(len=*), intent(in) :: name
  character(len=*), intent(out) :: mystring
  integer :: dims(2), i ,j

  dims = shape(matrix)
  associate(ni=>dims(1), nj=>dims(2))
    write(mystring,*) ( (matrix(i,j),j=1,nj), new_line("A"), i=1,ni )
  end associate
  mystring = trim(adjustl(name)) // ": "//new_line("A") // trim(mystring)

end subroutine


subroutine printMatrixInteger(matrix, name)
  use log_file_mod, only : write_detailed_log_to_screen
  use mpi_variables_mod, only : i_am_mpi_master
  implicit none
  integer, intent(in), dimension(:,:) :: matrix
  character(len=*), intent(in) :: name
  integer :: dims(2), i ,j

  if (i_am_mpi_master.and.write_detailed_log_to_screen) then
    dims = shape(matrix)
    associate(ni=>dims(1), nj=>dims(2))
      write(*,*) ""
      write(*,*) name, ": "
      write(*,*) ( (matrix(i,j),j=1,nj), new_line("A"), i=1,ni )
    end associate
  endif
end subroutine

subroutine printVector(vector, name)
  use log_file_mod, only : write_detailed_log_to_screen
  use mpi_variables_mod, only : i_am_mpi_master
  implicit none
  real(kind=k_real), intent(in), dimension(:) :: vector
  character(len=*), intent(in) :: name
  integer :: ni, i

  if (i_am_mpi_master.and.write_detailed_log_to_screen) then
    ni = size(vector)
    write(*,*) ""
    write(*,*) name, ": "
    write(*,*) ( vector(i), new_line("A"), i=1,ni )
  endif
end subroutine

subroutine printVectorInteger(vector, name)
  use log_file_mod, only : write_detailed_log_to_screen
  use mpi_variables_mod, only : i_am_mpi_master
  implicit none
  integer, intent(in), dimension(:) :: vector
  character(len=*), intent(in) :: name
  integer :: ni, i

  if (i_am_mpi_master.and.write_detailed_log_to_screen) then
    ni = size(vector)
    write(*,*) ""
    write(*,*) name, ": "
    write(*,*) ( vector(i), new_line("A"), i=1,ni )
  endif
end subroutine

subroutine printTensor4(t4, name)
  use log_file_mod, only : write_detailed_log_to_screen
  use mpi_variables_mod, only : i_am_mpi_master
  implicit none
  real(kind=k_real), intent(in), dimension(3,3,3,3) :: t4
  character(len=*), intent(in) :: name
  integer :: i ,j

  if (i_am_mpi_master.and.write_detailed_log_to_screen) then
    write(*,*) ""
    write(*,*) name, ": "
    do i=1,3
    do j=1,3
      write(*,*) "(:,:,",j,",",i,")"
      call printMatrix(t4(:,:,j,i), "")
      write(*,*)
    enddo
    enddo
  endif
end subroutine


end module
