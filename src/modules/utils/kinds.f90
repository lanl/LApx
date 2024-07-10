module kinds
use, intrinsic :: ISO_FORTRAN_ENV
implicit none

! changing the defintion of of
! basic types used by the program for:
integer, parameter :: k_real = REAL64 !-> real and complex data
integer, parameter :: k_int = INT32   !-> integers, usch as phase, and grainid


! additional integer types for loops and other variables
integer, parameter :: k_shor_int = INT8 !-> max 256
integer, parameter :: dp = REAL64 !-> real and complex data
integer, parameter :: qp = REAL128 !-> real and complex data
contains

  function int2real(myint) result(newReal)
    real(k_real) :: newReal
    integer, intent(in) :: myint
    newReal = real(myint, k_real)
  end function

  function real2int(myReal) result(newInt)
    real(k_real), intent(in) :: myReal
    integer(k_int):: newInt

    newInt = int(myReal, k_int)
  end function

  function logical2Int(mylogical) result(newint)
    integer(k_int) :: newint
    logical, intent(in) :: mylogical
    newint = 0
    if (mylogical) newint = 1
  end function


SUBROUTINE stderr(message)
  ! "@(#) stderr writes a message to standard error using a standard f2003 method"
      USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : ERROR_UNIT ! access computing environment
      IMPLICIT NONE
      CHARACTER(LEN=*),INTENT(IN) :: message
      WRITE(ERROR_UNIT,'(a)')trim(message) ! write message to standard error
  END SUBROUTINE stderr

end module kinds
