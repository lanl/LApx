module number_to_string_mod
use kinds
implicit none

contains

  function int2string(myint) result(newString)
    character(len=10) :: newString
    integer, intent(in) :: myint
    write(newString,'(i10)') myint
  end function

  function real2string(myreal) result(newString)
    character(len=100) :: newString
    real(k_real), intent(in) :: myreal
    write(newString, fmt='(*(G0))') myreal
  end function

end module number_to_string_mod
