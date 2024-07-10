module smooth_function_mod
use kinds
use test_utils_mod

implicit none
real(k_real), parameter, private :: p_value = 30._k_real
real(k_real), parameter, private :: symm_smooth_p_value = -30._k_real
contains

subroutine pNorm(x, p, pNorm_val)
  real(k_real), intent(in) :: x(:), p
  real(k_real), intent(out) :: pNorm_val
  pNorm_val = (sum(x**p))**(1._k_real/p)
end subroutine

subroutine dPNormdX(x, p, dPnormdX_vals)
  real(k_real), intent(in) :: x(:), p
  real(k_real), intent(out) :: dPnormdX_vals(:)
  real(k_real) :: temp

  temp = sum(x**p)**(1._k_real/p - 1._k_real)
  dPnormdX_vals = x**(p-1._k_real)*temp

end subroutine

function smoothMax(x1, x2) result(s_max)
  real(k_real), intent(in) :: x1, x2
  real(k_real) :: x(2), s_max

  x(1) = x1
  x(2) = x2

  call pNorm(x, p_value, s_max)
end function

function dsmoothMax_dX(x1, x2) result(dp_norm_dx)
  real(k_real), intent(in) :: x1, x2
  real(k_real) :: x(2), dp_norm_dx(2)

  x(1) = x1
  x(2) = x2

  call dPNormdX(x, p_value, dp_norm_dx)
end function


real (k_real) function symmSmoothClamp(y, y_lim) result(y_clamped)
  real(k_real), intent(in) :: y, y_lim
  real(k_real), dimension(2) :: x

  x(1) = y
  x(2) = y_lim

  call pNorm(x, symm_smooth_p_value, y_clamped)
  y_clamped = y_clamped*sign(1._k_real, y)
end function

real (k_real) function dSymmSmoothClampDY(y, y_lim) result(dyclamped_dy)
  real(k_real), intent(in) :: y, y_lim
  real(k_real), dimension(2) :: x, dpnorm_dy

  x(1) = y
  x(2) = y_lim
  if (y.ne.0._k_real) then
    call dPNormdX(x, symm_smooth_p_value, dpnorm_dy)
    dyclamped_dy = dpnorm_dy(1)*sign(1._k_real, y)

#ifdef __DEBUG__
    call test_is_finite(dpnorm_dy, "dpnorm_dy")
    call test_is_finite(dyclamped_dy, "dyclamped_dy")
#endif

  else
    dyclamped_dy = 1._k_real
  endif
end function

end module
