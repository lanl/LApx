module units_conversion_mod
use kinds
use math_constants, only : PI
real(k_real), parameter :: eV2Joule = 1.602176634e-19
real(k_real), parameter :: Joule2eV = 1._k_real/eV2Joule
real(k_real), parameter :: MPA2Pa = 1e6
real(k_real), parameter :: Pa2MPa = 1e-6
real(k_real), parameter :: deg2rad = PI/180._k_real
real(k_real), parameter :: rad2deg = 180._k_real/PI

end module
