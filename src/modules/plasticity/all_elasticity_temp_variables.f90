module all_elasticity_temp_variables
use stiffness_base_mod

! in this module we decleare variables used to initalize the phase elastic stiffness
class(stiffness_isotropic), pointer :: stiffness_model_isotropic_temp => null()
class(stiffness_isotropic_shear), pointer :: stiffness_model_isotropic_shear_temp => null()
class(stiffness_cubic), pointer :: stiffness_model_cubic_temp => null()
class(stiffness_hexagonal), pointer :: stiffness_model_hexagonal_temp => null()

end module
