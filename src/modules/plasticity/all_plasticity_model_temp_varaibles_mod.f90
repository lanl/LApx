module all_plasticity_models_temp_variables
use climb_mod
use glide_mod
use isotropic_diffusion_mod
use hardening_model_mod

! in this module we decleare varaible used to initalize a particualar model
! the variables below are helper and shall not be used to do proper calcualtion

! GLIDE MODELS
class(hutchinson_glide), pointer :: hutchinson_glide_temp  => null()
class(wkkct_glide), pointer :: wkkct_glide_temp  => null()

!CLIMB MODELS
class(transient_climb), pointer :: transient_climb_temp  => null()
class(exponential_climb), pointer :: exponential_climb_temp  => null()
class(wkkct_climb), pointer :: wkkct_climb_temp  => null()
class(wkkct_climb_irradiation), pointer :: wkkct_climb_irradiation_temp  => null()
! class(wkkct_climb_chemo_mech), pointer :: wkkct_climb_chemo_mech_temp

!ISOTROPIC DIFFUSION MODELS
! class(dev_stress_diffusion_gb), pointer :: dev_stress_diffusion_gb_temp  => null()
class(nabarro_herring_plus_coble_diffusion), pointer :: nabarro_herring_plus_coble_diffusion_temp  => null()

end module
