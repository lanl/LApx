module all_porosity_models_temp_variables
use porosity_base_mod

! in this module we decleare varaible used to initalize a particualar model
! the variables below are helper and shall not be used to do proper calcualtion

! POROSITY MODELS
class(porosity_base), pointer :: porosity_base_temp => null()

end module
