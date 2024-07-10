module mechanics_module_init
use grid_data_var_type_mod
use all_grid_data_mod
#include "macro_debug.fpp"

implicit none

contains
subroutine getMechanicsGridVariables(grid_data)
  type(all_grid_data), intent(inout) :: grid_data

  ! here we define all the variables we need for the mechanics module
  ! we add tehm directly to the input vector
  call grid_data%addVar("velocity_gradient", tensor2, stateful_level=2)               
  call grid_data%addVar("total_displacement_gradient", tensor2, stateful_level=2)     
  ! call grid_data%addVar("stress", tensor2, stateful_level=2)
  ! call grid_data%addVar("stress_rate", tensor2, stateful_level=2)
  ! call grid_data%addVar("elastic_strain", tensor2, stateful_level=2)                  
  ! call grid_data%addVar("elastic_strain_rate", tensor2, stateful_level=2)             
  ! call grid_data%addVar("plastic_strain", tensor2, stateful_level=2)  
  ! call grid_data%addVar("plastic_strain_rate", tensor2, stateful_level=2)       
  ! call grid_data%addVar("total_strain", tensor2, stateful_level=2)
  ! call grid_data%addVar("total_strain_rate", tensor2, stateful_level=2)
  ! call grid_data%addVar("stiffness", matrix66, stateful_level=2)   !->stiffness66


end subroutine

end module mechanics_module_init
