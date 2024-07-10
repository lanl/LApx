

!tensor2 quantitites
real(k_real), dimension(:,:,:,:,:), pointer :: vel_grad => null()
real(k_real), dimension(:,:,:,:,:), pointer :: vel_grad_old => null()
real(k_real), dimension(:,:,:,:,:), pointer :: disgradtot => null()
real(k_real), dimension(:,:,:,:,:), pointer :: disgradtot_old => null()

real(k_real), dimension(:,:,:,:,:), pointer :: grid_stress => null()
real(k_real), dimension(:,:,:,:,:), pointer :: grid_stress_old => null()
real(k_real), dimension(:,:,:,:,:), pointer :: grid_stress_rate => null()
real(k_real), dimension(:,:,:,:,:), pointer :: grid_stress_rate_old => null()


real(k_real), dimension(:,:,:,:,:), pointer :: eelfield_rate => null()
real(k_real), dimension(:,:,:,:,:), pointer :: eelfield => null()
real(k_real), dimension(:,:,:,:,:), pointer :: eelfield_old => null()

real(k_real), dimension(:,:,:,:,:), pointer :: edotp => null()
real(k_real), dimension(:,:,:,:,:), pointer :: edotp_old => null()
real(k_real), dimension(:,:,:,:,:), pointer :: ept => null()
real(k_real), dimension(:,:,:,:,:), pointer :: ept_old => null()
real(k_real), dimension(:,:,:,:,:), pointer :: edotp_eq => null() ! I'm not sure what this guy is or means
real(k_real), dimension(:,:,:,:,:), pointer :: total_strain => null()
real(k_real), dimension(:,:,:,:,:), pointer :: total_strain_old => null()
real(k_real), dimension(:,:,:,:,:), pointer :: total_strain_rate => null()
real(k_real), dimension(:,:,:,:,:), pointer :: total_strain_rate_old => null()

! matrix66 quantities
real(k_real), dimension(:,:,:,:,:), pointer :: stiffness66 => null()
real(k_real), dimension(:,:,:,:,:), pointer :: stiffness66_old => null()
