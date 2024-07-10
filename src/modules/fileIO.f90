#define MAXLINELENGTH 200
module fileIO
  use kinds
implicit none


contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!C ********************************************************************
!C     SUBROUTINE VPSC_INPUT      --->      VERSION 31/jan/99
!C
!C     READS CHARACTERISTICS OF THE RUN: # OF PHASES, NAMES OF INPUT FILES,
!C     DEFORMATION TO BE IMPOSED, CONVERGENCE PARAMETERS, ETC.
!C     READS SINGLE CRYSTAL PROPERTIES: DEFORMATION MODES, CRSS, HARDENING
!C     READS CRYSTAL AND MORPHOLOGIC TEXTURES.
!C     INITIALIZES ARRAYS REQUIRED TO RUN VPSC.
!C     OPENS AND CLOSES INPUT FILES.   OPENS OUTPUT FILES.
!C
!C     MODIFIED 21/07/98 by CNT:
!C     INITI34ALIZATION RELATED TO 'ELEMENTS' IS DONE INSIDE A SINGLE BLOCK.
!C *****************************************************************************
!C

SUBROUTINE INPUT
  USE global, only : all_mighty_grid_global, phase_material_array, microstrucutre_info, fft_grid_data_global, sim_all_macro_data, the_bc_object, &
  common_grid_data_global, voxel_size, wgt, initgriddata
  use read_microstructure_file_mod, only : read_microstructure_from_file
  ! use texture_mod, only : update_stiffness
  use tensor_math_mod, only: getSkewPart, getSymmetricPart, computeVMEquivalentStrainNonSymmetric
  use mpi_variables_mod, only : mpi_master_rank, mpi_rank, mpi_size
  use mpi_useful_routines_mod, only : MPISumMatrix, MPIBarrier
  use log_file_mod, only : writeToScreen, write_detailed_log_to_screen
  IMPLICIT NONE

  integer :: k, n_points_xyz(3)
  integer, pointer, dimension(:) :: n_ss => null(), n_gauss => null()
  ! let each process read the microstructure file in a sequential fashion


  ! read option file process by process in an ordered fashion
  do k = 0,mpi_size-1
    if (mpi_rank.eq.k) then
      call parse_options_file()
      write(*,*) "process ", mpi_rank, "finshed reading option file"
    endif
    call MPIBarrier()
  enddo
  call MPIBarrier()
  IF (write_detailed_log_to_screen) call writeToScreen("all process parsed the option file and the BC file")

  call microstrucutre_info%getGlobalDimensions(n_points_xyz)
  call all_mighty_grid_global%AMGSetDimension(n_points_xyz)
  call phase_material_array%linkBaseObjects(all_mighty_grid_global, sim_all_macro_data, the_bc_object)
  

    ! read option file process by process in an ordered fashion
  do k = 0,mpi_size-1
    if (mpi_rank.eq.k) then
      call phase_material_array%readFromFile(microstrucutre_info, n_ss, n_gauss )
      write(*,*) "process ", mpi_rank, "finshed reading material files"
    endif
    call MPIBarrier()
  enddo
  call MPIBarrier()
  IF (write_detailed_log_to_screen) call writeToScreen("all process parsed the option file and the BC file")

  call all_mighty_grid_global%AMGSetNumSlipSystems(n_ss)
  
  wgt = 1._k_real/int2real(all_mighty_grid_global%AMGgetNumVoxel())
  voxel_size = microstrucutre_info%getVoxelPhysicalDimension()
  ! we have the problem's dimension now we can initialize teh fftw arrays
  ! this will also provide the dimension of arrays for each mpi_rank
  ! call fft_grid_data_global%initFFTData(npts1, npts2, npts3, npts3_rank, npts3_offset, npts3_start, npts3_end)
  call all_mighty_grid_global%AMGgetFFTGridPtr(fft_grid_data_global)
  call all_mighty_grid_global%AMGGetCommonGridPtr(common_grid_data_global)
  call fft_grid_data_global%addFieldVariablesToGrid(common_grid_data_global)
  
  call phase_material_array%addFieldVariablesToGrid()
  ! finished popualting the  field variable list, now we can allocate stuff and set pointers
  call initGridData()
  ! grid allocated: read phase input files and store relevant information

  call MPIBarrier()
  IF (write_detailed_log_to_screen) call writeToScreen( "all processes succesfully initialized griddata" )

  do k = 0,mpi_size-1
    if (mpi_rank.eq.k) then
      write(*,*)  "this is the file", microstrucutre_info%getMicrostructureFileName()
      call read_microstructure_from_file(microstrucutre_info%getMicrostructureFileName(), all_mighty_grid_global)
      write(*,*) "process ", mpi_rank, "finshed reading microstructure file"
    endif
    call MPIBarrier()
  enddo
  CALL MPIBarrier()

  ! now that the grid is allocated we init all the grid pointers
  call fft_grid_data_global%initGridPointers()
  call fft_grid_data_global%initFrequencyVectorAndTensor(voxel_size)
  call sim_all_macro_data%sim_macro_field_averages%init(common_grid_data_global)
  call sim_all_macro_data%sim_macro_field_averages%initGridPointers()
  call the_bc_object%initGridPointers()
  call phase_material_array%initGridPointers()


end subroutine input

! subroutine data_crystal_andrea(fname, iph, phase_material, all_mighty_grid_in, n_ss, n_gauss)
!   use read_from_file_utils
!   use phase_material_mod, only : phase_material_type
!   use common_material_parameter_mod
!   use log_file_mod, only : write_detailed_log_to_screen
!   use embedded_slip_systems_mod, only : crystal_type_enum
!   use all_mighty_grid_mod, only : all_mighty_grid_type
!   use all_plasticity_models_temp_variables
!   use all_elasticity_temp_variables
!   use all_porosity_models_temp_variables
!   use mpi_useful_routines_mod, only : AmIMPIMaster
!   use global, only : elasticity_base_ptr, inelastic_strain_base_ptr, porosity_base_ptr, sim_all_macro_data, the_bc_object
!   integer, intent(in) :: iph
!   character(len=*), intent(in) :: fname
!   class(phase_material_type), pointer, intent(inout) :: phase_material
!   type(all_mighty_grid_type), pointer, intent(inout) :: all_mighty_grid_in
!   integer, intent(out) :: n_ss, n_gauss
!   type(string_array) :: dummy_string_array
!   type(string_type) :: material_model_name, elasticity_model_name
!   type(file_reader) :: matf_reader
!   integer :: phase_id
!   logical :: use_isotropic_plasticity, &
!              use_crystal_plasticity, &
!              use_glide, &
!              use_climb, &
!              use_diffusion, &
!              use_porosity, &
!              use_damage
!   ! variables used to read minimal CP data
!   real(k_real) :: ratio_edge_screw
!   integer, dimension(:), pointer :: n_ss_per_mode => null()
!   real(k_real), pointer, dimension(:) :: burgVectorL => null(), &
!                                          self_hardening => null(), &
!                                          latent_hardening => null()
!   real(k_real) :: latent_hardening_other_modes
!   class(common_material_parameter), pointer :: common_material_parameter_ptr => null()
!   integer(kind(crystal_type_enum)) :: crystal_type
!   ! integral formulation parameters
!   real(k_real) :: n_std_dev_gauss_integration
!   real(k_real), pointer, dimension(:,:) :: schmid_ca => null(), climb_ca => null()


!   n_ss = 1
!   n_gauss = 1
!   call matf_reader%openReadTextFile(fname)
!   allocate(common_material_parameter_ptr)

!   call matf_reader%readLineAndCheckStringsAreEqual("--Phase-Parameters", dummy_string_array)
!   call matf_reader%readParameter("phase-id", phase_id)
!   if (phase_id.ne.iph) then
!     if (AmIMPIMaster().and.write_detailed_log_to_screen) write(*,*) "I,m reading the phase file for phase ", iph, " but the provided file contains information for phase number ", phase_id
!     error stop " . Abort!"
!   end if

!   ! sanity checks
!   call matf_reader%readParameter("use-isotropic-plasticity", use_isotropic_plasticity)
!   call matf_reader%readParameter("use-crystal-plasticity", use_crystal_plasticity)
!   call matf_reader%readParameter("use-glide", use_glide)
!   call matf_reader%readParameter("use-climb", use_climb)
!   call matf_reader%readParameter("use-diffusion", use_diffusion)
!   call matf_reader%readParameter("use-porosity", use_porosity)
!   call matf_reader%readParameter("use-damage", use_damage)
!   call matf_reader%skipEmptyLine()

!   ! a few sanity checks before continuing
!   if (use_isotropic_plasticity.and.use_crystal_plasticity) error stop "you can't use both isotropic and crystal palsticity"
!   if (use_glide.and.(.not.use_crystal_plasticity)) error stop "you can't add glide without using crystal plasticity"
!   if (use_climb.and.(.not.use_crystal_plasticity)) error stop "you can't add diffusion without using crystal plasticity"
!   if (use_climb.and.(.not.use_glide)) error stop "you can't add climb without adding glide"
!   if (use_damage.and.(.not.(use_porosity))) error stop "you cannot add damage without adding porosity"

!   call matf_reader%readParameter("--Elasticity", elasticity_model_name)
!   select case(elasticity_model_name%getString())

!   case ("isotropic-linear")
!     allocate(stiffness_model_isotropic_temp)
!     elasticity_base_ptr => stiffness_model_isotropic_temp
!     call stiffness_model_isotropic_temp%linkBaseObjects(phase_id, all_mighty_grid_in, sim_all_macro_data, the_bc_object)
!     call stiffness_model_isotropic_temp%initParameters(phase_id, common_material_parameter_ptr)

!   case ("isotropic-shear-poisson-linear")
!     allocate(stiffness_model_isotropic_shear_temp)
!     elasticity_base_ptr => stiffness_model_isotropic_shear_temp
!     call stiffness_model_isotropic_shear_temp%linkBaseObjects(phase_id, all_mighty_grid_in, sim_all_macro_data, the_bc_object)
!     call stiffness_model_isotropic_shear_temp%initParameters(phase_id, common_material_parameter_ptr)

!   case ("cubic-linear")
!     allocate(stiffness_model_cubic_temp)
!     elasticity_base_ptr => stiffness_model_cubic_temp
!     call stiffness_model_cubic_temp%linkBaseObjects(phase_id, all_mighty_grid_in, sim_all_macro_data, the_bc_object)
!     call stiffness_model_cubic_temp%initParameters(phase_id, common_material_parameter_ptr)

!   case ("hexagonal-linear")
!     allocate(stiffness_model_hexagonal_temp)
!     elasticity_base_ptr => stiffness_model_hexagonal_temp
!     call stiffness_model_hexagonal_temp%linkBaseObjects(phase_id, all_mighty_grid_in, sim_all_macro_data, the_bc_object)
!     call stiffness_model_hexagonal_temp%initParameters(phase_id, common_material_parameter_ptr)

!   case default
!     if (AmIMPIMaster().and.write_detailed_log_to_screen) write(*,*) "unrecognized elasticity model name ", elasticity_model_name%getString()
!     if (AmIMPIMaster().and.write_detailed_log_to_screen) write(*,*) "available material models are: "
!     if (AmIMPIMaster().and.write_detailed_log_to_screen) write(*,*) "isotropic-linear, cubic-linear, and hexagonal-linear"
!     error stop "abort"
!   end select

!   call elasticity_base_ptr%readMaterialParametersFromFile(matf_reader)
!   call phase_material%setStiffness(elasticity_base_ptr)

!   select case(elasticity_model_name%getString())
!   case ("isotropic-linear")
!     nullify(stiffness_model_isotropic_temp)
!   case ("isotropic-shear-poisson-linear")
!     nullify(stiffness_model_isotropic_shear_temp)
!   case ("cubic-linear")
!     nullify(stiffness_model_cubic_temp)
!   case ("hexagonal-linear")
!     nullify(stiffness_model_hexagonal_temp)
!   case default
!     if (AmIMPIMaster().and.write_detailed_log_to_screen) write(*,*) "somthing is wrong when trying nullifying the pointer for ", elasticity_model_name%getString()
!     error stop "abort"
!   end select
!   nullify(elasticity_base_ptr)

!   if (use_isotropic_plasticity.or.use_crystal_plasticity) then
!    call matf_reader%skipEmptyLine()

!   call matf_reader%readLineAndCheckStringsAreEqual("--Common-Material-Parameters", dummy_string_array)
!   call common_material_parameter_ptr%readCommonMaterialParametersFromFile(matf_reader)
!   call matf_reader%skipEmptyLine()

!   !if we use crystal plasticity then
!   if (use_crystal_plasticity) then
!     call matf_reader%readLineAndCheckStringsAreEqual("--Crystal-Parameters", dummy_string_array)
!     ! call readCrystalParameters(matf_reader, &
!     !                                     ratio_edge_screw, &
!     !                                     n_ss_per_mode, &
!     !                                     burgVectorL, &
!     !                                     self_hardening, &
!     !                                     latent_hardening, &
!     !                                     latent_hardening_other_modes, &
!     !                                     n_gauss, &
!     !                                     n_std_dev_gauss_integration, &
!     !                                     schmid_ca, climb_ca, crystal_type, common_material_parameter_ptr )
    

!     if (n_gauss==0) error stop "n_gauss =0 after reading phase input file"
!     if (AmIMPIMaster().and.write_detailed_log_to_screen) write(*,*) "n_gauss ", n_gauss
!     n_ss = sum(n_ss_per_mode)
!     if (AmIMPIMaster().and.write_detailed_log_to_screen) write(*,*) "n_ss ", n_ss

!    call matf_reader%skipEmptyLine()

!   if (use_glide) then
!   !! read material specific parameters for Glide
!   call matf_reader%readParameter("--Glide-model", material_model_name)
!     ! select case(material_model_name%getString())
    ! case ("hutchinson-glide")
    !   allocate(hutchinson_glide_temp)
    !   inelastic_strain_base_ptr => hutchinson_glide_temp
    !   call hutchinson_glide_temp%linkBaseObjects(phase_id, all_mighty_grid_in, sim_all_macro_data, the_bc_object)
    !   call hutchinson_glide_temp%initParameters(phase_id, common_material_parameter_ptr, n_ss_per_mode, schmid_ca, climb_ca, ratio_edge_screw, &
    !                                           use_damage, n_gauss, n_std_dev_gauss_integration, &
    !                                           self_hardening, &
    !                                           latent_hardening, &
    !                                           latent_hardening_other_modes, &
    !                                           phase_material%stiffness_ptr)
    ! case ("wkkct-glide")
    !   allocate(wkkct_glide_temp)
    !   inelastic_strain_base_ptr => wkkct_glide_temp
    !   call wkkct_glide_temp%linkBaseObjects(phase_id, all_mighty_grid_in, sim_all_macro_data, the_bc_object)
    !   call wkkct_glide_temp%initParameters(phase_id, common_material_parameter_ptr,  n_ss_per_mode, schmid_ca, climb_ca, ratio_edge_screw, &
    !                                         use_damage, n_gauss, n_std_dev_gauss_integration, &
    !                                         self_hardening, &
    !                                         latent_hardening, &
    !                                         latent_hardening_other_modes, &
    !                                         burgVectorL, &
    !                                         phase_material%stiffness_ptr)

    ! case default
    !   if (AmIMPIMaster().and.write_detailed_log_to_screen) write(*,*) "unrecognized material model for glide", material_model_name%getString()
    !   if (AmIMPIMaster().and.write_detailed_log_to_screen) write(*,*) "available material models are: "
    !   if (AmIMPIMaster().and.write_detailed_log_to_screen) write(*,*) "wkkct-glide "
    !   if (AmIMPIMaster().and.write_detailed_log_to_screen) write(*,*) "hutchinson-glide "
    !   error stop "abort"
    ! end select

!     call inelastic_strain_base_ptr%readMaterialParametersFromFile(matf_reader)
!     call phase_material%addDeformationMode(inelastic_strain_base_ptr)

!     select case(material_model_name%getString())
!     case ("hutchinson-glide")
!       nullify(hutchinson_glide_temp)
!     case ("wkkct-glide")
!       nullify(wkkct_glide_temp)
!     case default
!       if (AmIMPIMaster().and.write_detailed_log_to_screen) write(*,*) "somthing is wrong when trying nullifying the pointer for ", material_model_name%getString()
!       error stop "abort"
!     end select
!     nullify(inelastic_strain_base_ptr)


!      call matf_reader%skipEmptyLine()
!   endif ! end glide

!   if (use_climb) then
!   ! call matf_reader%readParameter("--Climb-model", material_model_name)
!   !   select case(material_model_name%getString())
!   !   case ("wkkct-climb")
!   !     allocate(wkkct_climb_temp)
!   !     inelastic_strain_base_ptr => wkkct_climb_temp
!   !     call wkkct_climb_temp%linkBaseObjects(phase_id, all_mighty_grid_in, sim_all_macro_data, the_bc_object)
!   !     call wkkct_climb_temp%initParameters(phase_id, common_material_parameter_ptr, n_ss_per_mode, schmid_ca, climb_ca, ratio_edge_screw, &
!   !                                           use_damage, n_gauss, n_std_dev_gauss_integration, &
!   !                                           burgVectorL, phase_material%stiffness_ptr)
!   !   case ("wkkct-climb-irradiation")
!   !     allocate(wkkct_climb_irradiation_temp)
!   !     inelastic_strain_base_ptr => wkkct_climb_irradiation_temp
!   !     call wkkct_climb_irradiation_temp%linkBaseObjects(phase_id, all_mighty_grid_in, sim_all_macro_data, the_bc_object)
!   !     call wkkct_climb_irradiation_temp%initParameters(phase_id, common_material_parameter_ptr, n_ss_per_mode, schmid_ca, climb_ca, ratio_edge_screw, &
!   !                                           use_damage, n_gauss, n_std_dev_gauss_integration, &
!   !                                           burgVectorL,phase_material%stiffness_ptr)
!   !   case ("exponential-climb")
!   !     allocate(exponential_climb_temp)
!   !     inelastic_strain_base_ptr => exponential_climb_temp
!   !     call exponential_climb_temp%linkBaseObjects(phase_id, all_mighty_grid_in, sim_all_macro_data, the_bc_object)
!   !     call exponential_climb_temp%initParameters(phase_id, common_material_parameter_ptr, n_ss_per_mode, schmid_ca, climb_ca, ratio_edge_screw, &
!   !                                           use_damage, n_gauss, n_std_dev_gauss_integration, phase_material%stiffness_ptr)

!   !   case ("transient-climb")
!   !     allocate(transient_climb_temp)
!   !     inelastic_strain_base_ptr => transient_climb_temp
!   !     call transient_climb_temp%linkBaseObjects(phase_id, all_mighty_grid_in, sim_all_macro_data, the_bc_object)
!   !     call transient_climb_temp%initParameters(phase_id, common_material_parameter_ptr, n_ss_per_mode, schmid_ca, climb_ca, ratio_edge_screw, &
!   !                                           use_damage, n_gauss, n_std_dev_gauss_integration, phase_material%stiffness_ptr)
!   !   case default
!   !     if (AmIMPIMaster().and.write_detailed_log_to_screen) write(*,*) "unrecognized material model for climb", material_model_name%getString()
!   !     if (AmIMPIMaster().and.write_detailed_log_to_screen) write(*,*) "available material models are: "
!   !     if (AmIMPIMaster().and.write_detailed_log_to_screen) write(*,*) "transient-climb, exponential-climb, wkkct-climb, and  wkkct-climb-irradiation"
      
!   !     error stop "abort"
!   !   end select

!     call inelastic_strain_base_ptr%readMaterialParametersFromFile(matf_reader)
!     call phase_material%addDeformationMode(inelastic_strain_base_ptr)
!     select case(material_model_name%getString())
!     case ("wkkct-climb")
!       nullify(wkkct_climb_temp)
!     case ("wkkct-climb-irradiation")
!       nullify(wkkct_climb_irradiation_temp)
!     case ("exponential-climb")
!       nullify(exponential_climb_temp)
!     case ("transient-climb")
!       nullify(transient_climb_temp)
!     case default
!       if (AmIMPIMaster().and.write_detailed_log_to_screen) write(*,*) "somthing is worng when trying nullifying the pointer for ", material_model_name%getString()
!       error stop "abort"
!     end select

!     nullify(inelastic_strain_base_ptr)
!     deallocate(schmid_ca, climb_ca); nullify(schmid_ca, climb_ca)
!      call matf_reader%skipEmptyLine()
!   endif ! end climb
!   endif ! end crystal plasticity model

!   if (use_diffusion)then
!     call matf_reader%readParameter("--Diffusion-model", material_model_name)
!     select case(material_model_name%getString())
!     case ("coble-nabarro-diffusion")
!       allocate(nabarro_herring_plus_coble_diffusion_temp)
!       inelastic_strain_base_ptr => nabarro_herring_plus_coble_diffusion_temp
!       call nabarro_herring_plus_coble_diffusion_temp%linkBaseObjects(phase_id, all_mighty_grid_in, sim_all_macro_data, the_bc_object)
!       call nabarro_herring_plus_coble_diffusion_temp%initParameters(phase_id, common_material_parameter_ptr, use_damage)

!     case default
!       if (AmIMPIMaster().and.write_detailed_log_to_screen) write(*,*) "unrecognized material model for climb", material_model_name%getString()
!       if (AmIMPIMaster().and.write_detailed_log_to_screen) write(*,*) "available material models are: "
!       if (AmIMPIMaster().and.write_detailed_log_to_screen) write(*,*) "nabarro_herring_plus_coble_diffusion_temp "
!       error stop "abort"
!     end select

!     call inelastic_strain_base_ptr%readMaterialParametersFromFile(matf_reader)
!     call phase_material%addDeformationMode(inelastic_strain_base_ptr)
!     select case(material_model_name%getString())
!     case ("coble-nabarro-diffusion")
!       nullify(nabarro_herring_plus_coble_diffusion_temp)
!     case default
!       if (AmIMPIMaster().and.write_detailed_log_to_screen) write(*,*) "somthing is worng when trying nullifying the pointer for ", material_model_name%getString()
!       error stop "abort"
!     end select
!      call matf_reader%skipEmptyLine()
!   endif


!   if (use_porosity) then
!     call matf_reader%readParameter("--Porosity-model", material_model_name)
!     select case(material_model_name%getString())
!     case ("BTKTLC")
!       allocate(porosity_base_temp)
!       porosity_base_ptr => porosity_base_temp
!       call porosity_base_temp%linkBaseObjects(phase_id, all_mighty_grid_in, sim_all_macro_data, the_bc_object)
!       call porosity_base_temp%initParameters(phase_id, common_material_parameter_ptr)

!     case default
!       if (AmIMPIMaster().and.write_detailed_log_to_screen) write(*,*) "unrecognized material model for porosity ", material_model_name%getString()
!       if (AmIMPIMaster().and.write_detailed_log_to_screen) write(*,*) "available porosity models are: "
!       if (AmIMPIMaster().and.write_detailed_log_to_screen) write(*,*) "BTKTLC "
!       error stop "abort"
!     end select

!     call porosity_base_ptr%readMaterialParametersFromFile(matf_reader)
!     call phase_material%setPorosity(porosity_base_ptr)
!     select case(material_model_name%getString())
!     case ("BTKTLC")
!       nullify(porosity_base_temp)
!     case default
!       if (AmIMPIMaster().and.write_detailed_log_to_screen) write(*,*) "somthing is worng when trying nullifying the pointer for ", material_model_name%getString()
!       error stop "abort"
!     end select
!   endif
!   endif
!   nullify(common_material_parameter_ptr)
!   call matf_reader%closeTextFile()
! end subroutine


! subroutine readCrystalParameters(matf_reader, &
!                                  ratio_edge_screw, &
!                                  n_ss_per_mode, &
!                                  burgVectorL, &
!                                  self_hardening, &
!                                  latent_hardening, &
!                                  latent_hardening_other_modes, &
!                                  n_gauss, &
!                                  n_std_dev_gauss_integration, &
!                                  ss_direction_ca, ss_normal_ca, crystal_type_int, &
!                                  common_material_parameter_ptr )
!   use common_material_parameter_mod
!   use read_from_file_utils
!   use print_utils_mod, only : printToScreen
!   use tensor_math_mod, only : vector3OuterProduct, doubleContraction,doubleContractionBetweenTwoMatrix
!   use embedded_slip_systems_mod, only : getNormalizedSlipSystemDirectionandNormal
!   use embedded_slip_systems_mod, only : getNormalizedLoopsNormal
!   use log_file_mod, only : write_detailed_log_to_screen
!   use embedded_slip_systems_mod, only : crystal_type_enum, stringToEnumCrystalType, &
!                                         FCC, BCC, HCP
!   use mpi_useful_routines_mod, only : AmIMPIMaster
!   implicit none
!   class(file_Reader), intent(inout) :: matf_reader
!   real(k_real), intent(out) :: ratio_edge_screw
!   integer, dimension(:), pointer, intent(inout) :: n_ss_per_mode
!   real(k_real), pointer, dimension(:), intent(inout) :: burgVectorL, &
!                                                         self_hardening, &
!                                                         latent_hardening
!   real(k_real), intent(out) :: latent_hardening_other_modes

!   integer, intent(out) :: n_gauss
!   real(k_real), intent(out) :: n_std_dev_gauss_integration
!   real(k_real), pointer, dimension(:,:), intent(inout) :: ss_direction_ca, ss_normal_ca
!   real(k_real), allocatable, dimension(:,:) ::  ss_direction_ca_temp, ss_normal_ca_temp

!   ! temporary variables only used by this subroutine
!   type(string_type) :: crystal_type
!   type(string_array) :: slip_modes,dislocation_loops
!   real(k_real), pointer, dimension(:,:) :: temporary_ss_direction_ca=>null(), temporary_ss_normal_ca=>null()
!   real(k_real), pointer, dimension(:,:) :: temporary_loop_normal_ca=>null()
!   type(string_array) :: dummy_string_array
!   type(string_type) :: hardening_matrix_filename
!   integer :: n_slip_modes, i, n_ss_total, n_ss
!   integer :: n_loop, j, ss_idx
!   real(k_real) :: ca_ratio
!   integer(kind(crystal_type_enum)), intent(inout) :: crystal_type_int
!   class(common_material_parameter), pointer, intent(inout) :: common_material_parameter_ptr
!   real(k_real), dimension(3,3) :: ss_normal_ca_tensor,loop_normal_ca_tensor
!   logical :: read_hardening_matrix_from_file

!   read_hardening_matrix_from_file = .False.
!   if (associated(ss_direction_ca)) error stop "readCrystalParameters ss_direction_ca already associated"
!   if (associated(ss_normal_ca)) error stop "readCrystalParameters ss_normal_ca already associated"

!   call matf_reader%readParameter("crystal-type", crystal_type)
!   call stringToEnumCrystalType(crystal_type%getString(),  crystal_type_int)
!   call common_material_parameter_ptr%setCrystalType(crystal_type_int)

!   ! 
!   select case (crystal_type_int)
!   case (FCC)
!     ! nothing to do
!   case (BCC)
!     ! nothing to do
!   case (HCP)
!     call matf_reader%readParameter("c/a-ratio", ca_ratio)
!     call common_material_parameter_ptr%setCAratio(ca_ratio)
!   case default
!     if (AmIMPIMaster().and.write_detailed_log_to_screen) write(*,*) "the selected crystal type ", crystal_type%getString(), "is not implemented"
!     if (AmIMPIMaster().and.write_detailed_log_to_screen) write(*,*) "available options are: FCC!"
!     stop
!   end select

!   call matf_reader%readParameter("ratio-edge-screw", ratio_edge_screw)
!   call matf_reader%readParameter("n-slip-modes", n_slip_modes)
!   !sanity check before other improvements
!   ! if(n_slip_modes>1) error stop "For now only one slip mode is supported. Remove when support for more slip modes is implemented"
!   if(n_slip_modes<1) error stop "you need at least 1 slip modes to use crystal plasticity. Abort"
!   allocate(n_ss_per_mode(n_slip_modes))
!   call matf_reader%readVectorParameter("slip-modes", n_slip_modes, slip_modes)

!   do i =1,n_slip_modes
!     if(.not.(slip_modes%strings(i)%startsWith(crystal_type%getString()))) then
!       if (AmIMPIMaster().and.write_detailed_log_to_screen) write(*,*) " you can't select the slip mode ", slip_modes%getStringByIndex(i), " for a ", crystal_type%getString(), " crystal type"
!       error stop "Abort!"
!     endif
!     call common_material_parameter_ptr%addSSTypeString(slip_modes%getStringByIndex(i))
!   enddo

!   n_ss_total = 0
!   do i = 1,n_slip_modes
!     select case (crystal_type%getString())
!     case ("HCP")
!       call getNormalizedSlipSystemDirectionAndNormal(slip_modes%getStringByIndex(i), n_ss, temporary_ss_direction_ca, temporary_ss_normal_ca, ca_ratio)
!     case default
!       call getNormalizedSlipSystemDirectionAndNormal(slip_modes%getStringByIndex(i), n_ss, temporary_ss_direction_ca, temporary_ss_normal_ca)
!     end select
!     n_ss_per_mode(i) = n_ss
!     if (i==1) then
!       allocate(ss_direction_ca(3,n_ss), ss_normal_ca(3,n_ss))
!     else
!       allocate(ss_direction_ca_temp(3,n_ss_total), ss_normal_ca_temp(3,n_ss_total))
!       ss_direction_ca_temp = ss_direction_ca
!       ss_normal_ca_temp = ss_normal_ca
!       deallocate(ss_direction_ca, ss_normal_ca)
!       allocate(ss_direction_ca(3,n_ss_total+n_ss), ss_normal_ca(3,n_ss_total+n_ss))
!       ss_direction_ca(:,1:n_ss_total) = ss_direction_ca_temp
!       ss_normal_ca(:,1:n_ss_total) = ss_normal_ca_temp
!       deallocate(ss_direction_ca_temp, ss_normal_ca_temp)
!     endif
!     ss_direction_ca(:,n_ss_total+1:n_ss_total+n_ss) = temporary_ss_direction_ca
!     ss_normal_ca(:,n_ss_total+1:n_ss_total+n_ss) = temporary_ss_normal_ca

!     deallocate(temporary_ss_direction_ca, temporary_ss_normal_ca)
!     nullify(temporary_ss_direction_ca, temporary_ss_normal_ca)
!     n_ss_total = n_ss_total + n_ss
!   enddo

!   ! call printToScreen(transpose(ss_direction_ca), "ss_direction_ca")
!   ! call printToScreen(transpose(ss_normal_ca), "ss_normal_ca")
!   call matf_reader%readParameter("read-hardening-matrix-from-file[TRUE/FALSE]", read_hardening_matrix_from_file)
!   if (read_hardening_matrix_from_file) then
!     call matf_reader%readParameter("hardening-matrix-filename", hardening_matrix_filename)
!     call common_material_parameter_ptr%setHardeningMatrixFromFile(hardening_matrix_filename%getString(), n_ss_total)
!   else
!     call matf_reader%readVectorParameter("self-hardening-coeff", n_slip_modes, self_hardening)
!     call matf_reader%readVectorParameter("latent-hardening-coeff", n_slip_modes, latent_hardening)
!     call matf_reader%readParameter("latent-hardening-other-modes", latent_hardening_other_modes)
!   endif 
  
!   call matf_reader%readVectorParameter("burger-vector-length", n_slip_modes, burgVectorL)

!   allocate(common_material_parameter_ptr%number_disloc_loop_type)
!   call matf_reader%readParameter("n-dislocation-loop-type", common_material_parameter_ptr%number_disloc_loop_type)
!   write(*,*) common_material_parameter_ptr%number_disloc_loop_type, "n_dislocation_loop_type"

!   if(common_material_parameter_ptr%number_disloc_loop_type.eq.0) then
!     call matf_reader%readLineAndCheckStringsAreEqual("dislocation-loops", dummy_string_array)
!   else
!     if (common_material_parameter_ptr%number_disloc_loop_type.lt.0) then
!       if (AmIMPIMaster().and.write_detailed_log_to_screen) write(*,*) " the number of dislocation loop types must be >=0, instead I have ", common_material_parameter_ptr%number_disloc_loop_type
!       stop
!     endif
!     call matf_reader%readVectorParameter("dislocation-loops", common_material_parameter_ptr%number_disloc_loop_type, dislocation_loops)
!     allocate(common_material_parameter_ptr%loop_crystallography_factor(n_ss_total,common_material_parameter_ptr%number_disloc_loop_type))
!     common_material_parameter_ptr%loop_crystallography_factor = 0._k_real
!     do i = 1,common_material_parameter_ptr%number_disloc_loop_type
!       call getNormalizedLoopsNormal(dislocation_loops%getStringByIndex(i), n_loop, temporary_loop_normal_ca)
!       do ss_idx = 1, n_ss_total
!         ss_normal_ca_tensor = vector3OuterProduct(ss_normal_ca(:,ss_idx),ss_normal_ca(:,ss_idx))
!         do j=1,n_loop
!           loop_normal_ca_tensor = vector3OuterProduct(temporary_loop_normal_ca(:,j),temporary_loop_normal_ca(:,j))
!           common_material_parameter_ptr%loop_crystallography_factor(ss_idx,i) = common_material_parameter_ptr%loop_crystallography_factor(ss_idx,i) + (1._k_real/n_loop)*doubleContractionBetweenTwoMatrix(ss_normal_ca_tensor,loop_normal_ca_tensor)
!         end do
!        end do
!       deallocate(temporary_loop_normal_ca)
!       nullify(temporary_loop_normal_ca)
!     end do
!   endif

!   call matf_reader%skipEmptyLine()

!   call matf_reader%readLineAndCheckStringsAreEqual("--Integral-Formulation", dummy_string_array)
!   call matf_reader%readParameter("n-gauss", n_gauss)
!   call matf_reader%readParameter("n-standard-deviation-gaussian-integration", n_std_dev_gauss_integration)

!   if (n_gauss<=0) then
!     if (AmIMPIMaster().and.write_detailed_log_to_screen) write(*,*) " the number of gauss points must be >=0, instead I have ", n_gauss
!     stop
!   endif
!   if (n_std_dev_gauss_integration<0) then
!     if (AmIMPIMaster().and.write_detailed_log_to_screen) write(*,*) " the gauss point standard deviation must be >0, instead I have ", n_std_dev_gauss_integration
!     stop
!   endif
! end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine parse_options_file()

  use GLOBAL, only: all_mighty_grid_global, bc_file_name, the_bc_object, microstrucutre_info, &
  sim_all_macro_data, field_file_base_name, dump_file_base_name, restart_file_name, igas, imposed_stressBC_abs_tol, material_stress_abs_tol, &
  max_material_iterations, max_nl_iterations, min_meaningful_strain_rate_to_check_convergence, nl_strain_rate_abs_tol, nl_strain_rate_rel_tol, &
  nl_stress_rel_tol, num_restart_file_to_keep, phase_hardening, restart_from_dump, strain_rate_value_for_using_abs_tol, tdot, tdot_max, &
  tdot_pre, tdot_ref, time, update_texture, write_field_files_every_n_steps, write_non_converged_mp_info_to_file, write_restart_file_every_n_steps, &
  the_bc_object, tdot_min
  use read_from_file_utils
  use termination_criterion_mod, only : readTerminationCriteriaOptions
  use string_module
  use time_march_mod , only : readGlobalTimemarchTolerances
  use log_file_mod, only : write_detailed_log_to_screen
  use all_mighty_grid_mod, only : all_mighty_grid_type
  use all_grid_data_mod, only : all_grid_data
  implicit none

  real(k_real) :: t0, dt, dtmin, dtmax
  integer, parameter, dimension(6) :: component=(/1,6,5, &
                                                    2,4, &
                                                      3/)
  integer :: num_phases, n_points_xyz(3)
  type(string_type) :: read_string_buffer
  type(string_array) :: dummy_string_array
  type(file_reader) :: optf_reader
  
  call optf_reader%openReadTextFile('options.in')

  call optf_reader%readLineAndCheckStringsAreEqual("--Boundary-Conditions", dummy_string_array)
  call optf_reader%readParameter("-boundary-conditions-file-name", read_string_buffer)
  call bc_file_name%setString(read_string_buffer%getString())
  call the_bc_object%readBCFile(bc_file_name%getString())
  call the_bc_object%linkBaseObjects(sim_all_macro_data)
  call read_string_buffer%resetString()
  call optf_reader%skipEmptyLine()

  call optf_reader%readLineAndCheckStringsAreEqual("--Microstructure-info", dummy_string_array)
  call microstrucutre_info%readMicorstructureInfoFromOptionFile(optf_reader)
  num_phases = microstrucutre_info%getNumPhases()
  
  call all_mighty_grid_global%AMGinit(microstrucutre_info%getNumPhases())
  
  call microstrucutre_info%getGlobalDimensions(n_points_xyz)

  call microstrucutre_info%getIsGasPhaseArrayPtr(igas)

  call optf_reader%skipEmptyLine()


  call optf_reader%readLineAndCheckStringsAreEqual("--Simulation-options", dummy_string_array)
  call the_bc_object%getTimeSimulationBegins(t0)
  call optf_reader%readParameter("-initial-dt", dt)
  call optf_reader%readParameter("-dtmin", dtmin)
  call optf_reader%readParameter("-dtmax", dtmax)
  call optf_reader%readParameter("-update-texture", update_texture)
  call optf_reader%readParameter("-phase-hardening", phase_hardening)
  call optf_reader%skipEmptyLine()

  !initialize time object and set pointers to global variables
  call sim_all_macro_data%sim_time%init(t0, dt, dtmin, dtmax)
  call sim_all_macro_data%sim_time%getDeltaTimePointer(tdot)
  call sim_all_macro_data%sim_time%getTimePointer(time)
  call sim_all_macro_data%sim_time%getDeltaTimeOldPointer(tdot_pre)
  call sim_all_macro_data%sim_time%getDeltaTimeReferencePointer(tdot_ref)
  call sim_all_macro_data%sim_time%getDeltaTMinPointer(tdot_min)
  call sim_all_macro_data%sim_time%getDeltaTMaxPointer(tdot_max)

  call optf_reader%readLineAndCheckStringsAreEqual("--Numerical-options", dummy_string_array)
  call optf_reader%readParameter("-max-inner-loop-iterations", max_material_iterations)
  call optf_reader%readParameter("-material-stress-absolute-tolerance", material_stress_abs_tol)
  call optf_reader%readParameter("-max-outer-loop-iterations", max_nl_iterations)
  call optf_reader%readParameter("-imposed-stressBC-absolute-tolerance", imposed_stressBC_abs_tol)
  call optf_reader%readParameter("-non-linear-stress-relative-tolerance", nl_stress_rel_tol)
  call optf_reader%readParameter("-strain-rate-to-switch-from-realtive-to-absolute-tolerance[1/s,absolute-below-this-value]", strain_rate_value_for_using_abs_tol)
  call optf_reader%readParameter("-non-linear-strain-rate-relative-tolerance[1/s]", nl_strain_rate_rel_tol)
  call optf_reader%readParameter("-non-linear-strain-rate-absolute-tolerance[1/s]", nl_strain_rate_abs_tol)
  call optf_reader%skipEmptyLine()

  if(nl_strain_rate_abs_tol.ge.strain_rate_value_for_using_abs_tol) then
    error stop "strain-rate-absolute-tolerance[1/s] > strain-rate-to-switch-from-realtive-to-absolute-tolerance. Abort!"
  endif
  if(nl_strain_rate_abs_tol.lt.min_meaningful_strain_rate_to_check_convergence) then
    error stop "strain-rate-absolute-tolerance[1/s] < machine precision . Abort!"
  endif


  call readGlobalTimemarchTolerances(optf_reader)
  call optf_reader%skipEmptyLine()

  call readTerminationCriteriaOptions(optf_reader)
  call optf_reader%skipEmptyLine()

  call optf_reader%readLineAndCheckStringsAreEqual( "--Output-Options", dummy_string_array)
  call optf_reader%readParameter("-write-fields-every-n-steps", write_field_files_every_n_steps)
  call optf_reader%readParameter("-field-file-base-name", read_string_buffer)
  call field_file_base_name%setString(read_string_buffer%getString())
  call optf_reader%readParameter("-write-detailed-log-to-screen[TRUE/FALSE]", write_detailed_log_to_screen)
  call optf_reader%skipEmptyLine()

  call optf_reader%readLineAndCheckStringsAreEqual("--Restart-Options", dummy_string_array)
  call optf_reader%readParameter("-write-restart-files-every-n-steps", write_restart_file_every_n_steps)
  call optf_reader%readParameter("-num-restart-file-to-keep", num_restart_file_to_keep)
  call optf_reader%readParameter("-restart-file-base-name", read_string_buffer)
  call dump_file_base_name%setString(read_string_buffer%getString())
  call read_string_buffer%resetString()
  call optf_reader%readParameter("-is-this-a-restart[TRUE/FALSE]", restart_from_dump)
  call optf_reader%readParameter("-restart-file-name", read_string_buffer)
  call restart_file_name%setString(read_string_buffer%getString())
  call read_string_buffer%resetString()

  ! some sanity checks
  if(num_restart_file_to_keep.lt.0) error stop "-num-restart-file-to-keep must be a postive number"
  if(write_restart_file_every_n_steps.lt.0) error stop "-write-restart-files-every-n-steps must be a positive number. use 0 to never write a dump file"
  call optf_reader%skipEmptyLine()

  call optf_reader%readLineAndCheckStringsAreEqual("--Debug-Options", dummy_string_array)
  call optf_reader%readParameter("--dump-not-converged-material-point-info-to-file", write_non_converged_mp_info_to_file)

  call optf_reader%closeTextFile()


end subroutine parse_options_file


! subroutine parse_phase_material_file(n_ss, n_gauss )
!   use global, only : microstrucutre_info, phase_material_ptr, &
!   phase_material_array, all_mighty_grid_global, sim_all_macro_data, the_bc_object
!   implicit none

!   class(all_grid_data), pointer :: current_phase_ptr => null()
!   integer :: iph, n_points_xyz(3), num_phases
!   integer, allocatable, dimension(:), intent(out) :: n_ss, n_gauss 

!   num_phases = microstrucutre_info%getNumPhases()

!   do iph=1,num_phases
!     if  (microstrucutre_info%getIsGasPhaseByIndex(iph)) then
!       error stop "gas phase not implemented yet. Abort"
!     else 

!       allocate(phase_material_ptr)
!       call all_mighty_grid_global%AMGGetPhaseGridPtrByIdx(iph, current_phase_ptr)
!       call phase_material_ptr%linkBaseObjects(iph, all_mighty_grid_global, sim_all_macro_data, the_bc_object)
!       call phase_material_ptr%init()
!       call data_crystal_andrea(microstrucutre_info%getPhaseFileNameByIndex(iph), iph, phase_material_ptr, all_mighty_grid_global, n_ss_ph, n_gauss_ph)
!       call phase_material_array%addElement(phase_material_ptr)
!       nullify(phase_material_ptr, current_phase_ptr)
!       n_ss(iph) = n_ss_ph
!       n_gauss(iph) = n_gauss_ph
!     endif
!   enddo
!   call microstrucutre_info%getIsGasPhaseArrayPtr(igas)
! end subroutine

end module fileIO
