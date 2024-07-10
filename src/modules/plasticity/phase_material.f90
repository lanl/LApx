module phase_material_mod
#include "macro_debug.fpp"
use kinds
use diffusion_base_mod
use isotropic_diffusion_mod
use glide_base_mod
use glide_mod
use climb_base_mod
use climb_mod
use inelastic_strain_mod, only : inelastic_strain_base, inelastic_strain_array
use polymorphic_dtype_array_mod, only : dtype_array_ptr
use all_grid_data_mod, only : all_grid_data
use all_mighty_grid_mod, only : all_mighty_grid_type
use simumaltion_macro_data_mod, only : sim_all_macro_data_obj
use stiffness_base_mod
use porosity_base_mod
use bc_objects_mod, only : boundary_condition_array_type
use csv_writer_mod, only : csv_writer
use read_from_file_utils, only : file_reader
use common_material_parameter_mod

implicit none


type :: phase_material_type
  type(inelastic_strain_array) :: all_def_mod
  class(stiffness_base), pointer :: stiffness_ptr => null()
  class(porosity_base), pointer :: porosity_ptr => null()
  integer :: n_def_mode = 0

  integer :: phase_id = -1
  logical :: use_crystal_plasticity = .false.
  logical :: has_glide = .false.
  logical :: has_climb = .false.
  logical :: has_diffusion = .false.
  logical :: has_porosity = .false.
  logical :: is_gas_phase = .false.
  logical :: has_damage = .false.
  logical :: use_isotropic_plasticity = .false.

  ! variables required for output
  logical :: write_header = .true.

  ! simulation parameter
  logical :: macro_object_linked =.false.
  class(sim_all_macro_data_obj), pointer :: sim_all_macro_data => null()
  class(all_mighty_grid_type), pointer :: all_mighty_grid => null()
  class(boundary_condition_array_type), pointer :: the_bc_object => null()
  class(all_grid_data), pointer :: common_grid_data => null(), &
                                   phase_grid_data => null()

  real(k_real), pointer :: dt => null(), &
                           dt_max => null(), &
                           temperature => null(), &
                           time => null()
  
   integer(k_int), pointer, dimension(:,:,:) :: phase_grid => null() ! a pointer to the material phase_id
   real(k_real), pointer, dimension(:,:,:) :: phase_fraction => null() ! a pointer to the material phase_id
   real(k_real), pointer, dimension(:,:,:,:,:) :: stiffness66_ph => null(), &
                                                  R_crystal2sample_ph => null(),&
                                                  elastic_strain_ph => null(), &
                                                  elastic_strain_ph_old => null(), &
                                                  elastic_strain_rate_ph => null(), &
                                                  plastic_strain_ph => null(), &
                                                  plastic_strain_ph_old => null(), &
                                                  plastic_strain_rate_ph => null(), &
                                                  total_strain_ph => null(), &
                                                  total_strain_ph_old => null(), &
                                                  total_strain_rate_ph => null(), &
                                                  stress_ph => null(), &
                                                  stress_ph_old => null(), &
                                                  stress_rate_ph => null()
                                                  

   type(csv_writer) :: phase_material_stat_writer, strain_rate_contribution_writer
   type(file_reader) :: matf_reader
   class(common_material_parameter), pointer :: common_material_parameter_ptr => null()

   __DECL_CLASS_UNUSED_THIS__

contains
  procedure :: init => initPhaseMaterial
  procedure :: linkBaseObjects => linkBaseObjectsPhaseMat
  procedure :: addDeformationMode => addDeformationModePhaseMaterial
  procedure :: getElasticStrainRate => getElasticStrainRatePhaseMaterial
  procedure :: getStrainRateaAndStressJacobian => getStrainRateaAndStressJacobianPhaseMaterial
  procedure :: getConstitutiveStrainRate => getConstitutiveStrainRatePhaseMaterial
  procedure :: initGridPointers => initGridPointersPhaseMaterial
  procedure :: addFieldVariablesToGrid => addFieldVariablesToGridPhaseMaterial
  procedure :: initStateVariables => initStateVariablesPhaseMaterial
  procedure :: updateStateVariablesInnerLoop => updateStateVariablesInnerLoopPhaseMat
  procedure :: updateStateVariablesOuterLoop => updateStateVariablesOuterLoopPhaseMat
  procedure :: updateStateVariablesStaggered => updateStateVariablesStaggeredPhaseMat
  procedure :: isInitialized => isInitializedPhaseMat
  procedure :: writeInelastiStrainRateToFile => writeInelastiStrainRateToFilePhaseMat
  procedure :: getInealsticRotationRateContributionAtMaterialPoint => getRotationRateContributionPhaseMat
  procedure :: initTexture => initTexturePhaseMaterial
  procedure :: setStiffness => setStiffnessPhaseMaterial
  procedure :: setPorosity => setPorosityPhaseMaterial
  procedure :: checkPorosityUpdate
  procedure :: acceptRejectSolution => acceptRejectSolutionPhaseMaterial
  procedure :: readFromFile => readFromFilePhaseMaterial
  procedure :: writeAverageQuantitiesCSV => writeAverageQuantitiesCSVPhaseMaterial
  procedure :: updateRotationMatrix => updateRotationMatrixPhaseMaterial
  procedure :: updateStiffness => updateStiffnessPhaseMaterial
  procedure :: updatePhaseParameters => updatePhaseParametersPhasePhaseMaterial
end type

contains

  subroutine writeAverageQuantitiesCSVPhaseMaterial(this, write_headers)
    implicit none
    class(phase_material_type), intent(inout) :: this
    class(inelastic_strain_base), pointer :: def_mode =>null()
    integer :: i
    logical, intent(in) :: write_headers

    call this%phase_material_stat_writer%openToAppendCSV()
    call this%phase_material_stat_writer%AppendScalar(this%time, "time", write_headers)
    call this%phase_material_stat_writer%AppendScalar(this%temperature,"temperature", write_headers)

    do i=1,this%n_def_mode
      nullify(def_mode)
      call this%all_def_mod%getElementPtr(i, def_mode)
      call def_mode%writeAverageQuantitiesCSV(this%phase_material_stat_writer, write_headers)
    enddo

    if (this%has_porosity) then
      call this%porosity_ptr%writeAverageQuantitiesCSV(this%phase_material_stat_writer, write_headers)
    endif
    call this%phase_material_stat_writer%closeCSVFile()

  end subroutine

  subroutine initPhaseMaterial(this)
    use number_to_string_mod, only : int2string
    implicit none
    class(phase_material_type), intent(inout) :: this
    if (this%phase_id<0 ) error stop "phase_id not initialized"

    call this%phase_material_stat_writer%createNewCSV("material_properties_stat_ph"//trim(adjustl(int2string(this%phase_id)))//".csv")
    call this%strain_rate_contribution_writer%createNewCSV("contribution_phase"//trim(adjustl(int2string(this%phase_id)))//".csv")
    
  end subroutine

  subroutine linkBaseObjectsPhaseMat(this, phase_id, all_mighty_grid_in, sim_all_macro_data, the_bc_object)
    implicit none
    class(phase_material_type), intent(inout) :: this
    integer, intent(in) :: phase_id
    class(all_mighty_grid_type), intent(in), target :: all_mighty_grid_in
    class(sim_all_macro_data_obj), intent(in), target :: sim_all_macro_data
    class(boundary_condition_array_type), intent(in), target :: the_bc_object

    
    if (this%macro_object_linked) error stop "linkBaseObjects already intialized. abort! "
    


    this%all_mighty_grid => all_mighty_grid_in
    this%phase_id = phase_id
    call this%all_mighty_grid%AMGGetCommonGridPtr(this%common_grid_data) 
    if (this%phase_id<0 ) error stop "phase_id not initialized"
    call this%all_mighty_grid%AMGGetPhaseGridPtrByIdx(this%phase_id, this%phase_grid_data)  

    this%sim_all_macro_data => sim_all_macro_data
    this%the_bc_object => the_bc_object
    call the_bc_object%getTemperaturePointer(this%temperature)
    call sim_all_macro_data%sim_time%getTimePointer(this%time)
    call sim_all_macro_data%sim_time%getDeltaTimePointer(this%dt)
    call sim_all_macro_data%sim_time%getDeltaTMaxPointer(this%dt_max)
    this%macro_object_linked =.true.

  end subroutine

  subroutine updateRotationMatrixPhaseMaterial(this)
    use bunge_mod, only : bungeAnlgesFromRotationMatrixCrystal2Sample, bungeRotationMatrixCrystal2Sample, RodriguesRotationFromDeltaRotationMatrix
    use tensor_math_mod, only : getSkewPart
    implicit none
    class(phase_material_type), intent(inout) :: this
    real(k_real), pointer, dimension(:,:,:,:,:) :: R_crystal2sample=>null(), &
                                                   vel_grad=>null()
                                                   
    integer :: x_s, x_e, y_s, y_e, z_s, z_e, &
               ix, iy, iz
    real(k_real), dimension(3,3) :: DeltaR, Rdot_vel_grad, Rdot_inelastic, RodriguesRotation
    real(k_real) :: ph, th, tm
    nullify(R_crystal2sample, vel_grad)

    if (.not.(this%is_gas_phase)) then
      ! get pointer
    call this%phase_grid_data%getTensor2DataPointerByName("R_crystal2sample", R_crystal2sample)
    call this%common_grid_data%getTensor2DataPointerByName("velocity_gradient", vel_grad)
    call this%all_mighty_grid%AMGGetLoopLimitsRank(x_s, x_e, y_s, y_e, z_s, z_e)

    
    do  iz = z_s,z_e
      do  iy = y_s,y_e
         do  ix = x_s, x_e
             ! get the rotation contribution from the velocity gradient
             Rdot_vel_grad = getSkewPart(vel_grad(:,:,ix,iy,iz))
 
             ! get the rotation contribution from inelastic strains
             call this%getInealsticRotationRateContributionAtMaterialPoint(Rdot_inelastic, ix,iy,iz)
 
             ! compute the total Delta rotation
             DeltaR = ( Rdot_vel_grad - Rdot_inelastic)* this%dt
 
             !delta R to compute a rotation matrix update using rodrigues formula
             call RodriguesRotationFromDeltaRotationMatrix(DeltaR, RodriguesRotation)
 
             ! update the actual rotation matrix
             R_crystal2sample(:,:,ix,iy,iz) = matmul(RodriguesRotation,R_crystal2sample(:,:,ix,iy,iz))
 
             call bungeAnlgesFromRotationMatrixCrystal2Sample(R_crystal2sample(:,:,ix,iy,iz), ph, th, tm )
             call bungeRotationMatrixCrystal2Sample(ph, th, tm, R_crystal2sample(:,:,ix,iy,iz))
          
         end do
      end do
   end do
   
   nullify(R_crystal2sample, vel_grad)
  endif

  end subroutine

  subroutine updateStiffnessPhaseMaterial(this)
    use tensor_math_mod, only : getSymmetricPart
    use change_tensor_basis, only : chg_basis_tensor4_to_matrix66
    use tensor_math_mod, only : rotateTensor4
    implicit none
    class(phase_material_type), intent(inout) :: this
    real(k_real), pointer :: R_crystal2sample(:,:,:,:,:)=>null(), &
                             stiffness_ptr(:,:,:,:,:)=>null(), &
                             effective_porosity_grid(:,:,:)=> null()
    integer :: x_s, x_e, y_s, y_e, z_s, z_e, ix, iy, iz
    real(k_real) :: rotated_stiffness(3,3,3,3), stiffness66_temp(6,6)
    nullify(R_crystal2sample, stiffness_ptr, effective_porosity_grid )

    if (.not.(this%is_gas_phase)) then
      ! get pointer
    call this%phase_grid_data%getTensor2DataPointerByName("R_crystal2sample", R_crystal2sample)
    call this%phase_grid_data%getMatrix66DataPointerByName("stiffness", stiffness_ptr)

    if (this%has_porosity) then
      call this%phase_grid_data%getScalarDataPointerByName("effective_porosity", effective_porosity_grid)
    endif

    call this%all_mighty_grid%AMGGetLoopLimitsRank(x_s, x_e, y_s, y_e, z_s, z_e)
    
    do  iz = z_s,z_e
      do  iy = y_s,y_e
         do  ix = x_s, x_e
          associate(R=>R_crystal2sample(:,:,ix,iy,iz), &
                    stiffness=>stiffness_ptr(:,:,ix,iy,iz) )
              call rotateTensor4(this%stiffness_ptr%stiffness_crystal_axes, R, rotated_stiffness)
              call chg_basis_tensor4_to_matrix66(rotated_stiffness, stiffness66_temp)
              stiffness = getSymmetricPart(stiffness66_temp)

              if (this%has_porosity) then
                stiffness = stiffness*(1._k_real-effective_porosity_grid(ix,iy,iz))
              endif

          end associate
         end do
      end do
    end do
   
    nullify(R_crystal2sample, stiffness_ptr )
    
    if (this%has_porosity) then
       nullify(effective_porosity_grid)
    endif
   endif

  end subroutine

  subroutine setStiffnessPhaseMaterial(this, stiffness)
    implicit none
    class(phase_material_type), intent(inout) :: this
    class(stiffness_base), pointer, intent(inout) :: stiffness

    if (.not.(associated(stiffness))) error stop "setStiffnessPhaseMaterial stiffness is not associated!"
    if (associated(this%stiffness_ptr)) error stop "setStiffnessPhaseMaterial this%stiffness already associated!"

    this%stiffness_ptr => stiffness
    nullify(stiffness)

  end subroutine

  subroutine setPorosityPhaseMaterial(this, porosity_model)
    implicit none
    class(phase_material_type), intent(inout) :: this
    class(porosity_base), pointer, intent(inout) :: porosity_model

    if (.not.(associated(porosity_model))) error stop "setPorosityPhaseMaterial stiffness is not associated!"
    if (associated(this%porosity_ptr)) error stop "setPorosityPhaseMaterial this%porosity_ptr already associated!"

    this%porosity_ptr => porosity_model
    nullify(porosity_model)
    this%has_porosity = .true.

  end subroutine

  subroutine addDeformationModePhaseMaterial(this, new_def_mode)
    class(phase_material_type), intent(inout) :: this
    class(inelastic_strain_base), pointer, intent(inout) :: new_def_mode


    if (.not.(associated(new_def_mode))) error stop "addDeformationModePhaseMaterial new_def_mode is not associated!"

    if (new_def_mode%phase_id.ne.this%phase_id) error stop "addDeformationModePhaseMaterial the provided phase_id is != to the stored phase_id"

    select type(new_def_mode)
    class is(cp_glide_base)
      if (this%has_glide) error stop "addDeformationModePhaseMaterial already has a glide mode"
      this%has_glide = .true.
      this%use_crystal_plasticity = .true.

    class is(cp_climb_base)
      if (this%has_climb) error stop "addDeformationModePhaseMaterial already has a climb mode"
      this%has_climb = .true.
      this%use_crystal_plasticity = .true.

    class is(diffusion_base)
      if (this%has_diffusion) error stop "addDeformationModePhaseMaterial already has a diffusion mode"
      this%has_diffusion = .true.

    class default
      error stop "addDeformationModePhaseMaterial unercognized strain type"
    end select

    call this%all_def_mod%addElementFirst(new_def_mode)
    nullify(new_def_mode)
    this%n_def_mode = this%n_def_mode + 1
  end subroutine

  subroutine initGridPointersPhaseMaterial(this)
    implicit none
    class(phase_material_type), intent(inout) :: this
    class(inelastic_strain_base), pointer :: def_mode =>null()
    integer :: i
    real(k_real), pointer, dimension(:,:,:,:) :: temp_vector_pointer => null()

    associate (phase_grid => this%phase_grid_data)
      call phase_grid%getTensor2DataPointerByName("elastic_strain", this%elastic_strain_ph)
      call phase_grid%getTensor2DataPointerByNameOld("elastic_strain", this%elastic_strain_ph_old)
      call phase_grid%getTensor2DataPointerByName("elastic_strain_rate", this%elastic_strain_rate_ph)
      call phase_grid%getTensor2DataPointerByName("plastic_strain", this%plastic_strain_ph)
      call phase_grid%getTensor2DataPointerByNameOld("plastic_strain", this%plastic_strain_ph_old)
      call phase_grid%getTensor2DataPointerByName("plastic_strain_rate", this%plastic_strain_rate_ph)
      call phase_grid%getTensor2DataPointerByName("total_strain", this%total_strain_ph)
      call phase_grid%getTensor2DataPointerByNameOld("total_strain", this%total_strain_ph_old)
      call phase_grid%getTensor2DataPointerByName("total_strain_rate", this%total_strain_rate_ph)
      call phase_grid%getTensor2DataPointerByName("stress", this%stress_ph)
      call phase_grid%getTensor2DataPointerByNameOld("stress", this%stress_ph_old)
      call phase_grid%getTensor2DataPointerByName("stress_rate", this%stress_rate_ph)
      call phase_grid%getMatrix66DataPointerByName("stiffness", this%stiffness66_ph)
      call phase_grid%getTensor2DataPointerByName("R_crystal2sample", this%R_crystal2sample_ph)
    end associate
    call this%common_grid_data%getGenericVectorDataPointerByName("phase_fraction", temp_vector_pointer)
    this%phase_fraction => temp_vector_pointer(this%phase_id,:,:,:)
    nullify(temp_vector_pointer)
    do i=1,this%n_def_mode
      nullify(def_mode)
      call this%all_def_mod%getElementPtr(i, def_mode)
      call def_mode%initGridPointers()
    enddo

    if (this%has_porosity) then
      call this%porosity_ptr%initGridPointers()
    endif

  end subroutine

  subroutine initTexturePhaseMaterial(this)
    class(phase_material_type), intent(inout) :: this
    class(inelastic_strain_base), pointer :: def_mode =>null()
    integer :: i

    do i=1,this%n_def_mode
      nullify(def_mode)
      call this%all_def_mod%getElementPtr(i, def_mode)

      select type(def_mode)
      class is(cp_base)
        call def_mode%initTexture()
      end select
    enddo
  end subroutine

  subroutine addFieldVariablesToGridPhaseMaterial(this)
    use grid_data_var_type_mod
    implicit none
    class(phase_material_type), intent(inout) :: this
    class(inelastic_strain_base), pointer :: def_mode =>null()
    integer :: i

    associate (phase_grid => this%phase_grid_data)
      call phase_grid%addVar("elastic_strain", tensor2)
      call phase_grid%addVar("elastic_strain_rate", tensor2)
      call phase_grid%addVar("plastic_strain", tensor2)
      call phase_grid%addVar("plastic_strain_rate", tensor2)
      call phase_grid%addVar("total_strain", tensor2)
      call phase_grid%addVar("total_strain_rate", tensor2)
      call phase_grid%addVar("stress", tensor2)
      call phase_grid%addVar("stress_rate", tensor2)
      call phase_grid%addVar("stiffness", matrix66)
      call phase_grid%addVar("R_crystal2sample", tensor2)
    end associate

    do i=1,this%n_def_mode
      nullify(def_mode)
      call this%all_def_mod%getElementPtr(i, def_mode)
      call def_mode%addFieldVariablesToGrid()
    enddo

    if (this%has_porosity) then
      call this%porosity_ptr%addFieldVariablesToGrid()
    endif

  end subroutine

  subroutine initStateVariablesPhaseMaterial(this)
    use mpi_variables_mod, only : i_am_mpi_master
    implicit none
    class(phase_material_type), intent(inout) :: this
    class(inelastic_strain_base), pointer :: def_mode =>null()
    integer :: i

#ifdef __DEBUG__
    if(i_am_mpi_master) write(*,*)  "**************************************************************"
    if(i_am_mpi_master) write(*,*)  "*****PHASE MATERIAL STATE VARIABLES INITIALIZATION BEGIN******"
    if(i_am_mpi_master) write(*,*)  ""
#endif

    do i=1,this%n_def_mode

#ifdef __DEBUG__
    if(i_am_mpi_master) write(*,*)  "initializing variables for def mode number ", i, " of ", this%n_def_mode, " in phase ", this%phase_id
#endif

      nullify(def_mode)
      call this%all_def_mod%getElementPtr(i, def_mode)
      call def_mode%initStateVariables()
    enddo

    if (this%has_porosity) then
#ifdef __DEBUG__
  if(i_am_mpi_master) write(*,*)  "initializing porosity state varaibles in phase ", this%phase_id
#endif
      call this%porosity_ptr%initStateVariables()
    endif

#ifdef __DEBUG__
    if(i_am_mpi_master) write(*,*)  ""
    if(i_am_mpi_master) write(*,*)  "*****PHASE MATERIAL STATE VARIABLES INITIALIZATION BEGIN******"
    if(i_am_mpi_master) write(*,*)  "**************************************************************"
#endif

  end subroutine

  subroutine updateStateVariablesInnerLoopPhaseMat(this, ix_in, iy_in, iz_in)
    implicit none
    class(phase_material_type), intent(inout) :: this
    integer, intent(in) :: ix_in, iy_in, iz_in
    class(inelastic_strain_base), pointer :: def_mode =>null()
    integer :: i

    do i=1,this%n_def_mode
      nullify(def_mode)
      call this%all_def_mod%getElementPtr(i, def_mode)
      call def_mode%updateStateVariablesInnerLoop(ix_in, iy_in, iz_in)
    enddo

  end subroutine

  subroutine updateStateVariablesOuterLoopPhaseMat(this)
#ifdef __DEBUG__
    use mpi_variables_mod, only : i_am_mpi_master
#endif
    implicit none
    class(phase_material_type), intent(inout) :: this
    class(inelastic_strain_base), pointer :: def_mode =>null()
    integer :: i

    ! updating varaibles that are going top be present for each phase
    this%elastic_strain_ph = this%elastic_strain_rate_ph * this%dt + this%elastic_strain_ph_old
    this%plastic_strain_ph = this%plastic_strain_rate_ph * this%dt + this%plastic_strain_ph_old
    this%total_strain_ph = this%total_strain_rate_ph * this%dt + this%total_strain_ph_old
    this%stress_ph = this%stress_rate_ph * this%dt + this%stress_ph_old
    do i=1,this%n_def_mode
      nullify(def_mode)
      call this%all_def_mod%getElementPtr(i, def_mode)
      call def_mode%updateStateVariablesOuterLoop()
    enddo

    if (this%has_porosity) then
#ifdef __DEBUG__
  if(i_am_mpi_master) write(*,*)  "updating porosity in phase ", this%phase_id
#endif
      call this%porosity_ptr%updateStateVariablesOuterLoop()
      call this%updateStiffness()
    endif

  end subroutine

  subroutine updateStateVariablesStaggeredPhaseMat(this)
    use mpi_variables_mod, only : i_am_mpi_master
    implicit none
    class(phase_material_type), intent(inout) :: this
    class(inelastic_strain_base), pointer :: def_mode =>null()
    integer :: i

#ifdef __DEBUG__
    if(i_am_mpi_master) write(*,*)  "**************************************************************"
    if(i_am_mpi_master) write(*,*)  "***PHASE MATERIAL STATE VARIABLES STAGGER UPDATE BEGIN********"
    if(i_am_mpi_master) write(*,*)  ""
#endif

    do i=1,this%n_def_mode

#ifdef __DEBUG__
    if(i_am_mpi_master) write(*,*)  "updating variables for def mode number ", i, " of ", this%n_def_mode, " in phase ", this%phase_id
#endif

      nullify(def_mode)
      call this%all_def_mod%getElementPtr(i, def_mode)
      call def_mode%updateStateVariablesStaggered()
    enddo


    if (this%has_porosity) then
#ifdef __DEBUG__
  if(i_am_mpi_master) write(*,*)  "updating porosity in phase ", this%phase_id
#endif
      call this%porosity_ptr%updateStateVariablesStaggered()
      call this%updateStiffness()
    endif

#ifdef __DEBUG__
    if(i_am_mpi_master) write(*,*)  ""
    if(i_am_mpi_master) write(*,*)  "***PHASE MATERIAL STATE VARIABLES STAGGER UPDATE END**********"
    if(i_am_mpi_master) write(*,*)  "**************************************************************"
#endif

  end subroutine


  subroutine updatePhaseParametersPhasePhaseMaterial(this)
    implicit none
    class(phase_material_type), intent(inout) :: this
    class(inelastic_strain_base), pointer :: def_mode =>null()
    integer :: i

    do i=1,this%n_def_mode
      nullify(def_mode)
      call this%all_def_mod%getElementPtr(i, def_mode)
      call def_mode%updatePhaseParameters()
    enddo


    if (this%has_porosity) then
      call this%porosity_ptr%updatePhaseParameters()
      call this%updateStiffness()
    endif

  end subroutine

  subroutine checkPorosityUpdate(this, need_to_update)
    implicit none
    class(phase_material_type), intent(inout) :: this
    logical, intent(out) :: need_to_update

    need_to_update =.false.
    if (this%has_porosity) then
#ifdef __DEBUG__
  write(*,*)  "checking porosity update", this%phase_id
#endif
      call this%porosity_ptr%updateStateVariablesStaggered()
      if (this%porosity_ptr%total_porosity.gt.0._k_real) need_to_update=.true.
      call this%porosity_ptr%initStateVariables()
      call this%updateStiffness()
    endif

  end subroutine

  subroutine isInitializedPhaseMat(this)
    implicit none
    class(phase_material_type), intent(inout) :: this
    class(inelastic_strain_base), pointer :: def_mode =>null()
    integer :: i

    do i=1,this%n_def_mode
      nullify(def_mode)
      call this%all_def_mod%getElementPtr(i, def_mode)
      call def_mode%isInitialized()
    enddo

  end subroutine

  subroutine writeInelastiStrainRateToFilePhaseMat(this)
    use write_to_file_utils_mod
    use number_to_string_mod
    use mpi_useful_routines_mod, only : MPIAverageVoxelWeightGridMatrix, MPIAverageGridMatrix, MPIAverageGridScalar
    use change_tensor_basis, only : chg_basis_vector6_to_tensor2
    implicit none
    class(phase_material_type), intent(inout) :: this
    class(inelastic_strain_base), pointer :: def_mode =>null()
    integer :: i
    real(k_real) ::avg_tensor(3,3), avg_scalar
    integer :: num_write, idx

    num_write = 1
    if (this%write_header) num_write = num_write +1
  
    do idx=1,num_write
      call this%strain_rate_contribution_writer%openToAppendCSV()
      call this%strain_rate_contribution_writer%AppendScalar(this%time, "time", this%write_header)
      call this%strain_rate_contribution_writer%AppendScalar(this%temperature, "temperature", this%write_header)
      call MPIAverageGridScalar(this%phase_fraction, avg_scalar)
      call this%strain_rate_contribution_writer%AppendScalar(avg_scalar, "phase_fraction", this%write_header)

      call MPIAverageGridMatrix(this%stress_ph, avg_tensor)
      call this%strain_rate_contribution_writer%AppendTensor(avg_tensor, "stress","stress_ph", this%write_header)
      call MPIAverageGridMatrix(this%stress_rate_ph, avg_tensor)
      call this%strain_rate_contribution_writer%AppendTensor(avg_tensor, "stress","stress_rate_ph",this%write_header)

      call MPIAverageVoxelWeightGridMatrix(this%total_strain_ph, this%phase_fraction, avg_tensor)
      call this%strain_rate_contribution_writer%AppendTensor(avg_tensor, "strain", "total_strain_ph",this%write_header)
      call MPIAverageVoxelWeightGridMatrix(this%total_strain_rate_ph, this%phase_fraction, avg_tensor)
      call this%strain_rate_contribution_writer%AppendTensor(avg_tensor, "strain", "total_strain_rate_ph",this%write_header)

      call MPIAverageVoxelWeightGridMatrix(this%elastic_strain_ph, this%phase_fraction, avg_tensor)
      call this%strain_rate_contribution_writer%AppendTensor(avg_tensor, "strain", "elastic_strain_ph",this%write_header)
      call MPIAverageVoxelWeightGridMatrix(this%elastic_strain_rate_ph, this%phase_fraction, avg_tensor)
      call this%strain_rate_contribution_writer%AppendTensor(avg_tensor, "strain", "elastic_strain_rate_ph",this%write_header)

      if (this%n_def_mode.ge.1) then
        call MPIAverageVoxelWeightGridMatrix(this%plastic_strain_ph, this%phase_fraction, avg_tensor)
        call this%strain_rate_contribution_writer%AppendTensor(avg_tensor, "strain", "plastic_strain_ph",this%write_header)
        call MPIAverageVoxelWeightGridMatrix(this%plastic_strain_rate_ph, this%phase_fraction, avg_tensor)
        call this%strain_rate_contribution_writer%AppendTensor(avg_tensor, "strain", "plastic_strain_rate_ph",this%write_header)
        do i=1,this%n_def_mode
          call this%all_def_mod%getElementPtr(i, def_mode)
          call def_mode%writeInelastiStrainRateToFile(this%strain_rate_contribution_writer, write_headers=this%write_header)
        enddo
      endif

      call this%strain_rate_contribution_writer%closeCSVFile()
      this%write_header = .false.
    enddo

  end subroutine

  subroutine getElasticStrainRatePhaseMaterial(this, stress6, edot_el, dedot_el_dsigma,  ix, iy, iz)
    use change_tensor_basis, only : chg_basis_tensor2_to_vector6_fun, chg_basis_vector6_to_tensor2
    implicit none
    class(phase_material_type), intent(inout) :: this
    real(k_real), target, intent(in) ::stress6(6)
    real(k_real), dimension(6), intent(out) :: edot_el
    real(k_real), dimension(6,6), intent(out) :: dedot_el_dsigma
    integer, intent(in) :: ix, iy, iz
    real(k_real), dimension(6,6) ::S

    edot_el =0_k_real
    dedot_el_dsigma = 0_k_real
    
    call lapackInverseSymmetric(this%stiffness66_ph(:,:,ix,iy,iz) ,S)

    edot_el = matmul(S, (stress6-chg_basis_tensor2_to_vector6_fun(this%stress_ph_old(:,:, ix,iy,iz)))/this%dt)
    dedot_el_dsigma = S/this%dt
  end subroutine

  subroutine getStrainRateaAndStressJacobianPhaseMaterial(this, stress6, epsilon_dot_total, depsilon_dot_dstress_total, ix, iy, iz)
    use change_tensor_basis, only : chg_basis_vector6_to_tensor2
    implicit none
    class(phase_material_type), intent(inout) :: this
    real(k_real), target, intent(in) ::stress6(6)
    real(k_real), dimension(6), intent(out) :: epsilon_dot_total
    real(k_real), dimension(6,6), intent(out) :: depsilon_dot_dstress_total
    integer, intent(in) :: ix, iy, iz

    real(k_real), dimension(6) :: epsilon_dot_mode
    real(k_real), dimension(6,6) :: depsilon_dot_dstress_mode
    integer :: i
    class(inelastic_strain_base), pointer :: def_mode =>null()

    epsilon_dot_total =0_k_real
    depsilon_dot_dstress_total = 0_k_real
    do i=1,this%n_def_mode
      call this%all_def_mod%getElementPtr(i, def_mode)
      call def_mode%getStrainRateaAndStressJacobian( stress6, epsilon_dot_mode, depsilon_dot_dstress_mode, ix, iy, iz )
      epsilon_dot_total = epsilon_dot_total + epsilon_dot_mode
      depsilon_dot_dstress_total = depsilon_dot_dstress_total + depsilon_dot_dstress_mode
    enddo
  end subroutine

  subroutine getConstitutiveStrainRatePhaseMaterial(this, stress6, edot_total, dedot_total_dstress, edot_el, edot_inel, ix, iy, iz)
    use change_tensor_basis, only : chg_basis_vector6_to_tensor2, chg_basis_tensor2_to_vector6_fun
    implicit none
    class(phase_material_type), intent(inout) :: this
    real(k_real), target, intent(in) ::stress6(6)
    real(k_real), dimension(6), intent(out) :: edot_total, edot_el, edot_inel
    real(k_real), dimension(6,6), intent(out) :: dedot_total_dstress
    integer, intent(in) :: ix, iy, iz
    real(k_real), dimension(6,6) :: dedot_el_dstress, dedot_inel_dstress 
    real(k_real), dimension(6) :: stress_rate

    call this%getElasticStrainRate(stress6, edot_el, dedot_el_dstress,  ix, iy, iz)
    call this%getStrainRateaAndStressJacobian(stress6, edot_inel, dedot_inel_dstress, ix, iy, iz)

    stress_rate = (stress6-chg_basis_tensor2_to_vector6_fun( this%stress_ph_old(:,:,ix,iy,iz) ))/this%dt
    edot_total = edot_el + edot_inel
    dedot_total_dstress = dedot_el_dstress+dedot_inel_dstress

    call chg_basis_vector6_to_tensor2(stress_rate, this%stress_rate_ph(:,:,ix,iy,iz))
    call chg_basis_vector6_to_tensor2(edot_el, this%elastic_strain_rate_ph(:,:,ix,iy,iz) )
    call chg_basis_vector6_to_tensor2(edot_inel, this%plastic_strain_rate_ph(:,:,ix,iy,iz) )
    call chg_basis_vector6_to_tensor2(edot_total, this%total_strain_rate_ph(:,:,ix,iy,iz) )

  end subroutine

  subroutine acceptRejectSolutionPhaseMaterial(this, dt_max, accept_solution_flag)
    use mpi_variables_mod, only : i_am_mpi_master
    use log_file_mod, only : writeToScreen, write_detailed_log_to_screen
    implicit none
    class(phase_material_type), intent(inout) :: this
    real(k_real), intent(out) :: dt_max
    logical, intent(out) :: accept_solution_flag
    real(k_real) :: dt_max_mode
    logical :: accept_solution_flag_mode
    class(inelastic_strain_base), pointer :: def_mode =>null()
    integer :: i

    call writeToScreen("***************************************************************")
    call writeToScreen("***PHASE MATERIAL STATE VARIABLES ACCEPT/REJECT SOLUTION*******")
    call writeToScreen("")

    dt_max = this%dt_max
    accept_solution_flag = .true.
    ! dt_max = this%dtmax
    do i=1,this%n_def_mode

    if(i_am_mpi_master.and.write_detailed_log_to_screen) write(*,*)  "checking solution quality for def mode number ", i, " of ", this%n_def_mode, " in phase ", this%phase_id

      nullify(def_mode)
      call this%all_def_mod%getElementPtr(i, def_mode)
      call def_mode%acceptRejectSolution(dt_max_mode, accept_solution_flag_mode)
      accept_solution_flag = accept_solution_flag.and.accept_solution_flag_mode
      dt_max = min(dt_max_mode, dt_max)

    enddo

    if (this%has_porosity) then
      if(i_am_mpi_master.and.write_detailed_log_to_screen) write(*,*)  "checking solution quality for porosity in phase ", this%phase_id

      call this%porosity_ptr%acceptRejectSolution(dt_max_mode, accept_solution_flag_mode)
      accept_solution_flag = accept_solution_flag.and.accept_solution_flag_mode
      dt_max = min(dt_max_mode, dt_max)

    endif

    call writeToScreen("")
    call writeToScreen("***PHASE MATERIAL STATE VARIABLES ACCEPT/REJECT SOLUTION END***")
    call writeToScreen("***************************************************************")
    call writeToScreen("")

  end subroutine

  subroutine getRotationRateContributionPhaseMat(this, R_dot, ix, iy, iz)
    implicit none
    class(phase_material_type), intent(inout) :: this
    real(k_real), target, intent(out) :: R_dot(3,3)
    real(k_real) :: R_dot_mode(3,3)
    integer, intent(in) :: ix, iy, iz
    integer :: i
    class(inelastic_strain_base), pointer :: def_mode =>null()

    R_dot = 0
    do i=1,this%n_def_mode
      call this%all_def_mod%getElementPtr(i, def_mode)
      call def_mode%getRotationRateContribution(R_dot_mode, ix, iy, iz )
      R_dot = R_dot + R_dot_mode
    enddo

  end subroutine

  subroutine readFromFilePhaseMaterial(this, fname, n_ss, n_gauss)
    use string_module, only : string_array
    use stiffness_base_mod, only : readElasticityFromFile
    use inelastic_strain_mod, only : inelastic_strain_base
    use cp_base_mod, only : crystal_paremeters_type
    use glide_mod, only : readMaterialParametersFromFileCPGlide
    use climb_mod, only : readMaterialParametersFromFileCPClimb
    use isotropic_diffusion_mod, only : readMaterialParametersFromFileISOdiffusion
    use porosity_base_mod, only : porosity_base, readMaterialParametersFromFilePorosity
    implicit none
    class(phase_material_type), intent(inout) :: this
    integer, intent(inout) :: n_ss, n_gauss
    character(len=*), intent(in) :: fname
    integer :: phase_id
    type(string_array) :: dummy_string_array
    class(inelastic_strain_base), pointer :: inelastic_strain_base_ptr => null()
    class(porosity_base), pointer :: porosity_base_ptr => null()
    class(crystal_paremeters_type), pointer :: crystal_paremeters_ptr => null()
    logical :: use_crystal_plasticity = .false. , &
               use_isotropic_plasticity = .false., &
               use_glide =.false. , &
               use_climb =.false. , &
               use_diffusion =.false., &
               use_porosity = .false., &
               use_damage = .false.

    n_ss = 1
    n_gauss = 1
    call this%matf_reader%openReadTextFile(fname)
    
    

    call this%matf_reader%readLineAndCheckStringsAreEqual("--Phase-Parameters", dummy_string_array)
    call this%matf_reader%readParameter("phase-id", phase_id)
    if (this%phase_id.ne.phase_id) then
       write(*,*) "readFromFilePhaseMaterial: I'm reading the phase file for phase ", this%phase_id, " but the provided file contains information for phase number ", phase_id
      error stop " . Abort!"
    end if

  ! sanity checks
    call this%matf_reader%readParameter("use-isotropic-plasticity", use_isotropic_plasticity)
    call this%matf_reader%readParameter("use-crystal-plasticity", use_crystal_plasticity)
    call this%matf_reader%readParameter("use-glide",  use_glide)
    call this%matf_reader%readParameter("use-climb",  use_climb)
    call this%matf_reader%readParameter("use-diffusion",  use_diffusion)
    call this%matf_reader%readParameter("use-porosity",  use_porosity)
    call this%matf_reader%readParameter("use-damage",  use_damage)
    call this%matf_reader%skipEmptyLine()

  ! a few sanity checks before continuing
    if (use_isotropic_plasticity.and.use_crystal_plasticity) error stop "you can't use both isotropic and crystal palsticity"
    if (use_glide.and.(.not.use_crystal_plasticity)) error stop "you can't add glide without using crystal plasticity"
    if (use_climb.and.(.not.use_crystal_plasticity)) error stop "you can't add diffusion without using crystal plasticity"
    if (use_climb.and.(.not.use_glide)) error stop "you can't add climb without adding glide"
    if (use_damage.and.(.not.(use_porosity))) error stop "you cannot add damage without adding porosity"


    this%has_damage = use_damage
    call readElasticityFromFile(this%matf_reader, this%phase_id, this%all_mighty_grid, &
    this%sim_all_macro_data, this%the_bc_object, this%stiffness_ptr, this%common_material_parameter_ptr)

    if (use_isotropic_plasticity.or.use_crystal_plasticity) then
      call this%matf_reader%skipEmptyLine()
      allocate(this%common_material_parameter_ptr)
      call this%common_material_parameter_ptr%readCommonMaterialParametersFromFile(this%matf_reader)
      call this%matf_reader%skipEmptyLine()

      if (use_crystal_plasticity) then
        allocate(crystal_paremeters_ptr)
        crystal_paremeters_ptr%use_damage = this%has_damage
        call crystal_paremeters_ptr%readCrystalParameters(this%matf_reader)
        call this%matf_reader%skipEmptyLine()
        
        n_ss = crystal_paremeters_ptr%n_ss_total
        n_gauss = crystal_paremeters_ptr%n_gauss
          if (use_glide) then
            call readMaterialParametersFromFileCPGlide(this%matf_reader, this%phase_id, & 
                      this%has_damage, crystal_paremeters_ptr%n_gauss, crystal_paremeters_ptr%n_std_dev_gauss_integration, & 
                      this%all_mighty_grid, this%sim_all_macro_data, this%the_bc_object, this%stiffness_ptr, &
                      this%common_material_parameter_ptr, crystal_paremeters_ptr,&
                      inelastic_strain_base_ptr)
            call this%addDeformationMode(inelastic_strain_base_ptr)
            nullify(inelastic_strain_base_ptr)
            call this%matf_reader%skipEmptyLine()
          endif

          if (use_climb) then
            call readMaterialParametersFromFileCPClimb(this%matf_reader, this%phase_id, & 
                      this%has_damage, crystal_paremeters_ptr%n_gauss, crystal_paremeters_ptr%n_std_dev_gauss_integration, & 
                      this%all_mighty_grid, this%sim_all_macro_data, this%the_bc_object, this%stiffness_ptr, &
                      this%common_material_parameter_ptr, crystal_paremeters_ptr,&
                      inelastic_strain_base_ptr)
            call this%addDeformationMode(inelastic_strain_base_ptr)
            nullify(inelastic_strain_base_ptr)
            call this%matf_reader%skipEmptyLine()
          endif
      endif
    endif


    if (use_diffusion)then
      call readMaterialParametersFromFileISODiffusion(this%matf_reader, this%phase_id, this%has_damage, this%all_mighty_grid, &
                      this%sim_all_macro_data, this%the_bc_object, this%common_material_parameter_ptr , &
                      inelastic_strain_base_ptr)
      call this%addDeformationMode(inelastic_strain_base_ptr)
      nullify(inelastic_strain_base_ptr)
      call this%matf_reader%skipEmptyLine()
    endif


    if (use_porosity)then
      call readMaterialParametersFromFilePorosity(this%matf_reader, this%phase_id, this%all_mighty_grid, &
                                                  this%sim_all_macro_data, this%the_bc_object, this%common_material_parameter_ptr , &
                                                  porosity_base_ptr)
      call this%setPorosity(porosity_base_ptr)
      nullify(porosity_base_ptr)
    endif

    call this%matf_reader%closeTextFile()
  end subroutine


end module
