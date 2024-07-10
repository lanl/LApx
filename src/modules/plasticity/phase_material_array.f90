module phase_material_array_mod
#include "macro_debug.fpp"
use kinds
use inelastic_strain_mod, only : inelastic_strain_base, inelastic_strain_array
use polymorphic_dtype_array_mod, only : dtype_array_ptr
use all_grid_data_mod, only : all_grid_data
use all_mighty_grid_mod, only : all_mighty_grid_type, Matrix66MultiPhase
use simumaltion_macro_data_mod, only : sim_all_macro_data_obj
use stiffness_base_mod
use bc_objects_mod, only : boundary_condition_array_type
use csv_writer_mod, only : csv_writer
use phase_material_mod, only : phase_material_type
implicit none



type, extends(dtype_array_ptr) :: phase_material_array_type
    ! simulation parameter
  logical :: macro_object_linked =.false.
  class(all_mighty_grid_type), pointer :: all_mighty_grid => null()
  class(sim_all_macro_data_obj), pointer :: sim_all_macro_data => null()
  class(all_grid_data), pointer :: common_grid_data => null()
  class(boundary_condition_array_type), pointer :: the_bc_object => null()

  real(k_real), pointer :: dt => null(), &
                           dt_max => null(), &
                           temperature => null(), &
                           time => null()
  
  logical :: write_header = .true.
  real(k_real), dimension(:,:,:,:), pointer :: phase_fraction_grid => null()
  real(k_real), dimension(:,:,:,:,:), pointer :: average_stiffness_grid => null(), &
                                               elastic_strain => null(), &
                                               elastic_strain_old => null(), &
                                               elastic_strain_rate => null(), &
                                               plastic_strain => null(), &
                                               plastic_strain_old => null(), &
                                               plastic_strain_rate => null(), &
                                               total_strain => null(), &
                                               total_strain_old => null(), &
                                               total_strain_rate => null(), &
                                               cauchy_stress => null(), &
                                               cauchy_stress_old => null(), &
                                               cauchy_stress_rate => null(), &
                                               test_generic_matrix => null()
  
  type(Matrix66MultiPhase) :: multi_phase_stiffness
  type(csv_writer) :: simulation_averages_writer
  integer :: num_phases
contains
  procedure :: linkBaseObjects => linkBaseObjectsPhaseMatArray
  procedure :: addElement => addPhaseMaterial
  procedure :: getElementPtr => getPhaseMaterialPtr
  procedure :: initGridPointers => initGridPointersPhaseMatArray
  procedure :: addFieldVariablesToGrid => addFieldVariablesToGridPhaseMatArray
  procedure :: initStateVariables => initStateVariablesPhaseMaterialArray
  procedure :: updateStateVariablesInnerLoop => updateStateVariablesInnerLoopPhaseMatArray
  procedure :: updateStateVariablesOuterLoop => updateStateVariablesOuterLoopPhaseMatArray
  procedure :: updateStateVariablesStaggered => updateStateVariablesStaggeredPhaseMatArray
  procedure :: initTexture => initTexturePhaseMatArray
  procedure :: getStiffnessArray
  procedure :: updatePhaseParameters => updatePhaseParametersPhaseMatArray
  procedure :: acceptRejectSolution => acceptRejectSolutionPhaseMatArray
  procedure :: writeInelastiStrainRateToFile => writeInelastiStrainRateToFilePhaseMatArray
  procedure :: writeAverageQuantitiesCSV => writeAverageQuantitiesCSVPhaseMaterialArray
  procedure :: getConstitutiveStrainRate => getConstitutiveStrainRatePhaseMatArray
  procedure :: updateRotationMatrix => updateRotationMatrixPhaseMatArray
  procedure :: updateStiffness => updateStiffnessPhaseMatArray
  procedure :: updateAverageStiffness => updateAvgStiffnessPhaseMatArray
  procedure :: readFromFile => readFromFilePhaseMatArray
  end type

contains

  !------------------------------------------------------------------------------!
  !                  PHASE MATERIAL ARRAY SUBROUTINE                          !
  !------------------------------------------------------------------------------!
  subroutine addPhaseMaterial(this, new_element)
    class(phase_material_array_type), intent(inout) :: this
    class(phase_material_type), pointer, intent(inout) :: new_element

    call this%extend()
    this%all_pt(this%n_elements)%pt => new_element
    nullify(new_element)
    this%num_phases = this%n_elements
  end subroutine

  subroutine getPhaseMaterialPtr(this, idx, ptr)
    class(phase_material_array_type), intent(in) :: this
    integer, intent(in) :: idx
    class(phase_material_type), pointer, intent(out) :: ptr

    call this%checkElementExist(idx)

    associate(elem => this%all_pt(idx)%pt)
    select type(elem)
      class is (phase_material_type)
        ptr => elem
      class default
        error stop "getPhaseMaterialPtr wrong class"
    end select
    end associate
  end subroutine

  subroutine initGridPointersPhaseMatArray(this)
    implicit none
    class(phase_material_array_type), intent(inout) :: this
    class(phase_material_type), pointer :: ph_mat => null()
    integer :: i

    if (.not.this%macro_object_linked ) error stop "initGridPointersPhaseMatArray: macro_object_linked have not been linked"

    call this%common_grid_data%getGenericVectorDataPointerByName("phase_fraction", this%phase_fraction_grid)
    call this%common_grid_data%getTensor2DataPointerByName("elastic_strain", this%elastic_strain)
    call this%common_grid_data%getTensor2DataPointerByNameOld("elastic_strain", this%elastic_strain_old)
    call this%common_grid_data%getTensor2DataPointerByName("elastic_strain_rate", this%elastic_strain_rate)
    call this%common_grid_data%getTensor2DataPointerByName("plastic_strain", this%plastic_strain)
    call this%common_grid_data%getTensor2DataPointerByNameOld("plastic_strain", this%plastic_strain_old)
    call this%common_grid_data%getTensor2DataPointerByName("plastic_strain_rate", this%plastic_strain_rate)
    call this%common_grid_data%getTensor2DataPointerByName("total_strain", this%total_strain)
    call this%common_grid_data%getTensor2DataPointerByNameOld("total_strain", this%total_strain_old)
    call this%common_grid_data%getTensor2DataPointerByName("total_strain_rate", this%total_strain_rate)
    call this%common_grid_data%getTensor2DataPointerByName("cauchy_stress", this%cauchy_stress)
    call this%common_grid_data%getTensor2DataPointerByNameOld("cauchy_stress", this%cauchy_stress_old)
    call this%common_grid_data%getTensor2DataPointerByName("cauchy_stress_rate", this%cauchy_stress_rate)
    call this%common_grid_data%getMatrix66DataPointerByName("stiffness", this%average_stiffness_grid)
    call this%common_grid_data%getGenericMatrixDataPointerByName("test_generic_matrix", this%test_generic_matrix)

    call this%all_mighty_grid%AMGGetAllPhaseMatrix66("stiffness", this%multi_phase_stiffness)
    do i=1,this%n_elements
      call this%getElementPtr(i, ph_mat)
      call ph_mat%initGridPointers()
    enddo
  end subroutine

  subroutine addFieldVariablesToGridPhaseMatArray(this)
    use grid_data_var_type_mod
    implicit none
    class(phase_material_array_type), intent(inout) :: this
    class(phase_material_type), pointer :: ph_mat => null()
    integer :: i

    associate (common_grid_data => this%common_grid_data)
      call common_grid_data%addVar("elastic_strain", tensor2)
      call common_grid_data%addVar("elastic_strain_rate", tensor2)
      call common_grid_data%addVar("plastic_strain", tensor2)
      call common_grid_data%addVar("plastic_strain_rate", tensor2)
      call common_grid_data%addVar("total_strain", tensor2)
      call common_grid_data%addVar("total_strain_rate", tensor2)
      call common_grid_data%addVar("cauchy_stress", tensor2)
      call common_grid_data%addVar("cauchy_stress_rate", tensor2)
      call common_grid_data%addVar("stiffness", matrix66)
      call common_grid_data%addVar("test_generic_matrix", generic_matrix, additional_var_dimensions=(/3,2/))
    end associate

    do i=1,this%n_elements
      call this%getElementPtr(i, ph_mat)
      call ph_mat%addFieldVariablesToGrid()
    enddo
  end subroutine

  subroutine initTexturePhaseMatArray(this)
    implicit none
    class(phase_material_array_type), intent(inout) :: this
    class(phase_material_type), pointer :: ph_mat => null()
    integer :: i
    do i=1,this%n_elements
      call this%getElementPtr(i, ph_mat)
      call ph_mat%initTexture()
    enddo
  end subroutine

  subroutine initStateVariablesPhaseMaterialArray(this)
    implicit none
    class(phase_material_array_type), intent(inout) :: this
    class(phase_material_type), pointer :: ph_mat => null()
    integer :: i,j
    do,j=1,2
      do i=1,3
    this%test_generic_matrix(i,j,:,:,:) = i*10+j
      enddo
    enddo
    do i=1,this%n_elements
      call this%getElementPtr(i, ph_mat)
      call ph_mat%initStateVariables()
    enddo
  end subroutine

  subroutine updateStateVariablesInnerLoopPhaseMatArray(this, ix_in, iy_in, iz_in)
    implicit none
    class(phase_material_array_type), intent(inout) :: this
    integer, intent(in) :: ix_in, iy_in, iz_in
    class(phase_material_type), pointer :: ph_mat => null()
    integer :: i
    do i=1,this%n_elements
      call this%getElementPtr(i, ph_mat)
      call ph_mat%updateStateVariablesInnerLoop(ix_in, iy_in, iz_in)
    enddo
  end subroutine

  subroutine updateStateVariablesOuterLoopPhaseMatArray(this)
    implicit none
    class(phase_material_array_type), intent(inout) :: this
    class(phase_material_type), pointer :: ph_mat => null()
    integer :: i

    ! updating varaibles that are going top be present for each phase
    this%cauchy_stress = this%cauchy_stress_rate * this%dt + this%cauchy_stress_old
    this%elastic_strain = this%elastic_strain_rate * this%dt + this%elastic_strain_old
    this%plastic_strain = this%plastic_strain_rate * this%dt + this%plastic_strain_old
    this%total_strain = this%total_strain_rate * this%dt + this%total_strain_old

    do i=1,this%n_elements
      call this%getElementPtr(i, ph_mat)
      call ph_mat%updateStateVariablesOuterLoop()
    enddo
  end subroutine

  subroutine updateStateVariablesStaggeredPhaseMatArray(this)
    implicit none
    class(phase_material_array_type), intent(inout) :: this
    class(phase_material_type), pointer :: ph_mat => null()
    integer :: i
    do i=1,this%n_elements
      call this%getElementPtr(i, ph_mat)
      call ph_mat%updateStateVariablesStaggered()
    enddo
  end subroutine

  subroutine getStiffnessArray(this, stiffness_array)
    implicit none
    class(phase_material_array_type), intent(in) :: this
    class(phase_material_type), pointer :: ph_mat => null()
    real(k_real), dimension(:,:,:,:,:), pointer :: stiffness_array
    integer :: i
    if (associated(stiffness_array)) error stop "getStiffnessArray stiffness_array already allocated. Abort."
    allocate(stiffness_array(3,3,3,3,this%n_elements))
    do i=1,this%n_elements
      stiffness_array(:,:,:,:,i) = 0._k_real
      call this%getElementPtr(i, ph_mat)
      stiffness_array(:,:,:,:,i) = ph_mat%stiffness_ptr%stiffness_crystal_axes
    enddo
  end subroutine



  subroutine acceptRejectSolutionPhaseMatArray(this, dt_max, accept_solution_flag)
    implicit none
    class(phase_material_array_type), intent(inout) :: this
    real(k_real), intent(out) :: dt_max
    logical, intent(out) :: accept_solution_flag
    logical :: accept_solution_phase
    class(phase_material_type), pointer :: ph_mat => null()
    integer :: i

    accept_solution_flag = .true.
    do i=1,this%n_elements
      call this%getElementPtr(i, ph_mat)
      call ph_mat%acceptRejectSolution(dt_max, accept_solution_phase)
      accept_solution_flag = accept_solution_flag.and.accept_solution_phase
    enddo
  end subroutine

  subroutine updatePhaseParametersPhaseMatArray(this)
    implicit none
    class(phase_material_array_type), intent(inout) :: this
    class(phase_material_type), pointer :: ph_mat => null()
    integer :: i
    do i=1,this%n_elements
      call this%getElementPtr(i, ph_mat)
      if (.not.associated(ph_mat%stiffness_ptr)) error stop "updatePhaseParametersPhaseMatArray ph_mat%stiffness_ptr not associated!"
      call ph_mat%stiffness_ptr%updatePhaseParameters()
      call ph_mat%updatePhaseParameters()
    enddo
  end subroutine

  subroutine writeInelastiStrainRateToFilePhaseMatArray(this)
    use mpi_useful_routines_mod, only : MPIAverageGridMatrix
    implicit none
    class(phase_material_array_type), intent(inout) :: this
    class(phase_material_type), pointer :: ph_mat => null()
    real(k_real) :: avg_tensor(3,3)
    integer :: idx, num_write

    num_write = 1
    if (this%write_header) num_write = num_write +1
  
    do idx=1,num_write
      call this%simulation_averages_writer%openToAppendCSV()
      call this%simulation_averages_writer%AppendScalar(this%time, "time", this%write_header)
      call this%simulation_averages_writer%AppendScalar(this%temperature, "temperature", this%write_header)

      call MPIAverageGridMatrix(this%cauchy_stress, avg_tensor)
      call this%simulation_averages_writer%AppendTensor(avg_tensor, "stress", "cauchy_stress", this%write_header)
      call MPIAverageGridMatrix(this%cauchy_stress_rate, avg_tensor)
      call this%simulation_averages_writer%AppendTensor(avg_tensor, "stress", "cauchy_stress_rate", this%write_header)

      call MPIAverageGridMatrix(this%total_strain, avg_tensor)
      call this%simulation_averages_writer%AppendTensor(avg_tensor, "strain", "total_strain", this%write_header)
      call MPIAverageGridMatrix(this%total_strain_rate, avg_tensor)
      call this%simulation_averages_writer%AppendTensor(avg_tensor, "strain", "total_strain_rate", this%write_header)
      
      call MPIAverageGridMatrix(this%elastic_strain, avg_tensor)
      call this%simulation_averages_writer%AppendTensor(avg_tensor, "strain", "elastic_strain", this%write_header)
      call MPIAverageGridMatrix(this%elastic_strain_rate, avg_tensor)
      call this%simulation_averages_writer%AppendTensor(avg_tensor, "strain", "elastic_strain_rate", this%write_header)

      call MPIAverageGridMatrix(this%plastic_strain, avg_tensor)
      call this%simulation_averages_writer%AppendTensor(avg_tensor, "strain", "plastic_strain", this%write_header)
      call MPIAverageGridMatrix(this%plastic_strain_rate, avg_tensor)
      call this%simulation_averages_writer%AppendTensor(avg_tensor, "strain", "plastic_strain_rate", this%write_header)
      call this%simulation_averages_writer%closeCSVFile()
      this%write_header = .false.
    enddo

    do idx=1,this%n_elements
      call this%getElementPtr(idx, ph_mat)
      call ph_mat%writeInelastiStrainRateToFile()
    enddo

  end subroutine

  subroutine writeAverageQuantitiesCSVPhaseMaterialArray(this, write_headers)
    implicit none
    class(phase_material_array_type), intent(inout) :: this
    class(phase_material_type), pointer :: ph_mat => null()
    logical, intent(in) :: write_headers
    integer :: i
    do i=1,this%n_elements
      call this%getElementPtr(i, ph_mat)
      call ph_mat%writeAverageQuantitiesCSV(write_headers)
    enddo
  end subroutine

  subroutine linkBaseObjectsPhaseMatArray(this, all_mighty_grid_in, sim_all_macro_data, the_bc_object)
    implicit none
    class(phase_material_array_type), intent(inout) :: this
    class(all_mighty_grid_type), intent(in), pointer :: all_mighty_grid_in
    class(sim_all_macro_data_obj), intent(in), pointer :: sim_all_macro_data
    class(boundary_condition_array_type), intent(in), pointer :: the_bc_object
  
    if (this%macro_object_linked) error stop "linkBaseObjects already intialized. abort! "
    
    this%all_mighty_grid => all_mighty_grid_in
    call this%all_mighty_grid%AMGGetCommonGridPtr(this%common_grid_data)
    this%sim_all_macro_data => sim_all_macro_data
    this%the_bc_object => the_bc_object

    call the_bc_object%getTemperaturePointer(this%temperature)
    call sim_all_macro_data%sim_time%getTimePointer(this%time)
    call sim_all_macro_data%sim_time%getDeltaTimePointer(this%dt)
    call sim_all_macro_data%sim_time%getDeltaTMaxPointer(this%dt_max)

    this%macro_object_linked =.true.

    call this%simulation_averages_writer%createNewCSV("simulation_macro_averages.csv")
    
  end subroutine


  subroutine getConstitutiveStrainRatePhaseMatArray(this, stress6, epsilon_dot_total, depsilon_dot_total_dstress, epsilon_dot_el, epsilon_dot_inel, ix, iy, iz)
    use change_tensor_basis, only :chg_basis_vector6_to_tensor2, chg_basis_vector6_to_tensor2_fun
    implicit none
    class(phase_material_array_type), intent(inout) :: this
    real(k_real), target, intent(in) ::stress6(6)
    real(k_real), dimension(6), intent(out) :: epsilon_dot_total, epsilon_dot_el, epsilon_dot_inel
    real(k_real), dimension(6,6), intent(out) :: depsilon_dot_total_dstress
    integer, intent(in) :: ix, iy, iz
    real(k_real), dimension(6) :: epsilon_dot_inel_ph, epsilon_dot_el_ph, epsilon_dot_total_ph
    real(k_real), dimension(6,6) :: depsilon_dot_dstress_ph
    integer :: iph
    class(phase_material_type), pointer :: ph_mat_ptr =>null()
    real(k_real), pointer, dimension(:,:) :: phase_stiffness

    epsilon_dot_total =0._k_real
    depsilon_dot_total_dstress = 0._k_real
    epsilon_dot_inel  = 0._k_real

    ! call homegenization
    do iph=1,this%n_elements
      associate(ph_weight => this%phase_fraction_grid(iph, ix, iy, iz))

      phase_stiffness => this%multi_phase_stiffness%ph(iph)%data(:,:,ix,iy,iz)

      call this%getElementPtr(iph, ph_mat_ptr)
      call ph_mat_ptr%getConstitutiveStrainRate(stress6, epsilon_dot_total_ph, depsilon_dot_dstress_ph, epsilon_dot_el_ph, epsilon_dot_inel_ph, ix, iy, iz)
      

     

      epsilon_dot_el = epsilon_dot_el + epsilon_dot_el_ph * ph_weight
      epsilon_dot_inel = epsilon_dot_inel+ epsilon_dot_inel_ph * ph_weight
      epsilon_dot_total =  epsilon_dot_total + epsilon_dot_total_ph * ph_weight
      depsilon_dot_total_dstress = depsilon_dot_total_dstress + depsilon_dot_dstress_ph * ph_weight

      nullify(ph_mat_ptr)
      end associate
    enddo

    call chg_basis_vector6_to_tensor2((chg_basis_vector6_to_tensor2_fun(stress6) - this%cauchy_stress_old(:,:, ix, iy, iz))/this%dt, this%cauchy_stress_rate(:,:,ix,iy,iz))
    call chg_basis_vector6_to_tensor2(epsilon_dot_el, this%elastic_strain_rate(:,:,ix,iy,iz) )
    call chg_basis_vector6_to_tensor2(epsilon_dot_inel, this%plastic_strain_rate(:,:,ix,iy,iz) )
    call chg_basis_vector6_to_tensor2(epsilon_dot_total, this%total_strain_rate(:,:,ix,iy,iz) )
    
  end subroutine

  subroutine updateRotationMatrixPhaseMatArray(this)
    implicit none
    class(phase_material_array_type), intent(inout) :: this
    class(phase_material_type), pointer :: ph_mat_ptr =>null()
    integer :: i

    do i=1,this%n_elements
      call this%getElementPtr(i, ph_mat_ptr)
      call ph_mat_ptr%updateRotationMatrix()
      nullify(ph_mat_ptr)
    enddo
  end subroutine

  subroutine updateStiffnessPhaseMatArray(this)
    implicit none
    class(phase_material_array_type), intent(inout) :: this
    class(phase_material_type), pointer :: ph_mat_ptr =>null()
    integer :: i
    integer :: x_s, x_e, y_s, y_e, z_s, z_e, ix, iy, iz
    real(k_real), pointer :: phase_stiffness_grid(:,:,:,:,:) => null()
    class(all_grid_data), pointer :: phase_grid_ptr => null()
    do i=1,this%n_elements
      call this%getElementPtr(i, ph_mat_ptr)
      call ph_mat_ptr%updateStiffness()
      nullify(ph_mat_ptr)
    enddo

    ! get pointer

    call this%all_mighty_grid%AMGGetLoopLimitsRank(x_s, x_e, y_s, y_e, z_s, z_e)

    this%average_stiffness_grid = 0._k_real
    do i=1,this%n_elements
      call this%all_mighty_grid%AMGGetPhaseGridPtrByIdx(i, phase_grid_ptr )
      call phase_grid_ptr%getMatrix66DataPointerByName("stiffness", phase_stiffness_grid)

      do  iz = z_s,z_e
        do  iy = y_s,y_e
          do  ix = x_s, x_e
            this%average_stiffness_grid(:,:,ix,iy,iz) = this%average_stiffness_grid(:,:,ix,iy,iz) + &
            phase_stiffness_grid(:,:,ix,iy,iz) * this%phase_fraction_grid(i,ix,iy,iz)
          enddo
        enddo
      enddo
      nullify(phase_stiffness_grid, phase_grid_ptr)
    enddo

  end subroutine

  subroutine updateAvgStiffnessPhaseMatArray(this, C0_BB, S0_BB, C0_T4, S0_T4)
    use mpi_variables_mod, only : i_am_mpi_master
    use mpi_useful_routines_mod, only : MPIAverageGridMatrix
    use change_tensor_basis, only : chg_basis_matrix66_to_tensor4
    use voigt_indicial_conversion_mod, only : Tensor4ToMatrixVoigt
    use print_utils_mod, only : printToScreen
    use matrix_inversion, only : matrixInverseSymmetric
    implicit none
    class(phase_material_array_type), intent(inout) :: this
    real(k_real), intent(out) :: C0_BB(6,6), S0_BB(6,6), &
                                 C0_T4(3,3,3,3), S0_T4(3,3,3,3)
    real(k_real) :: CO_Voigt(6,6)

    call MPIAverageGridMatrix(this%average_stiffness_grid, C0_BB)
    call chg_basis_matrix66_to_tensor4(C0_BB, C0_T4)
    CO_Voigt = Tensor4ToMatrixVoigt(C0_T4)

    if (i_am_mpi_master) call printToScreen(C0_BB, "average stiffness C066 NEW")
    if (i_am_mpi_master) call printToScreen(CO_Voigt, "average stiffness c0Voigt NEW")
    
  call matrixInverseSymmetric(C0_BB, S0_BB)
  call chg_basis_matrix66_to_tensor4(S0_BB, S0_T4)

  end subroutine

  subroutine readFromFilePhaseMatArray(this, microstrucutre_info, n_ss, n_gauss)
    use microstructure_info_mod
    implicit none
    class(phase_material_array_type), intent(inout) :: this
    class(phase_material_type), pointer :: ph_mat_ptr =>null()
    type(microstructure_info_type) :: microstrucutre_info
    integer :: iph, nph
    integer, pointer, dimension(:), intent(inout) :: n_ss, n_gauss 
    integer ::  n_ss_ph, n_gauss_ph
    nph = microstrucutre_info%getNumPhases()
    allocate(n_ss(nph), &
             n_gauss(nph))

    do iph=1,nph
        allocate(ph_mat_ptr)
        call ph_mat_ptr%linkBaseObjects(iph, this%all_mighty_grid, this%sim_all_macro_data, this%the_bc_object)
        call ph_mat_ptr%init()
        call ph_mat_ptr%readFromFile(microstrucutre_info%getPhaseFileNameByIndex(iph), n_ss_ph, n_gauss_ph)
        call this%addElement(ph_mat_ptr)
        nullify(ph_mat_ptr)
        n_ss(iph) = n_ss_ph
        n_gauss(iph) = n_gauss_ph
    enddo
  end subroutine
end module phase_material_array_mod
