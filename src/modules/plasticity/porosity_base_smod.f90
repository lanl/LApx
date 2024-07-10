submodule(porosity_base_mod) porosity_base_smod
use kinds
use mpi_variables_mod, only : i_am_mpi_master

#include "macro_debug.fpp"

implicit none

contains

  module subroutine readMaterialParametersFromFilePorosityBase(this, matf_reader)
    use read_from_file_utils, only : file_reader
    use string_module , only : string_array
    use units_conversion_mod, only : deg2rad

    class(porosity_base), intent(inout) :: this
    class(file_reader), intent(inout) :: matf_reader
    type(string_array) :: dummy_string_array
    real(k_real) :: theta

    call matf_reader%readLineAndCheckStringsAreEqual("-pore-nucleation-parameters", dummy_string_array)
    call matf_reader%readParameter("use-porosity-nucleation[TRUE/FALSE]", this%use_porosity_nucleation)
    call matf_reader%readParameter("use-preseeded-porosity[TRUE/FALSE]", this%use_preseeded_porosity)
    call matf_reader%readParameter("use-bulk-porosity-evolution[TRUE/FALSE]", this%use_bulk_porosity_evolution)
    call matf_reader%readParameter("initial-porosity", this%initial_porosity)
    
    call matf_reader%readParameter("pore-nucleation-standard-deviation[mm/mm]", this%nucleation_std_dev)
    call matf_reader%readParameter("critical-pore-nucleation-strain[mm/mm]", this%critical_nucleation_strain)
    call matf_reader%readParameter("max-nucleation-porosity-fraction[%]", this%max_porosity_nucleation)
    this%max_porosity_nucleation = this%max_porosity_nucleation/100._k_real
    call matf_reader%skipEmptyLine()

    if (this%use_porosity_nucleation.and.this%use_preseeded_porosity) error stop "Porosity module: both use nucleationa dn use preseed porosity are true"
    if (.not.(this%use_porosity_nucleation).and.(.not.(this%use_preseeded_porosity))) error stop "Porosity module: both use nucleationa dn use preseed porosity are false"
    if ((this%use_preseeded_porosity.and.(this%initial_porosity.le.0._k_real))) error stop "Porosity module: use preseed porosity TRUE, but initial porosity <= 0."
    if ((this%use_preseeded_porosity.and.(this%initial_porosity.gt.1._k_real))) error stop "Porosity module: use preseed porosity TRUE, but initial porosity > 1."

    call matf_reader%readLineAndCheckStringsAreEqual("-pore-growth-parameters", dummy_string_array)
    call matf_reader%readParameter("use-plasticity-constrained-growth[TRUE/FALSE]", this%use_plasticity_constrained_growth)
    call matf_reader%readParameter("use-diffusive-growth[TRUE/FALSE]", this%use_diffusive_growth)
    call matf_reader%readParameter("pore-half-tip-angle[deg]", theta)
    if (theta<=0._k_real) error stop "PorosityBase theta <=0."
    if (theta>90._k_real) error stop "PorosityBase theta > 90."
    theta = theta*deg2rad
    this%sin_theta=sin(theta)
    this%h = (1._k_real/(1._k_real+cos(theta)) - 0.5_k_real*cos(theta))/this%sin_theta
    call matf_reader%readParameter("min-allowed-void-size[m]", this%min_void_size)

    this%omega => this%common_material_parameter_ptr%atomic_volume
    this%E_formation => this%common_material_parameter_ptr%vacancy_formation_energy
    this%E_migration => this%common_material_parameter_ptr%vacancy_migration_energy
    this%delta_gb => this%common_material_parameter_ptr%gb_thickness
    this%vacancy_diffusivity_coefficient => this%common_material_parameter_ptr%vacancy_diffusivity_coefficient
    
    call matf_reader%readParameter("grain-boundary-surface-tension[J/m2]", this%surface_tension)
    this%surface_tension = this%surface_tension
    call matf_reader%skipEmptyLine()


    call matf_reader%readLineAndCheckStringsAreEqual("-void-coalescence-parameters", dummy_string_array)
    call matf_reader%readParameter("use-coalescence[TRUE/FALSE]", this%use_coalescence)
    call matf_reader%readParameter("coalescence-critical-porosity[fc,unitless]", this%coalescence_critical_porosity)
    call matf_reader%readParameter("saturation-porosity[fs,unitless]", this%saturation_porosity)
    call matf_reader%readParameter("effective-porosity-at-fracture[fu,unitless]", this%fustar)
    call matf_reader%skipEmptyLine()
    
    call matf_reader%readLineAndCheckStringsAreEqual("-time-march-tolerances", dummy_string_array)
    call matf_reader%readVectorParameter( "porosity-absolute-and-relative-tol[volume-fraction,unitless]", 2, this%porosity_abs_rel_tol)
  

    this%initial_void_radius = 3e-8_k_real
    call this%computeVoidVolume(this%initial_void_radius, this%initial_void_volume)
  end subroutine

  module subroutine addFieldVariablesToGridPorosityBase(this)
    use grid_data_var_type_mod
    implicit none
    class(porosity_base), intent(inout) :: this

    if (.not.(this%macro_object_linked)) error stop "addFieldVariablesToGridPorosityBase: you can't add field varaibles to the grid without first linking the global obejcts"
    if ( .not.(this%parameters_initialized)) error stop "addFieldVariablesToGridPorosityBase: you can't add field varaibles to the grid without first initializing a material parameters"

    call this%material_base%addFieldVariablesToGrid()
    associate (all_grid_data_vars => this%grid_data)
    call all_grid_data_vars%addVar("porosity", scalar, stateful_level=2)
    call all_grid_data_vars%addVar("porosity_rate", scalar, stateful_level=2)
    call all_grid_data_vars%addVar("porosity_nucleation", scalar, stateful_level=2)
    call all_grid_data_vars%addVar("effective_porosity", scalar, stateful_level=2)
    call all_grid_data_vars%addVar("effective_porosity_rate", scalar, stateful_level=2)
    call all_grid_data_vars%addVar("porosity_growth", scalar, stateful_level=2)
    call all_grid_data_vars%addVar("stress", tensor2)
    call all_grid_data_vars%addVar("void_radius", generic_vector, additional_var_dimensions=(/this%n_bins/), stateful_level=2)
    call all_grid_data_vars%addVar("void_density", generic_vector, additional_var_dimensions=(/this%n_bins/), stateful_level=2)

    end associate

    this%grid_variables_provided =.true.
  end subroutine

  module subroutine initParameters(this, phase_id, common_material_parameter_ptr)
    implicit none
    class(porosity_base), intent(inout) :: this
    integer, intent(in) :: phase_id
    class(common_material_parameter), pointer, intent(in) :: common_material_parameter_ptr


    if (this%parameters_initialized) error stop "you can initialize parameters only once"

    call this%initParametersMaterialBase(phase_id, common_material_parameter_ptr)

  end subroutine

  module subroutine initGridPointersPorosityBase(this)
    implicit none
    class(porosity_base), intent(inout) :: this

    if (.not.this%macro_object_linked) &
      error stop "initGridPointersPorosityBase you cannot initialized GridPointers without first linking base objects. Abort!"
    if (.not.this%grid_variables_provided) &
      error stop "initGridPointersPorosityBase you cannot initialized GridPointers without first providing the grid variables. Abort!"
    if (this%grid_pointers_linked) error stop "you can link grid pointers only once. Abort!"
    if (.not.(this%grid_data%initialized)) error stop "initGridPointersPorosityBase you cannot init grid pointer without first initializing the grid!"

    call this%material_base%initGridPointers()

    call this%grid_data%getTensor2DataPointerByName("stress", this%stress_grid)
    call this%grid_data%getTensor2DataPointerByName("plastic_strain", this%inelastic_strain_grid)
    call this%grid_data%getTensor2DataPointerByName("plastic_strain_rate", this%inelastic_strain_rate_grid)
    call this%grid_data%getScalarDataPointerByName("porosity", this%porosity_grid)
    call this%grid_data%getScalarDataPointerByNameOld("porosity", this%porosity_grid_old)
    call this%grid_data%getScalarDataPointerByName("porosity_rate", this%porosity_rate_grid)
    call this%grid_data%getScalarDataPointerByNameOld("porosity_rate", this%porosity_rate_grid_old)
    call this%grid_data%getScalarDataPointerByName("porosity_growth", this%porosity_growth_grid)
    call this%grid_data%getScalarDataPointerByNameOld("porosity_growth", this%porosity_growth_grid_old)
    call this%grid_data%getScalarDataPointerByName("porosity_nucleation", this%porosity_nucleation_grid)
    call this%grid_data%getScalarDataPointerByNameOld("porosity_nucleation", this%porosity_nucleation_grid_old)

    call this%grid_data%getScalarDataPointerByName("effective_porosity", this%effective_porosity_grid)
    call this%grid_data%getScalarDataPointerByNameOld("effective_porosity", this%effective_porosity_grid_old)
    call this%grid_data%getScalarDataPointerByName("effective_porosity_rate", this%effective_porosity_rate_grid)
    call this%grid_data%getScalarDataPointerByNameOld("effective_porosity_rate", this%effective_porosity_rate_grid_old)

    call this%common_grid_data%getRealVectorDataPointerByName("gb_normals", this%gb_normals_grid)
    call this%common_grid_data%getScalarIntegerDataPointerByName("gb_id", this%gb_id_grid)
    call this%grid_data%getGenericVectorDataPointerByName("void_radius", this%void_radius_grid)
    call this%grid_data%getGenericVectorDataPointerByName("void_density", this%void_density_grid)
    call this%grid_data%getGenericVectorDataPointerByNameOld("void_radius", this%void_radius_old_grid)
    call this%grid_data%getGenericVectorDataPointerByNameOld("void_density", this%void_density_old_grid)

    this%grid_pointers_linked =.true.
  end subroutine

module subroutine setPointersPorosityBase(this, ix, iy, iz)
  use change_tensor_basis, only : chg_basis_tensor2_to_vector6
  implicit none
  class(porosity_base), intent(inout) :: this
  integer, intent(in) ::  ix, iy, iz
  if (.not.(this%grid_data%initialized)) error stop "setPointersPorosityBase you cannot init grid pointer without first initializing the grid!"

  call this%material_base%setPointData(ix, iy, iz)
  this%stress_ptr => this%stress_grid(:,:, ix, iy, iz)
  this%inelastic_strain_rate_ptr => this%inelastic_strain_rate_grid(:, :, ix, iy, iz)
  this%inelastic_strain_ptr => this%inelastic_strain_grid(:,:, ix, iy, iz)
  this%gb_normals_ptr => this%gb_normals_grid(:, ix, iy, iz)
  this%gb_id_ptr => this%gb_id_grid( ix, iy, iz)
  this%void_radius_ptr => this%void_radius_grid( :, ix, iy, iz )
  this%void_density_ptr => this%void_density_grid( :, ix, iy, iz )
  this%void_radius_old_ptr => this%void_radius_old_grid( :, ix, iy, iz )
  this%void_density_old_ptr => this%void_density_old_grid( :, ix, iy, iz )

  call chg_basis_tensor2_to_vector6(this%stress_ptr, this%stress6)

  this%porosity_ptr => this%porosity_grid(ix, iy, iz)
  this%porosity_ptr_old  => this%porosity_grid_old(ix, iy, iz)
  this%porosity_rate_ptr => this%porosity_rate_grid(ix, iy, iz)
  this%porosity_rate_ptr_old  => this%porosity_rate_grid_old(ix, iy, iz)
  this%porosity_nucleation_ptr => this%porosity_nucleation_grid(ix, iy, iz)
  this%porosity_nucleation_ptr_old  => this%porosity_nucleation_grid_old(ix, iy, iz)
  this%porosity_growth_ptr => this%porosity_growth_grid(ix, iy, iz)
  this%porosity_growth_ptr_old  => this%porosity_growth_grid_old(ix, iy, iz)
  this%effective_porosity_ptr  => this%effective_porosity_grid(ix, iy, iz)
  this%effective_porosity_ptr_old  => this%effective_porosity_grid_old(ix, iy, iz)
  this%effective_porosity_rate_ptr => this%effective_porosity_rate_grid(ix, iy, iz)
  this%effective_porosity_rate_ptr_old  => this%effective_porosity_rate_grid_old(ix, iy, iz)
end subroutine

module subroutine updatePorosityBase(this)
  use tensor_math_mod, only : computeVMEquivalentStrain, computeTrace, computeVMEquivalentStress, computePressure
  use change_tensor_basis, only : chg_basis_vector6_to_tensor2
  use math_constants, only : PI
  use test_utils_mod, only : test_is_finite
  use print_utils_mod, only : printToScreen
  use units_conversion_mod, only : MPA2Pa
  class(porosity_base), intent(inout) :: this
  ! real(k_real) :: edotp_VM, ep_VM, edotp_H, s_VM, &
  !                 Tn, D, kc, fadj, da_needleman

  ! integer :: b_idx, b_idx2, new_void_bin_idx
  ! real(k_real) :: a_new_temp(this%n_bins), total_density(this%n_bins), new_size(this%n_bins), &
  !                 dn_nucleation, n_old(this%n_bins), &
  !                 df_growth, df_growth_gurson, &
  !                 da, df_nucleation, L, lambda_f, void_sintering_stress
  ! integer :: new_bin_idx(this%n_bins)
  logical  :: gb_voxel

  REAL(k_real) :: D, df_growth, df_growth_gurson, df_nucleation, dn_nucleation, edotp_H, edotp_VM, ep_VM, KC, s_VM, TN
  ! there is a bug in the GB normal calculation and some points get 0. normal
  ! even if the point is classified as GB. We shut down those points by skipping all points
  ! with a zero normal
  
  associate(edotp33 => this%inelastic_strain_rate_ptr, &
            ep33 => this%inelastic_strain_ptr, &
            a => this%void_radius_ptr, &
            a_old => this%void_radius_old_ptr, &
            n => this%void_density_ptr, &
            n_old => this%void_density_old_ptr, &
            fc => this%coalescence_critical_porosity, &
            fs =>this%saturation_porosity, &
            fustar => this%fustar, &
            f => this%porosity_ptr, &
            f_old => this%porosity_ptr_old, &
            f_nucleation => this%porosity_nucleation_ptr, &
            f_nucleation_old => this%porosity_nucleation_ptr_old, &
            f_growth => this%porosity_growth_ptr, &
            f_growth_old => this%porosity_growth_ptr_old, &
            fstar => this%effective_porosity_ptr, &
            fstar_old => this%effective_porosity_ptr_old, &
            h => this%h &
            )

  
  df_growth = 0._k_real
  df_growth_gurson = 0._k_real
  df_nucleation = 0._k_real
  f = f_old
  f_growth = f_growth_old
  f_nucleation = f_nucleation_old
  this%porosity_rate_ptr = 0._k_real
  this%effective_porosity_rate_ptr = 0._k_real
  
  
  gb_voxel = this%gb_id_ptr.gt.1.and.(any(this%gb_normals_ptr.ne.0._k_real))
  
  if (gb_voxel.or.this%use_bulk_porosity_evolution) then
  if (f_old.lt.fs) then ! check we are not at saturation porosity yet

#ifdef __DEBUG__
  call test_is_finite(Tn, "normal traction")
#endif
  
    s_VM = computeVMEquivalentStress(this%stress_ptr)*MPA2Pa
    ep_VM = computeVMEquivalentStrain(ep33)
    edotp_VM = computeVMEquivalentStrain(edotp33)
    edotp_H = computeTrace(edotp33)

    Tn = -1._k_real ! negative value to not trigger VoidNucleation
    if (gb_voxel) then
      Tn = sum((matmul(this%stress_ptr, this%gb_normals_ptr)* this%gb_normals_ptr))
    else if (this%use_bulk_porosity_evolution) then
      Tn = computePressure(this%stress_ptr)
    endif
    Tn = Tn * MPA2Pa

    ! growth
    df_growth = 0._k_real  
    df_growth_gurson = 0._k_real  
    if (f_old.gt.0._k_real) then
      if (this%use_plasticity_constrained_growth)  then 
        df_growth_gurson = (1._k_real-f)*edotp_H*this%dt
        df_growth = df_growth_gurson
      endif
      if (this%use_diffusive_growth.and.gb_voxel) then 
        call this%computeDiffusivityPorosityBase(D)
        call this%VoidGrowthGursonPlusDiffusion(a_old(1), n_old(1), D, s_VM, Tn,  edotp_VM, edotp_H, df_growth_gurson, df_growth)
      endif
    endif



    f_growth = f_growth_old + df_growth

    df_nucleation = 0._k_real
    dn_nucleation = 0._k_real
    if(this%use_porosity_nucleation) then
      ! compute nucelation stress
      if (this%use_diffusive_growth.and.gb_voxel) then
        call this%VoidNucleationChuNeedlemannNSites(Tn, this%critical_nucleation_strain, this%nucleation_std_dev, ep_VM, edotp_VM, dn_nucleation, df_nucleation)
        df_nucleation = max(df_nucleation, 0._k_real)
        dn_nucleation = max(dn_nucleation, 0._k_real)
        f_nucleation = f_nucleation_old + df_nucleation
        n(1) = n_old(1) + dn_nucleation
      else 
        call this%VoidNucleationChuNeedlemann(Tn, this%critical_nucleation_strain, this%nucleation_std_dev, ep_VM, edotp_VM, df_nucleation)

        df_nucleation = max(df_nucleation, 0._k_real)
        f_nucleation = f_nucleation_old + df_nucleation
        ! if (this%ix*this%iy*this%iz.eq.1.and.i_am_mpi_master) then
        !   write(*,*) "stress"
        !   write(*,*) this%stress_ptr
        !   write(*,*) "inelastic_strain"
        !   write(*,*) this%inelastic_strain_ptr
        !   write(*,*) "inelastic_strain_rate"
        !   write(*,*) this%inelastic_strain_rate_ptr
        !   write(*,*) this%critical_nucleation_strain, this%nucleation_std_dev, ep_VM, edotp_VM, df_nucleation
        !   write(*,*) "this%critical_nucleation_strain, this%nucleation_std_dev, ep_VM, edotp_VM, df_nucleation"
        !   write(*,*) this%critical_nucleation_strain, this%nucleation_std_dev, ep_VM, edotp_VM, df_nucleation
        !   write(*,*) "df_nucleation, f_nucleation"
        !   write(*,*) df_nucleation, f_nucleation
        ! endif
      endif
    endif
    
    ! update current f 
    f = max(f_old+df_nucleation+df_growth,0._k_real) ! f cannot go negative
    f = min(f,fs) ! limit f to fs                    ! f cannot be > f saturation
    this%porosity_rate_ptr = (f-f_old)/this%dt
    this%effective_porosity_rate_ptr = this%porosity_rate_ptr
    fstar = f ! effective porosity

    ! homegnize void size
    if (this%use_diffusive_growth) then
      a(1) = (f/(n(1)*4._k_real/3._k_real*PI*this%h))**(1._k_real/3._k_real)
    endif

  ! **************************************
  ! void coalescence
  ! **************************************
  if (this%use_coalescence) then
    if (fstar>fc) then ! once we get to the saturation porosity we stop
      kc = (fustar-fc)/(fs-fc)
      fstar = fc + kc*(f-fc)
    endif
    fstar = min(fstar,fustar) ! upper limit is fustar
    this%effective_porosity_rate_ptr = (fstar-fstar_old)/this%dt
  endif


#ifdef __DEBUG__
  call test_is_finite(fstar, "fstar")
#endif

endif !f_old.lt.fs
endif ! gb_voxel.or.this%use_bulk_porosity_evolution
  if (f.gt.0._k_real) this%total_porosity = this%total_porosity + fstar

  end associate



end subroutine

module subroutine computeDiffusivityPorosityBase(this, D)
  use math_constants, only : kBSI
  class(porosity_base), intent(inout) :: this
  real(k_real), intent(out) :: D
  real(k_real) :: D0
  associate(delta_gb=>this%delta_gb, &
            omega => this%omega, &
            T => this%temperature, &
            Diff => this%vacancy_diffusivity_coefficient(min(this%gb_id_ptr,2)), &
            E_migration => this%E_migration(min(this%gb_id_ptr,2)), & 
            E_formation => this%E_formation(min(this%gb_id_ptr,2)) & 
            )

  D0 = Diff*exp(-E_migration/(kBSI* T))

  D = D0 * delta_gb * omega/(kBSI* T)*&
      exp(-E_formation/(kBSI* T))
  end associate
end subroutine

module subroutine updateStateVarsAtMaterialPointStaggeredPorosityBase(this)
  implicit none
  class(porosity_base), intent(inout) :: this
  __DECL_CLASS_UNUSED_THIS__
  __SUPPRESS_CLASS_UNUSED_THIS__
end subroutine

module subroutine updateStateVarsAtMaterialPointOuterLoopPorosityBase(this)
  implicit none
  class(porosity_base), intent(inout) :: this

  call this%updatePorosity()
end subroutine

module subroutine VoidNucleationChuNeedlemann(this, Tn, e_critical, e_std_dev, eq_plastic_strain, eq_plastic_strain_rate, df)
  use probability_mod, only : gaussianProbabilityFromMeanAndStdDev
  class(porosity_base), intent(inout) :: this
  real(k_real), intent(in) :: Tn, e_critical, e_std_dev, eq_plastic_strain, eq_plastic_strain_rate
  real(k_real), intent(out) :: df
  real(k_real) :: F_nucleation, probability
  if (Tn>0._k_real) then
    call gaussianProbabilityFromMeanAndStdDev( e_critical, e_std_dev, eq_plastic_strain, probability)
    F_nucleation = this%max_porosity_nucleation * probability
    df = F_nucleation * eq_plastic_strain_rate * this%dt
  else
    df = 0._k_real
  endif

end subroutine

module subroutine VoidNucleationChuNeedlemannNSites(this, Tn, e_critical, e_std_dev, eq_plastic_strain, eq_plastic_strain_rate, dn, df)
  use probability_mod, only : gaussianProbabilityFromMeanAndStdDev
  class(porosity_base), intent(inout) :: this
  real(k_real), intent(in) :: Tn, e_critical, e_std_dev, eq_plastic_strain, eq_plastic_strain_rate
  real(k_real), intent(out) :: dn, df
  real(k_real) :: F_nucleation, probability


  if (Tn>0._k_real) then
    call gaussianProbabilityFromMeanAndStdDev( e_critical, e_std_dev, eq_plastic_strain, probability)
    F_nucleation = this%max_porosity_nucleation * probability
    df = F_nucleation * eq_plastic_strain_rate * this%dt
    dn = df / this%initial_void_volume
  else
    dn = 0._k_real
    df = 0._k_real
  endif

end subroutine

module subroutine VoidNucleationSize(this, Tn,  a0)
  class(porosity_base), intent(inout) :: this
  real(k_real), intent(in) :: Tn
  real(k_real), intent(out) :: a0

  associate( gamma_s => this%surface_tension, &
             sin_theta => this%sin_theta)

  if (Tn>0._k_real) then
    a0 = 2._k_real*gamma_s * this%sin_theta/(Tn)
  else
    a0= 0._k_real
  endif

  end associate

end subroutine

module subroutine VoidGrowthGursonPlusDiffusion(this, a_old, n_density, D, s_VM, Tn, edotp_VM, edotp_H, df_gurson, df_growth)
  use test_utils_mod, only : test_is_finite
  use print_utils_mod, only : printToScreen
  use math_constants, only : PI
  use, intrinsic :: IEEE_ARITHMETIC
  implicit none
  class(porosity_base), intent(inout) :: this
  real(k_real), intent(in) :: a_old, & !-> old void size
                              n_density, & !-> void number density
                              D, &     ! GB diffusivity
                              s_VM,  & ! von mises stress
                              Tn,  &   ! normal traction
                              edotp_VM, &  ! equivalent strain rate
                              edotp_H, & ! hidrostatic plastic strain rate
                              df_gurson ! df due to gurson
  real(k_real), intent(out) :: df_growth
  real(k_real) :: da_diff, &
                  da_vp,&
                  lambda_f, &
                  fadj, &
                  void_sintering_stress, &
                  L, &
                  Vdot, &
                  a_new
  __DECL_UNUSED_REAL__        
  __SUPPRESS_UNUSED_REAL__(edotp_H)
  __SUPPRESS_UNUSED_REAL_OUT__(Vdot)
  associate( a => a_old, &
             a_min => this%min_void_size, &
             gamma_s => this%surface_tension, &
             sin_theta => this%sin_theta, &
             f => 4._k_real/3._k_real * PI * this%h * a_old**3 * n_density, &
             h =>  this%h)

  

  da_vp = 0._k_real
  if (this%use_plasticity_constrained_growth) then
    ! let's compute the cahnge in radius for single void
    da_vp = df_gurson/(4._k_real * PI * h * a**2*n_density)
  endif
  
  da_diff = 0._k_real
  if (this%use_diffusive_growth.and.this%gb_id_ptr>1) then

    !! compute void interspacing
    lambda_f = 1._k_real/n_density**(1._k_real/3._k_real) ! void mean spacing

    !! compute adjusted porosity
    if (edotp_VM.ne.0._k_real) then
      L = ( D * s_VM/edotp_VM)**(1._k_real/3._k_real)
      fadj = max((a/lambda_f)**2, (a/(a+1.5_k_real*L))**2)
    else
      fadj = (a/lambda_f)**2
    endif

    if (fadj.le.0._k_real) error stop "fadj  is <=0."

#ifdef __DEBUG__
  call test_is_finite(fadj, "fadj")
#endif

    ! compute the sintering stress
    void_sintering_stress = 2._k_real*gamma_s * sin_theta/a
    da_diff = D / (a**2*h) * (Tn - (1._k_real-fadj)*void_sintering_stress)/( &
          log(1._k_real/fadj) - 0.5_k_real*(3._k_real-fadj)*(1._k_real-fadj)   ) * this%dt
  endif

!
#ifdef __DEBUG__
  call test_is_finite(da_diff, "da_diff")
#endif

  if (IEEE_IS_NAN(da_diff).or.(.not.IEEE_IS_FINITE(da_diff))) then
    call printToScreen(D, "D")
    call printToScreen(fadj, "fadj")
    call printToScreen(Tn, "Tn")
    call printToScreen(void_sintering_stress, "void_sintering_stress")
    call printToScreen(a, "a")
    error stop "something is not right"
  endif
  a_new = max(a_old + da_vp + da_diff, this%min_void_size)
  ! we compute deltaF as a change in void volume
  df_growth= 4._k_real/3._k_real*PI*h*(a_new**3-a**3)*n_density

  end associate
end subroutine

module subroutine FindVoidRadiusBin(this, a, bin_idx)
  class(porosity_base), intent(in) :: this
  real(k_real), intent(in) :: a
  integer, intent(out) :: bin_idx

  ! first we check if we are in the last bin
  if (a>this%void_radius_bin_limits(this%n_bins-1)) then
    bin_idx = this%n_bins
  else ! if we are not in the last bin then loop and search for the bin
    do bin_idx=1,(this%n_bins-1)
      if (a<this%void_radius_bin_limits(bin_idx)) exit
    enddo
  endif

end subroutine

module subroutine initStateVariablesAtMaterialPointPorosityBase(this)
  implicit none
  class(porosity_base), intent(inout) :: this

  this%porosity_ptr = 0._k_real
  if (this%use_preseeded_porosity) then
    if (((this%gb_id_ptr.eq.1).and.this%use_bulk_porosity_evolution).or.(this%gb_id_ptr.gt.1)) &
    this%porosity_ptr = this%initial_porosity
    this%effective_porosity_ptr = this%initial_porosity 
  endif
  this%void_radius_ptr = 0._k_real
  this%void_density_ptr = 0._k_real
end subroutine

module subroutine updateStateVarsStaggeredPorosityBase(this)
  class(porosity_base), intent(inout) :: this

  this%total_porosity = 0._k_real

  call this%updatePhaseParameters


  associate (ix=>this%ix, xs=>this%xs_rank, xe=>this%xe_rank, &
              iy=>this%iy, ys=>this%ys_rank, ye=>this%ye_rank, &
              iz=>this%iz, zs=>this%zs_rank, ze=>this%ze_rank)

    do iz=zs,ze
      do iy=ys,ye
        do ix=xs,xe
          call this%setPointData(ix,iy,iz)
          if (this%phase_fraction_grid(ix,iy,iz).gt.0._k_real) then
          call this%updateStateVariablesAtMaterialPointOuterLoop()
          endif
        enddo
      enddo
    enddo

  end associate


end subroutine

module subroutine updateStateVarsOuterLoopPorosityBase(this)
  class(porosity_base), intent(inout) :: this
  integer, dimension(3) :: nx_ny_nz

  this%total_porosity = 0._k_real
  nx_ny_nz = shape(this%phase_fraction_grid)

  call this%updatePhaseParameters


  associate (ix=>this%ix, xs=>this%xs_rank, xe=>this%xe_rank, &
            iy=>this%iy, ys=>this%ys_rank, ye=>this%ye_rank, &
            iz=>this%iz, zs=>this%zs_rank, ze=>this%ze_rank)

      do iz=zs,ze
        do iy=ys,ye
          do ix=xs,xe
          call this%setPointData(ix,iy,iz)
          if (this%phase_fraction_grid(ix,iy,iz).gt.0._k_real) then
          call this%updateStateVariablesAtMaterialPointOuterLoop()
          endif
        enddo
      enddo
    enddo

  end associate


end subroutine

module subroutine acceptRejectSolutionPorosityBase(this, dt_max, accept_solution_flag)
  use mpi_useful_routines_mod, only : MPIComputeQualityScalar, MPIMaxIncrementGridSSScalar
  use mpi_variables_mod, only : i_am_mpi_master
  use time_march_mod, only : computeMaxAllowedTimeStepQuality
  use number_to_string_mod
  implicit none
  class(porosity_base), intent(inout) :: this
  real(k_real), intent(out) :: dt_max
  logical, intent(out) :: accept_solution_flag
  real(k_real) :: max_dt_var
  real(k_real) :: Q, max_porosity_increment
  integer :: npoints
  logical :: accept_var

  if (i_am_mpi_master) write(*,*) "*********************************"
  if (i_am_mpi_master) write(*,*) "acceptRejectSolutionPorosityBase"

  accept_solution_flag = .true.
  dt_max = this%dt_max


  npoints = this%grid_data%getGlobalGridNPoints()
  call MPIComputeQualityScalar(this%porosity_grid, this%porosity_grid_old, this%porosity_rate_grid, this%porosity_rate_grid_old, &
                               npoints, this%dt, this%porosity_abs_rel_tol(2), this%porosity_abs_rel_tol(1), Q, max_porosity_increment)
  call computeMaxAllowedTimeStepQuality(Q, this%dt, "porosity ",  accept_var, max_porosity_increment, max_dt_var)
  accept_solution_flag = accept_solution_flag.and.accept_var
  dt_max = min(dt_max, max_dt_var)

  call MPIComputeQualityScalar(this%effective_porosity_grid, this%effective_porosity_grid_old, this%effective_porosity_rate_grid, this%effective_porosity_rate_grid_old, &
                               npoints, this%dt, this%porosity_abs_rel_tol(2), this%porosity_abs_rel_tol(1), Q, max_porosity_increment)
  call computeMaxAllowedTimeStepQuality(Q, this%dt, "effective porosity ",  accept_var, max_porosity_increment, max_dt_var)
  accept_solution_flag = accept_solution_flag.and.accept_var
  dt_max = min(dt_max, max_dt_var)

  if (accept_solution_flag.and.i_am_mpi_master) write(*,*) "solution might be ACCEPTED"
  if (.not.accept_solution_flag.and.i_am_mpi_master) write(*,*) "solution will be REJECTED"
  if (i_am_mpi_master) write(*,*) "new allowable dt is ", dt_max
  if (i_am_mpi_master) write(*,*) "*********************************"
  if (i_am_mpi_master) write(*,*) ""


end subroutine

module subroutine computeVoidVolume(this, r, V)
  use math_constants, only : PI
  implicit none
  class(porosity_base), intent(in) :: this
  real(k_real), intent(in) :: r
  real(k_real), intent(out) :: V
  V = 4._k_real/3._k_real * PI * r**3 * this%h
end subroutine

module subroutine computeVoidVolumeRate(this, r, rdot, Vdot)
  use math_constants, only : PI
  implicit none
  class(porosity_base), intent(in) :: this
  real(k_real), intent(in) :: r, rdot
  real(k_real), intent(out) :: Vdot
  Vdot = 4._k_real * PI * r**2 *rdot * this%h
end subroutine


module subroutine writeAverageQuantitiesPorosityBase(this, csv_writer_obj, write_headers)
  use mpi_useful_routines_mod, only : MPIMaximumGridScalarMasked, MPIAverageGridScalarMasked
  use csv_writer_mod, only : csv_writer
  use number_to_string_mod, only : int2string

  implicit none
  class(porosity_base), intent(inout) :: this
  class(csv_writer), intent(inout) :: csv_writer_obj
  logical, intent(in) :: write_headers
  real(k_real) :: max_local_porosity, max_local_por, &
                  max_porosity_nucleation, max_porosity_growth, &
                  average_porosity, average_por, &
                  average_porosity_nucleation, average_porosity_growth, &
                  average_GB_porosity, average_GB_por
  
  call MPIMaximumGridScalarMasked(this%effective_porosity_grid, this%phase_fraction_grid.gt.0._k_real, max_local_porosity)
  call MPIMaximumGridScalarMasked(this%porosity_grid, this%phase_fraction_grid.gt.0._k_real, max_local_por)
  call MPIMaximumGridScalarMasked(this%porosity_nucleation_grid, this%phase_fraction_grid.gt.0._k_real, max_porosity_nucleation)
  call MPIMaximumGridScalarMasked(this%porosity_growth_grid, this%phase_fraction_grid.gt.0._k_real, max_porosity_growth)
  call MPIAverageGridScalarMasked(this%effective_porosity_grid, this%phase_fraction_grid.gt.0._k_real, average_porosity)
  call MPIAverageGridScalarMasked(this%porosity_grid, this%phase_fraction_grid.gt.0._k_real , average_por)
  call MPIAverageGridScalarMasked(this%porosity_nucleation_grid, this%phase_fraction_grid.gt.0._k_real, average_porosity_nucleation)
  call MPIAverageGridScalarMasked(this%porosity_growth_grid, this%phase_fraction_grid.gt.0._k_real, average_porosity_growth)
  call MPIAverageGridScalarMasked(this%porosity_grid, (this%phase_fraction_grid.gt.0._k_real).and.(this%gb_id_grid.gt.1), average_GB_porosity)
  call MPIAverageGridScalarMasked(this%effective_porosity_grid, (this%phase_fraction_grid.gt.0._k_real).and.(this%gb_id_grid.gt.1), average_GB_por)

  call csv_writer_obj%AppendScalar(max_local_porosity,"max_effective_local_porosity", write_headers)
  call csv_writer_obj%AppendScalar(max_local_por,"max_local_porosity", write_headers)
  call csv_writer_obj%AppendScalar(max_porosity_nucleation,"max_porosity_nucleation", write_headers)
  call csv_writer_obj%AppendScalar(max_porosity_growth,"max_porosity_growth", write_headers)
  call csv_writer_obj%AppendScalar(average_porosity,"effective_average_porosity", write_headers)
  call csv_writer_obj%AppendScalar(average_por,"average_porosity", write_headers)
  call csv_writer_obj%AppendScalar(average_porosity_nucleation,"average_porositynucleation", write_headers)
  call csv_writer_obj%AppendScalar(average_porosity_growth,"average_porosity_growth", write_headers)
  call csv_writer_obj%AppendScalar(average_GB_porosity,"average_GB_porosity", write_headers)
  call csv_writer_obj%AppendScalar(average_GB_por,"effective_average_GB_porosity", write_headers)

end subroutine

end submodule porosity_base_smod
