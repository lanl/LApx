module common_material_parameter_mod
use kinds
use string_module, only : string_type, string_array
use embedded_slip_systems_mod, only : crystal_type_enum, NONE

implicit none


type :: common_material_parameter
  real(k_real), pointer :: grain_diameter  => null()
  real(k_real), pointer :: atomic_volume  => null()
  real(k_real), pointer :: gb_thickness  => null()
  real(k_real), pointer :: cellwall_thickness  => null()
  real(k_real), pointer :: mass_density  => null()
  real(k_real), pointer, dimension(:) :: vacancy_migration_energy => null(), &
                                         vacancy_formation_energy => null(), &
                                         vacancy_diffusivity_coefficient => null(), &
                                         pipe_diffusivity_coefficient => null(), &
                                         pipe_diffusion_migration_energy => null()
  contains
    procedure :: readCommonMaterialParametersFromFile

end type


contains

subroutine readCommonMaterialParametersFromFile(this, matf_reader)
  use read_from_file_utils
  use units_conversion_mod, only : eV2Joule
  use string_module, only : string_array
  implicit none
  class(common_material_parameter), intent(inout) :: this
  class(file_reader), intent(inout) :: matf_reader
  type(string_array) :: dummy_string_array
  real(k_real) :: temp

  allocate(this%atomic_volume, &
           this%gb_thickness, &
           this%grain_diameter, &
           this%cellwall_thickness, &
           this%mass_density)

  ! this%crystal_type="   "
  call matf_reader%readLineAndCheckStringsAreEqual("--Common-Material-Parameters", dummy_string_array)
  call matf_reader%readParameter("atomic-volume[m^3]", temp); this%atomic_volume = temp
  call matf_reader%readParameter("gb-thickness[m]", temp); this%gb_thickness = temp
  call matf_reader%readParameter("grain-diameter[m]", temp); this%grain_diameter = temp
  call matf_reader%readParameter("cellwall-thickness[m]", temp); this%cellwall_thickness = temp
  call matf_reader%readParameter("mass-density[kg/m^3]", temp); this%mass_density = temp
  call matf_reader%readVectorParameter("vacancy-diffusivity-coefficient[bulk,GB-[m^2/s]]", 2, this%vacancy_diffusivity_coefficient)
  call matf_reader%readVectorParameter("vacancy-formation-energy[bulk,GB-[eV]]", 2, this%vacancy_formation_energy)
  this%vacancy_formation_energy = this%vacancy_formation_energy*eV2Joule
  call matf_reader%readVectorParameter("vacancy-migration-energy[bulk,GB-[eV]]", 2, this%vacancy_migration_energy)
  this%vacancy_migration_energy = this%vacancy_migration_energy*eV2Joule
  call matf_reader%readVectorParameter("pipe-diffusivity-coefficient[bulk,GB-[m^2/s]]", 2, this%pipe_diffusivity_coefficient)
  call matf_reader%readVectorParameter("pipe-diffusion-activation-energy[bulk,GB-[eV]]", 2, this%pipe_diffusion_migration_energy)
  this%pipe_diffusion_migration_energy = this%pipe_diffusion_migration_energy*eV2Joule
end subroutine

end module common_material_parameter_mod
