module embedded_slip_systems_mod
#include "macro_debug.fpp"
use kinds

implicit none
! CRYSTAL TYPE ENUM 
integer, parameter :: number_of_crystal_types = 5
! the grid_data_variable_types
! ALERT!!! if you add a griddata variable type you need also to add it to the enum definition,
! to the all_crystal_type_enum array and to all_crystal_type_enum
public :: crystal_type_enum, &
          FCC, BCC, HCP, CUSTOM, NONE


! enum of grid_data variable types
enum, bind(c)
enumerator :: crystal_type_enum = 0 !-> used as kind for utilizing enumeration
enumerator :: NONE=1, FCC=2, BCC=3, HCP=4, CUSTOM=5
end enum

! a vector parameter including all variable types
integer(kind(crystal_type_enum)), dimension(number_of_crystal_types), parameter :: all_crystal_type_enum= &
              (/NONE, FCC, BCC, HCP, CUSTOM/)

! this array is just for error printying error purposes, do not use it anywhere else
! character(len=6), dimension(number_of_crystal_types), parameter, private ::&
!   all_crystal_type_names=(/&
!   "NONE  ",&
!   "FCC   ",&
!   "BCC   ",&
!   "HCP   ",&
!   "CUSTOM"&
!   /)
! END OF CRYSTAL TYPE ENUM DECLARATION
  
contains

subroutine stringToEnumCrystalType(ct_string, ct_type)
  character(len=*), intent(in) :: ct_string
  integer(kind(crystal_type_enum)), intent(out) :: ct_type
  select case (trim(adjustl(ct_string)))
    case("FCC")
      ct_type = FCC
    case("BCC")
      ct_type = BCC
    case("HCP")
      ct_type = HCP
    case("NONE")
      ct_type = NONE
    case("CUSTOM")
      ct_type = CUSTOM
    case default
      error stop "stringToEnumCrytalType wrong type"
  end select

end subroutine

subroutine computeSchmidTensorFromNormalAndDirection(ss_normal, ss_direction, schmid_tensor)
  use tensor_math_mod, only : vector3OuterProduct, normalizeVector
  real(k_real), intent(in), dimension(3) :: ss_normal, ss_direction
  real(k_real), dimension(3) :: unit_ss_normal, unit_ss_direction
  real(k_real), dimension(3,3), intent(out) :: schmid_tensor

  ! first we normalize the slip normal and direction
  call normalizeVector(ss_normal, unit_ss_normal)
  call normalizeVector(ss_direction, unit_ss_direction)
  schmid_tensor = vector3OuterProduct(unit_ss_direction, unit_ss_normal)

end subroutine

subroutine computeClimbTensorFromNormalDirectionAndPsi(ss_normal, ss_direction, psi_climb, climb_tensor)
  use tensor_math_mod, only : vector3OuterProduct, normalizeVector, vectorCrossProduct, getSymmetricPart
  use change_tensor_basis, only :chg_basis_tensor2_to_vector6

  real(k_real), intent(in), dimension(3) :: ss_normal, ss_direction
  real(k_real), intent(in) :: psi_climb
  real(k_real), dimension(3) :: unit_ss_normal, unit_ss_direction, b_cross_n, chi
  real(k_real), intent(out), dimension(3,3) :: climb_tensor

  ! first we normalize the slip normal and direction
  call normalizeVector(ss_normal, unit_ss_normal)
  call normalizeVector(ss_direction, unit_ss_direction)
  ! we compute the third vector to complete the coordinate system (bxn)
  call vectorCrossProduct(unit_ss_direction, unit_ss_normal, b_cross_n) ! vector perpendicular to slip direction and normal

  ! chi is the climb velocity projected on the dislocation plane which is calcualte as:
  chi = unit_ss_direction*sin(psi_climb) - b_cross_n*cos(psi_climb) !projected vector

  climb_tensor = vector3OuterProduct(unit_ss_direction, chi)

end subroutine

subroutine computeClimbDirectionTensor(ss_direction, climb_tensor)
  use tensor_math_mod, only : vector3OuterProduct, normalizeVector, vectorCrossProduct, getSymmetricPart
  use change_tensor_basis, only :chg_basis_tensor2_to_vector6

  real(k_real), intent(in), dimension(3) :: ss_direction
  real(k_real), dimension(3) :: unit_ss_direction
  real(k_real), intent(out), dimension(3,3) :: climb_tensor

  ! first we normalize the slip normal and direction
  call normalizeVector(ss_direction, unit_ss_direction)

  climb_tensor = vector3OuterProduct(unit_ss_direction, unit_ss_direction)

end subroutine

subroutine getNormalizedSlipSystemDirectionandNormal(ss_type, n_slip_system, ss_direction, ss_normal, ca_ratio)
  use tensor_math_mod, only : normalizeVector
  use FCC_SS_mod, only : getFCC_111_110, getFCC_111_112
  use BCC_SS_mod, only : getBCC_110_111, getBCC_112_111
  use HCP_SS_mod, only : getHCP_0001_1120, getHCP_1100_1120, getHCP_1101_1120, getHCP_1101_1123
  use print_utils_mod, only : printToScreen
  use mpi_variables_mod, only : i_am_mpi_master
  use log_file_mod, only : write_detailed_log_to_screen
  character(len=*), intent(in) :: ss_type
  integer, intent(out) :: n_slip_system
  real(k_real), dimension(:,:), pointer, intent(inout) :: ss_direction, ss_normal
  real(k_real), intent(in), optional :: ca_ratio
  integer :: ss_idx

  if (associated(ss_direction)) error stop "getNormalizedSlipSystemDirectionandNormal: ss_direction already associated. Abort "
  if (associated(ss_normal)) error stop "getNormalizedSlipSystemDirectionandNormal: ss_normal already associated. Abort "

  select case(ss_type)
  case("FCC_{111}<110>")
    call getFCC_111_110(ss_normal, ss_direction, n_slip_system)
  case("FCC_{111}<112>")
    call getFCC_111_112(ss_normal, ss_direction, n_slip_system)
  case("BCC_{110}<111>")
    call getBCC_110_111(ss_normal, ss_direction, n_slip_system)
  case("BCC_{112}<111>")
    call getBCC_112_111(ss_normal, ss_direction, n_slip_system)

  case("HCP_{0001}<1120>")
    if (.not.present(ca_ratio)) error stop "you must provide a c/a ratio for a Hexagonal crystal"
    call getHCP_0001_1120(ss_normal, ss_direction, ca_ratio, n_slip_system)
  case("HCP_{1100}<1120>")
    if (.not.present(ca_ratio)) error stop "you must provide a c/a ratio for a Hexagonal crystal"
    call getHCP_1100_1120(ss_normal, ss_direction, ca_ratio, n_slip_system)
  case("HCP_{1101}<1120>")
    if (.not.present(ca_ratio)) error stop "you must provide a c/a ratio for a Hexagonal crystal"
    call getHCP_1101_1120(ss_normal, ss_direction, ca_ratio, n_slip_system)
  case("HCP_{1101}<1123>")
    if (.not.present(ca_ratio)) error stop "you must provide a c/a ratio for a Hexagonal crystal"
    call getHCP_1101_1123(ss_normal, ss_direction, ca_ratio, n_slip_system)

  case default
    if (i_am_mpi_master.and.write_detailed_log_to_screen) write(*,*) "computeSchmidVector slip system type ", ss_type, " not supported. "
    if (i_am_mpi_master.and.write_detailed_log_to_screen) write(*,*) "available options are: FCC_{111}<110>, FCC_{111}<112>, BCC_{110}<111>, BCC_{112}<111> or "
    if (i_am_mpi_master.and.write_detailed_log_to_screen) write(*,*) "HCP_{0001}<1120>, HCP_{1100}<1120>, HCP_{1101}<1120>, HCP_{1101}<1123>"
    if (i_am_mpi_master.and.write_detailed_log_to_screen) write(*,*) "Abort"
    error stop "Error when computing the the schmidVectors"
  end select

  if (i_am_mpi_master.and.write_detailed_log_to_screen) write(*,*) "n_ss ", n_slip_system
  do ss_idx = 1,n_slip_system
    call normalizeVector(ss_normal(:,ss_idx))
    call normalizeVector(ss_direction(:,ss_idx))
    ! call printToScreen(transpose(ss_normal), "ss_normal")
    ! call printToScreen(transpose(ss_direction), "ss_direction")
  enddo
end subroutine

subroutine getNormalizedLoopsNormal(loop_type, n_loop_system, loop_normal, ca_ratio)
  use tensor_math_mod, only : normalizeVector
  use BCC_SS_mod, only : getBCC_100_loop, getBCC_111_loop
  use print_utils_mod, only : printToScreen
  use mpi_variables_mod, only : i_am_mpi_master
  use log_file_mod, only : write_detailed_log_to_screen
  character(len=*), intent(in) :: loop_type
  integer, intent(out) :: n_loop_system
  real(k_real), dimension(:,:), pointer, intent(inout) :: loop_normal
  real(k_real), intent(in), optional :: ca_ratio
  integer :: loop_idx
  __DECL_UNUSED_REAL__


  if (present(ca_ratio)) then 
    __SUPPRESS_UNUSED_REAL__(ca_ratio)
  endif
  
  if (associated(loop_normal)) error stop "getNormalizedLoopsNormal: loop_normal already associated. Abort "

  select case(loop_type)
  case("BCCLoop_{100}")
    call getBCC_100_loop(loop_normal, n_loop_system)
  case("BCCLoop_{111}")
    call getBCC_111_loop(loop_normal, n_loop_system)

  case default
    if (i_am_mpi_master.and.write_detailed_log_to_screen) write(*,*) "compute normals of loop type ", loop_type, " not supported. "
    if (i_am_mpi_master.and.write_detailed_log_to_screen) write(*,*) "available options are: BCCLoop_{100}, BCCLoop_{111} "
    if (i_am_mpi_master.and.write_detailed_log_to_screen) write(*,*) "Abort"
    error stop "Error when computing the the loop normal"
  end select

  if (i_am_mpi_master.and.write_detailed_log_to_screen) write(*,*) "n_loops ", n_loop_system
  do loop_idx = 1,n_loop_system
    call normalizeVector(loop_normal(:,loop_idx))
  enddo
  call printToScreen(transpose(loop_normal), "loop_normal")
end subroutine

end module
