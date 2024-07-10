module HCP_SS_mod
  use kinds
  implicit none
  ! BASAL
  real(k_real), parameter, dimension(8,3) :: HCP_0001_1120=reshape((/ &
  0, 0, 0, 1,   -2,  1,  1, 0,&
  0, 0, 0, 1,    1, -2,  1, 0,&
  0, 0, 0, 1,    1,  1, -2, 0/), (/8,3/))

  ! PRISMATIC
  real(k_real), parameter, dimension(8,3) :: HCP_1100_1120=reshape((/ &
  0,  1, -1, 0,   -2,  1,  1, 0,&
 -1,  0,  1, 0,    1, -2,  1, 0,&
  1, -1,  0, 0,    1,  1, -2, 0/), (/8,3/))

  ! PYRAMIDAL<a>
  real(k_real), parameter, dimension(8,6) :: HCP_1101_1120=reshape((/ &
  0,  1, -1, 1,    -2,  1,  1, 0,&
 -1,  0,  1, 1,     1, -2,  1, 0,&
  1, -1,  0, 1,     1,  1, -2, 0,&
  0, -1,  1, 1,     2, -1, -1, 0,&
  1,  0, -1, 1,    -1,  2, -1, 0,&
 -1,  1,  0, 1,    -1, -1,  2, 0/), (/8,6/))

  ! PYRAMIDAL<c+a>first
  real(k_real), parameter, dimension(8,12) :: HCP_1101_1123=reshape((/ &
  1,  0, -1, 1,    -2,  1,  1, 3,&
 -1,  1,  0, 1,     1, -2,  1, 3,&
  0, -1,  1, 1,     1,  1, -2, 3,&
 -1,  0,  1, 1,     2, -1, -1, 3,&
  1, -1,  0, 1,    -1,  2, -1, 3,&
  0,  1, -1, 1,    -1, -1,  2, 3,&
  1, -1,  0, 1,    -2,  1,  1, 3,&
  0,  1, -1, 1,     1, -2,  1, 3,&
 -1,  0,  1, 1,     1,  1, -2, 3,&
 -1,  1,  0, 1,     2, -1, -1, 3,&
  0, -1,  1, 1,    -1,  2, -1, 3,&
  1, -1,  0, 1,    -2,  1,  1, 3/), (/8,12/))

  ! CLIMB DIRECTIONS A
  real(k_real), parameter, dimension(4,4) :: HCP_CLIMB_DIRECTIONS=reshape((/ &
  -2,  1,  1, 0,&
   1, -2,  1, 0,&
   1,  1, -2, 0,&
   0,  0,  0, 1 /), (/4,4/))


! HCP SLIP SYSTEM TYPE ENUM 
integer, parameter :: number_of_hcp_ss_types = 4
! the grid_data_variable_types
! ALERT!!! if you add a griddata variable type you need also to add it to the enum definition,
! to the all_hcp_ss_types_enum array and to all_hcp_ss_type_enum
public :: hcp_ss_type_enum, &
          HCP_0001_1120_ss, HCP_1100_1120_ss, HCP_1101_1120_ss, HCP_1101_1123_ss


! enum of grid_data variable types
enum, bind(c)
enumerator :: hcp_ss_type_enum = 0 !-> used as kind for utilizing enumeration
enumerator :: HCP_0001_1120_ss=1, HCP_1100_1120_ss=2, HCP_1101_1120_ss=3, HCP_1101_1123_ss=4
end enum

! a vector parameter including all variable types
integer(kind(hcp_ss_type_enum)), dimension(number_of_hcp_ss_types), parameter :: all_hcp_ss_type_enum= &
              (/HCP_0001_1120_ss, HCP_1100_1120_ss, HCP_1101_1120_ss, HCP_1101_1123_ss/)

! this array is just for error printying error purposes, do not use it anywhere else
! character(len=16), dimension(number_of_hcp_ss_types), parameter, private ::&
!   all_hcp_ss_type_names=(/&
!   "HCP_0001_1120_ss",&
!   "HCP_1100_1120_ss",&
!   "HCP_1101_1120_ss",&
!   "HCP_1101_1123_ss"&
!   /)
! END OF HCP SLIP SYSTEM TYPE ENUM DECLARATION

contains

subroutine StringToEnumHCPSS(hcp_ss_string, hcp_ss_type)
  character(len=*), intent(in) :: hcp_ss_string
  integer(kind(hcp_ss_type_enum)), intent(out) :: hcp_ss_type
  select case (trim(adjustl(hcp_ss_string)))
    case("HCP_{0001}<1120>")
      hcp_ss_type = HCP_0001_1120_ss
    case("HCP_{1100}<1120>")
      hcp_ss_type = HCP_1100_1120_ss
    case("HCP_{1101}<1120>")
      hcp_ss_type = HCP_1101_1120_ss
    case("HCP_{1101}<1123>")
      hcp_ss_type = HCP_1101_1123_ss
    case default
      error stop "stringToEnumCrytalType wrong type"
  end select
end subroutine


  subroutine getHCP_0001_1120(slip_normal, slip_direction, ca_ratio, n_ss_out)
    real(k_real), intent(in) :: ca_ratio
    real(k_real), dimension(:,:), intent(inout), pointer :: slip_normal, slip_direction
    integer, intent(out) :: n_ss_out
    integer, parameter :: n_ss = 3

    if (associated(slip_normal)) error stop "getHCP_0001_1120: slip_normal already associated. Abort "
    if (associated(slip_direction)) error stop "getHCP_0001_1120: slip_direction already associated. Abort "
    allocate(slip_normal(3,n_ss), slip_direction(3,n_ss))

    call Hex2CartesianSlipSystems(HCP_0001_1120(1:4,:), HCP_0001_1120(5:8,:), ca_ratio, slip_normal, slip_direction )

    n_ss_out = n_ss
  end subroutine

  subroutine getHCP_1100_1120(slip_normal, slip_direction, ca_ratio, n_ss_out)
    real(k_real), intent(in) :: ca_ratio
    real(k_real), dimension(:,:), intent(inout), pointer :: slip_normal, slip_direction
    integer, intent(out) :: n_ss_out
    integer, parameter :: n_ss = 3

    if (associated(slip_normal)) error stop "getHCP_1100_1120: slip_normal already associated. Abort "
    if (associated(slip_direction)) error stop "getHCP_1100_1120: slip_direction already associated. Abort "
    allocate(slip_normal(3,n_ss), slip_direction(3,n_ss))

    call Hex2CartesianSlipSystems(HCP_1100_1120(1:4,:), HCP_1100_1120(5:8,:), ca_ratio, slip_normal, slip_direction )

    n_ss_out = n_ss
  end subroutine

  subroutine getHCP_1101_1120(slip_normal, slip_direction, ca_ratio, n_ss_out)
    real(k_real), intent(in) :: ca_ratio
    real(k_real), dimension(:,:), intent(inout), pointer :: slip_normal, slip_direction
    integer, intent(out) :: n_ss_out
    integer, parameter :: n_ss = 6

    if (associated(slip_normal)) error stop "getHCP_1101_1120: slip_normal already associated. Abort "
    if (associated(slip_direction)) error stop "getHCP_1101_1120: slip_direction already associated. Abort "
    allocate(slip_normal(3,n_ss), slip_direction(3,n_ss))

    call Hex2CartesianSlipSystems(HCP_1101_1120(1:4,:), HCP_1101_1120(5:8,:), ca_ratio, slip_normal, slip_direction )

    n_ss_out = n_ss
  end subroutine

  subroutine getHCP_1101_1123(slip_normal, slip_direction, ca_ratio, n_ss_out)
    real(k_real), intent(in) :: ca_ratio
    real(k_real), dimension(:,:), intent(inout), pointer :: slip_normal, slip_direction
    integer, intent(out) :: n_ss_out
    integer, parameter :: n_ss = 12

    if (associated(slip_normal)) error stop "getHCP_1101_1123: slip_normal already associated. Abort "
    if (associated(slip_direction)) error stop "getHCP_1101_1123: slip_direction already associated. Abort "
    allocate(slip_normal(3,n_ss), slip_direction(3,n_ss))

    call Hex2CartesianSlipSystems(HCP_1101_1123(1:4,:), HCP_1101_1123(5:8,:), ca_ratio, slip_normal, slip_direction )

    n_ss_out = n_ss
  end subroutine

  subroutine getHCP_climb_directions(climb_direction, ca_ratio, n_ss_out)
    real(k_real), intent(in) :: ca_ratio
    real(k_real), dimension(:,:), intent(inout), pointer :: climb_direction
    integer, intent(out) :: n_ss_out
    integer, parameter :: n_ss = 4

    if (associated(climb_direction)) error stop "getHCP_climb_directions: climb_direction already associated. Abort "
    allocate(climb_direction(3,n_ss))

    call Hex2CartesianDirection(HCP_CLIMB_DIRECTIONS, ca_ratio, climb_direction )

    n_ss_out = n_ss
  end subroutine

  subroutine Hex2CartesianSlipSystems(ss_normal_hex, ss_direction_hex, ca_ratio, ss_normal, ss_direction)
    use math_constants, only : PI
    use units_conversion_mod, only : deg2rad
    use print_utils_mod, only : printToScreen
    real(k_real), intent(in), dimension(:,:) :: ss_normal_hex, ss_direction_hex
    real(k_real), intent(in) :: ca_ratio
    real(k_real), intent(out), dimension(:,:) :: ss_normal, ss_direction
    real(k_real) :: CDIM(3), Q(3,3), vector_3(3)
    real(k_real), dimension(3), parameter :: CANG=(/90._k_real*deg2rad, &
                                                    90._k_real*deg2rad, &
                                                    120._k_real*deg2rad/)
    !
    real(k_real) :: dot_val
    integer :: ss_dims(2), n_ss, ss

    CDIM = (/1._k_real, 1._k_real, ca_ratio/)

    call getProjectionMatrix(CANG, CDIM, Q)

    ss_dims = shape(ss_normal_hex)
    n_ss = ss_dims(2)

    do ss = 1, n_ss
      vector_3(1) =  ss_normal_hex(1, ss)
      vector_3(2) =  ss_normal_hex(2, ss)
      vector_3(3) =  ss_normal_hex(4, ss)
      ss_normal(3,ss)= vector_3(3)/Q(3,3)
      ss_normal(1,ss) = (vector_3(1)-Q(3,1)*ss_normal(3,ss))/Q(1,1)
      ss_normal(2,ss)=(vector_3(2)-Q(1,2)*ss_normal(1,ss)-Q(3,2)*ss_normal(3,ss))/Q(2,2)

      vector_3(1) =  ss_direction_hex(1, ss) - ss_direction_hex(3, ss)
      vector_3(2) =  ss_direction_hex(2, ss) - ss_direction_hex(3, ss)
      vector_3(3) =  ss_direction_hex(4, ss)

      ss_direction(:,ss) = matmul(Q, vector_3)

      dot_val = sum(ss_normal(:,ss)*ss_direction(:,ss))
      if (abs(dot_val).ge.1e-6) then
        write(*,*) " ss_normal_hex ", ss_normal_hex(:,ss)
        write(*,*) " ss_direction_hex ", ss_direction_hex(:,ss)
        call printToScreen(Q, "Q")
        write(*,*) " ss_normal ", ss_normal(:,ss)
        write(*,*) " ss_direction ", ss_direction(:,ss)
        write(*,*) " dot product: ", dot_val
        error stop "slip normal and lisp direction are not orthogonal"
      endif
    enddo
    
  end subroutine

  subroutine Hex2CartesianDirection(ss_direction_hex, ca_ratio, ss_direction)
    use math_constants, only : PI
    use units_conversion_mod, only : deg2rad
    use print_utils_mod, only : printToScreen
    real(k_real), intent(in), dimension(:,:) :: ss_direction_hex
    real(k_real), intent(in) :: ca_ratio
    real(k_real), intent(out), dimension(:,:) :: ss_direction
    real(k_real) :: CDIM(3), Q(3,3), vector_3(3)
    real(k_real), dimension(3), parameter :: CANG=(/90._k_real*deg2rad, &
                                                    90._k_real*deg2rad, &
                                                    120._k_real*deg2rad/)
    !
    integer :: ss_dims(2), n_ss, ss

    CDIM = (/1._k_real, 1._k_real, ca_ratio/)

    call getProjectionMatrix(CANG, CDIM, Q)

    ss_dims = shape(ss_direction_hex)
    n_ss = ss_dims(2)

    do ss = 1, n_ss

      vector_3(1) =  ss_direction_hex(1, ss) - ss_direction_hex(3, ss)
      vector_3(2) =  ss_direction_hex(2, ss) - ss_direction_hex(3, ss)
      vector_3(3) =  ss_direction_hex(4, ss)

      ss_direction(:,ss) = matmul(Q, vector_3)

    enddo

  end subroutine


  subroutine getProjectionMatrix(CANG, CDIM, Q)
    real(k_real), dimension(3), intent(in) :: CANG, CDIM
    real(k_real), dimension(3,3), intent(out) :: Q
    integer :: i,j

    Q(1,1)=SIN(CANG(2))
    Q(2,1)=0.
    Q(3,1)=COS(CANG(2))
    Q(1,2)=(COS(CANG(3))-COS(CANG(1))*COS(CANG(2)))/SIN(CANG(2))
    Q(3,2)=COS(CANG(1))
    Q(2,2)=SQRT(1.-Q(1,2)**2-Q(3,2)**2)
    Q(1,3)=0.
    Q(2,3)=0.
    Q(3,3)=1.

    DO J=1,3
    DO I=1,3
      Q(I,J)=CDIM(J)*Q(I,J)
      if (abs(Q(I,J)).le.1e-6_k_real) Q(I,J)=0._k_real
    ENDDO
    ENDDO

  end subroutine
end module
