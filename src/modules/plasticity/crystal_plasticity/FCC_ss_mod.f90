module FCC_SS_mod
  use kinds
  implicit none
  real(k_real), parameter, dimension(6,12) :: FCC_111_110=reshape((/ &
  1, 1, 1,   0, 1,-1,&
  1, 1, 1,   1, 0,-1,&
  1, 1, 1,   1,-1, 0,&
 -1, 1, 1,   0, 1,-1,&
 -1, 1, 1,   1, 0, 1,&
 -1, 1, 1,   1, 1, 0,&
 -1,-1, 1,   0, 1, 1,&
 -1,-1, 1,   1, 0, 1,&
 -1,-1, 1,   1,-1, 0,&
  1,-1, 1,   0, 1, 1,&
  1,-1, 1,   1, 0,-1,&
  1,-1, 1,   1, 1, 0/), (/6,12/))

  real(k_real), parameter, dimension(6,12) :: FCC_111_112=reshape((/ &
  1, 1, 1,  -2, 1, 1,&
  1, 1, 1,   1,-2, 1,&
  1, 1, 1,   1, 1,-2,&
 -1,-1, 1,   2,-1, 1,&
 -1,-1, 1,  -1, 2, 1,&
 -1,-1, 1,  -1,-1,-2,&
  1,-1,-1,  -1,-1,-2,&
  1,-1,-1,   1, 2,-1,&
  1,-1,-1,   1,-1, 2,&
 -1, 1,-1,   2, 1,-1,&
 -1, 1,-1,  -1,-2, 1,&
 -1, 1,-1,  -1, 1, 2/), (/6,12/))

contains
  subroutine getFCC_111_110(slip_normal, slip_direction, n_ss_out)
    real(k_real), dimension(:,:), intent(inout), pointer :: slip_normal, slip_direction
    integer, intent(out) :: n_ss_out
    integer, parameter :: n_ss = 12

    if (associated(slip_normal)) error stop "getFCC_111_110: slip_normal already associated. Abort "
    if (associated(slip_direction)) error stop "getFCC_111_110: slip_direction already associated. Abort "
    allocate(slip_normal(3,n_ss), slip_direction(3,n_ss))

    slip_normal = FCC_111_110(1:3,:)
    slip_direction = FCC_111_110(4:6,:)
    n_ss_out = n_ss
  end subroutine

  subroutine getFCC_111_112(slip_normal, slip_direction, n_ss_out)
    real(k_real), dimension(:,:), intent(inout), pointer :: slip_normal, slip_direction
    integer, intent(out) :: n_ss_out
    integer, parameter :: n_ss = 12

    if (associated(slip_normal)) error stop "getFCC_111_112: slip_normal already associated. Abort "
    if (associated(slip_direction)) error stop "getFCC_111_112: slip_direction already associated. Abort "
    allocate(slip_normal(3,n_ss), slip_direction(3,n_ss))

    slip_normal = FCC_111_112(1:3,:)
    slip_direction = FCC_111_112(4:6,:)
    n_ss_out = n_ss
  end subroutine
end
