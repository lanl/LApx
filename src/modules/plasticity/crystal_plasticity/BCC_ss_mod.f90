module BCC_SS_mod
  use kinds
  implicit none
  real(k_real), parameter, dimension(6,12) :: BCC_110_111=reshape((/ &
  0, 1,-1,   1, 1, 1,   &
  1, 0,-1,   1, 1, 1,   &
  1,-1, 0,   1, 1, 1,   &
  0, 1,-1,  -1, 1, 1,   &
  1, 0, 1,  -1, 1, 1,   &
  1, 1, 0,  -1, 1, 1,   &
  0, 1, 1,  -1,-1, 1,   &
  1, 0, 1,  -1,-1, 1,   &
  1,-1, 0,  -1,-1, 1,   &
  0, 1, 1,   1,-1, 1,   &
  1, 0,-1,   1,-1, 1,   &
  1, 1, 0,   1,-1, 1/), (/6,12/))

  real(k_real), parameter, dimension(6,12) :: BCC_112_111=reshape((/ &
 -2, 1, 1,    1, 1, 1,&
  1,-2, 1,    1, 1, 1,&
  1, 1,-2,    1, 1, 1,&
  2,-1, 1,   -1,-1, 1,&
 -1, 2, 1,   -1,-1, 1,&
 -1,-1,-2,   -1,-1, 1,&
 -1,-1,-2,    1,-1,-1,&
  1, 2,-1,    1,-1,-1,&
  1,-1, 2,    1,-1,-1,&
  2, 1,-1,   -1, 1,-1,&
 -1,-2, 1,   -1, 1,-1,&
 -1, 1, 2,   -1, 1,-1/), (/6,12/))

 real(k_real), parameter, dimension(3,4) :: BCC_111_loop=reshape((/ &
  1, 1, 1,&
 -1,-1, 1,&
  1,-1,-1,&
 -1, 1,-1/), (/3,4/))

 real(k_real), parameter, dimension(3,3) :: BCC_100_loop=reshape((/ &
  1, 0, 0,&
  0, 1, 0,&
  0, 0, 1/), (/3,3/))

contains
  subroutine getBCC_110_111(slip_normal, slip_direction, n_ss_out)
    real(k_real), dimension(:,:), intent(inout), pointer :: slip_normal, slip_direction
    integer, intent(out) :: n_ss_out
    integer, parameter :: n_ss = 12

    if (associated(slip_normal)) error stop "getBCC_110_111: slip_normal already associated. Abort "
    if (associated(slip_direction)) error stop "getBCC_110_111: slip_direction already associated. Abort "
    allocate(slip_normal(3,n_ss), slip_direction(3,n_ss))

    slip_normal = BCC_110_111(1:3,:)
    slip_direction = BCC_110_111(4:6,:)
    n_ss_out = n_ss
  end subroutine

  subroutine getBCC_112_111(slip_normal, slip_direction, n_ss_out)
    real(k_real), dimension(:,:), intent(inout), pointer :: slip_normal, slip_direction
    integer, intent(out) :: n_ss_out
    integer, parameter :: n_ss = 12

    if (associated(slip_normal)) error stop "getBCC_112_111: slip_normal already associated. Abort "
    if (associated(slip_direction)) error stop "getBCC_112_111: slip_direction already associated. Abort "
    allocate(slip_normal(3,n_ss), slip_direction(3,n_ss))

    slip_normal = BCC_112_111(1:3,:)
    slip_direction = BCC_112_111(4:6,:)
    n_ss_out = n_ss
  end subroutine

  subroutine getBCC_100_loop(loop_normal, n_loop_out)
    real(k_real), dimension(:,:), intent(inout), pointer :: loop_normal
    integer, intent(out) :: n_loop_out
    integer, parameter :: n_loop = 3

    if (associated(loop_normal)) error stop "getBCC_100_loop: loop_normal already associated. Abort "
    allocate(loop_normal(3,n_loop))

    loop_normal = BCC_100_loop(1:3,:)
    n_loop_out = n_loop
  end subroutine

  subroutine getBCC_111_loop(loop_normal, n_loop_out)
    real(k_real), dimension(:,:), intent(inout), pointer :: loop_normal
    integer, intent(out) :: n_loop_out
    integer, parameter :: n_loop = 4

    if (associated(loop_normal)) error stop "getBCC_111_loop: loop_normal already associated. Abort "
    allocate(loop_normal(3,n_loop))

    loop_normal = BCC_111_loop(1:3,:)
    n_loop_out = n_loop
  end subroutine
end
