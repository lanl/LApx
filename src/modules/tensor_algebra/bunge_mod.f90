module bunge_mod
use kinds
implicit none

contains

  subroutine bungeRotationMatrixSample2Crystal(ph, th, tm, r_mtx)
    implicit none
    real(k_real), intent(in) :: ph, th, tm
    real(k_real), intent(out) :: r_mtx(3,3)

    associate(cph =>cos(ph), sph => sin(ph), &
              cth =>cos(th), sth => sin(th), &
              ctm =>cos(tm), stm => sin(tm))

      r_mtx(1,1)=ctm*cph-sph*stm*cth
      r_mtx(2,1)=-stm*cph-sph*ctm*cth
      r_mtx(3,1)=sph*sth
      r_mtx(1,2)=ctm*sph+cph*stm*cth
      r_mtx(2,2)=-sph*stm+cph*ctm*cth
      r_mtx(3,2)=-sth*cph
      r_mtx(1,3)=sth*stm
      r_mtx(2,3)=ctm*sth
      r_mtx(3,3)=cth
    end associate
  end subroutine

  subroutine bungeRotationMatrixCrystal2Sample(ph, th, tm, r_mtx)
    implicit none
    real(k_real), intent(in) :: ph, th, tm
    real(k_real), intent(out) :: r_mtx(3,3)

    call bungeRotationMatrixSample2Crystal(ph, th, tm, r_mtx)
    r_mtx =  transpose(r_mtx)
  end subroutine

  subroutine bungeAnlgesFromRotationMatrixSample2Crystal(r_mtx, ph, th, tm)
    use math_constants, only : PI
    use matrix_determinant_mod, only : matrixDeterminat
    use print_utils_mod, only : printToScreen
    implicit none
    real(k_real), intent(in) :: r_mtx(3,3)
    real(k_real), intent(out) :: ph, th, tm
    real(k_real), parameter :: tol = 1e-6
    real(k_real) :: sth
#ifdef __DEBUG__
    real(k_real) :: detR

  if (abs(r_mtx(3,3)).gt.1.001_k_real) error stop "bungeAnlgesFromRotationMatrixSample2Crystal r_mtx(3,3) > 1.001. Abort!"
  call matrixDeterminat(r_mtx, detR)
  if (abs(detR).gt.1.001_k_real) error stop "bungeAnlgesFromRotationMatrixSample2Crystal Your matrix has a determinant larger than 1. Abort"
#endif

    th =  acos(max(min(r_mtx(3,3),1._k_real),-1._k_real))

    if (abs(r_mtx(3,3)) > 1._k_real-tol ) then
       tm=0.
       ph=atan2(r_mtx(1,2),r_mtx(1,1))
    else
       sth=sin(th)
       tm=atan2(r_mtx(1,3)/sth,r_mtx(2,3)/sth)
       ph=atan2(r_mtx(3,1)/sth,-r_mtx(3,2)/sth)
    endif
  end subroutine

  subroutine bungeAnlgesFromRotationMatrixCrystal2Sample(r_mtx, ph, th, tm)
    implicit none
    real(k_real), intent(in) :: r_mtx(3,3)
    real(k_real), intent(out) :: ph, th, tm

    call bungeAnlgesFromRotationMatrixSample2Crystal(transpose(r_mtx), ph, th, tm)
  end subroutine

  subroutine rotMatrix2AxisAnglePair(r_mtx, axis, angle)
    use tensor_math_mod, only : computeTrace, vectorNorm
    use math_constants, only : PI
    implicit none
    real(k_real), intent(in) :: r_mtx(3,3)
    real(k_real), intent(out) :: axis(3), angle
    real(k_real) :: r_trace
    real(k_real), parameter :: tol = 1e-6
    integer :: max_diag_idx, i
    r_trace = computeTrace(r_mtx)

    ! the axis angle rotation matrix is defined as
    ! | cth+v1^2(1-cth)    v1v2(1-cth)-v3sth v1v3(1-cth)+v2sth |
    ! | v2v1(1-cth)-v3sth  cth+v2^2(1-cth)   v2v3(1-cth)-v1sth |
    ! | v3v1(1-cth)-v2sth  v3v2(1-cth)+v1sth cth+v3^2(1-cth)   |

    ! let's first take care of special cases
    if (abs(r_trace-3._k_real).le.tol) then !-> case theta = 0 i.e. R=I33
      ! | cth+v1^2(1-cth)    v1v2(1-cth)-v3sth v1v3(1-cth)+v2sth |
      ! | v1v2(1-cth)-v3sth  cth+v2^2(1-cth)   v2v3(1-cth)-v1sth |
      ! | v1v3(1-cth)-v2sth  v3v2(1-cth)+v1sth cth+v3^2(1-cth)   |
      angle = 0._k_real
      axis = (/0._k_real,0._k_real,1._k_real/)
    else if (abs(r_trace+1._k_real).le.tol) then !-> case theta = +/- pi
      angle = PI
      ! we need to determine the rotation axis
      ! for this case the rotation matrices to
      ! | 2v1^2-1  2v1v2    2v1v3   |
      ! | 2v2v1    2v2^2-1  2v2v3   |
      ! | 2c3c1    2v1v2    2v3^2-1 |
      ! and hence we can determine the rotation axis using the largest diagonal value,
      ! adding to it 1 and normalizing the assocaited row

      !determine the largest diagonal element
      max_diag_idx = 1
      do i =2,3
        if (r_mtx(i,i) > r_mtx(i-1,i-1)) max_diag_idx = i
      enddo

      axis = r_mtx(max_diag_idx, :)
      axis(max_diag_idx) = axis(max_diag_idx) + 1._k_real
      ! now normalize axis
      axis = axis/vectorNorm(axis)
    else ! this is the general case
      angle = acos((r_trace-1._k_real) * 0.5_k_real)
      axis(1) = r_mtx(3,2) - r_mtx(2,3)
      axis(2) = r_mtx(1,3) - r_mtx(3,1)
      axis(3) = r_mtx(2,1) - r_mtx(1,2)
      axis = axis*1._k_real/(2._k_real*sin(angle))
    end if
  end subroutine

  subroutine RodriguesRotationFromAxisAngle(axis, angle, R)
    use tensor_math_mod, only : I33
    implicit none
    real(k_real), intent(in) :: axis(3), angle
    real(k_real), intent(out) :: R(3,3)
    real(k_real) :: K(3,3)
    integer :: i
    ! construct the matrix K from v
    K(1,2) = -axis(3)
    K(2,1) = -K(1,2);

    K(1,3) = axis(2)
    K(3,1) = -K(1,3);

    K(2,3) = -axis(1);
    K(3,2) = -K(2,3);

    do i=1,3
      K(i,i) = 0
    enddo

    ! compute the rotation matrix using rodrigues formula
    R = I33 +sin(angle)*K+(1-cos(angle))*matmul(K,transpose(K))
  end subroutine

  subroutine RodriguesRotationFromDeltaRotationMatrix(deltaR, R)
    implicit none
    real(k_real), intent(in) :: deltaR(3,3)
    real(k_real), intent(out) :: R(3,3)
    real(k_real) :: axis(3), angle

    call rotMatrix2AxisAnglePair(deltaR, axis, angle)
    call RodriguesRotationFromAxisAngle(axis, angle, R)
  end subroutine
end module
