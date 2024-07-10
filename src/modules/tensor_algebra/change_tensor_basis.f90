module change_tensor_basis
  use kinds
  implicit none

  real(k_real), parameter :: RSQ2=0.70710678118654744_k_real
  real(k_real), parameter :: RSQ3=0.57735026918962584_k_real
  real(k_real), parameter :: RSQ6=0.40824829046386304_k_real
  real(k_real), parameter :: O=0._k_real
  real(k_real), dimension(3,3,6), parameter :: B_change_base = reshape(&
  (/-RSQ2,O,   O,     O,    RSQ2, O,    O,   O,   O, &
    -RSQ6,O,   O,     O,   -RSQ6, O,    O,   O,  2._k_real*RSQ6, &
     O,   O,   O,     O,   O,    RSQ2,  O,   RSQ2,O, &
     O,   O,   RSQ2,  O,   O,    O,     RSQ2,O,   O, &
     O,   RSQ2,O,     RSQ2,O,    O,     O,   O,   O, &
     RSQ3,O,   O,     O,   RSQ3, O,     O,   O,   RSQ3/), (/3,3,6/))

  ! the derivative of the transforamtion from tensor2 to vector6
  real(k_real), parameter, dimension(3,3,6) :: dtensor2_dvector6 = B_change_base

  ! all the constants and base methods should not be accessible from outside these module
  private :: RSQ2, RSQ3, RSQ6, O, B_change_base, &
             chg_basis_vector_to_tensor2_fun, chg_basis_tensor2_to_vector_fun, &
             chg_basis_matrix_to_tensor4, chg_basis_tensor4_to_matrix
             

  ! public methods accessible from outside
  public :: chg_basis_tensor2_to_vector5, chg_basis_tensor2_to_vector6, &
            chg_basis_vector5_to_tensor2, chg_basis_vector6_to_tensor2, &
            chg_basis_tensor4_to_matrix55, chg_basis_tensor4_to_matrix66, &
            chg_basis_matrix55_to_tensor4, chg_basis_matrix66_to_tensor4, &
            dtensor2_dvector6, &
            chg_basis_vector6_to_tensor2_fun, chg_basis_tensor2_to_vector6_fun
contains

! all rank2 change of basis subroutine

  ! case IOPT=2, KDIM=5
  subroutine chg_basis_tensor2_to_vector5(t2, v5)
    integer, parameter :: KDIM=5
    real(k_real), intent(in) :: t2(3,3)
    real(k_real), intent(out) :: v5(KDIM)

    v5 = chg_basis_tensor2_to_vector_fun(t2, KDIM)
  end subroutine

  ! case IOPT=2, KDIM=6
  subroutine chg_basis_tensor2_to_vector6(t2, v6)
    integer, parameter :: KDIM=6
    real(k_real), intent(in) :: t2(3,3)
    real(k_real), intent(out) :: v6(KDIM)

    v6 = chg_basis_tensor2_to_vector_fun(t2, KDIM)
  end subroutine

  function chg_basis_tensor2_to_vector6_fun(t2) result(v6)
    integer, parameter :: KDIM=6
    real(k_real), intent(in) :: t2(3,3)
    real(k_real) :: v6(KDIM)

    v6 = chg_basis_tensor2_to_vector_fun(t2, KDIM)
  end function

  ! case IOPT=1, KDIM=5
  subroutine chg_basis_vector5_to_tensor2(v5, t2)
    integer, parameter :: KDIM=5
    real(k_real), intent(out) :: t2(3,3)
    real(k_real), intent(in) :: v5(KDIM)

    t2 = chg_basis_vector_to_tensor2_fun(v5,KDIM)
  end subroutine

  ! case IOPT=1, KDIM=6
  subroutine chg_basis_vector6_to_tensor2(v6, t2)
    integer, parameter :: KDIM=6
    real(k_real), intent(out) :: t2(3,3)
    real(k_real), intent(in) :: v6(KDIM)

    t2 = chg_basis_vector_to_tensor2_fun(v6,KDIM)
  end subroutine

  function chg_basis_vector6_to_tensor2_fun(v6) result(t2)
    integer, parameter :: KDIM=6
    real(k_real), intent(in) :: v6(KDIM)
    real(k_real) :: t2(3,3)

    t2 = chg_basis_vector_to_tensor2_fun(v6,KDIM)
  end function

! all rank4 change of basis subroutine
  ! case IOPT=4, KDIM=5
  subroutine chg_basis_tensor4_to_matrix55(t4, m55)
    integer, parameter :: KDIM=5
    real(k_real), intent(in) :: t4(3,3,3,3)
    real(k_real), intent(out) :: m55(KDIM,KDIM)

    call chg_basis_tensor4_to_matrix(t4,m55,KDIM)
  end subroutine

  ! case IOPT=4, KDIM=6
  subroutine chg_basis_tensor4_to_matrix66(t4, m66)
    integer, parameter :: KDIM=6
    real(k_real), intent(in) :: t4(3,3,3,3)
    real(k_real), intent(out) :: m66(KDIM,KDIM)

    call chg_basis_tensor4_to_matrix(t4,m66,KDIM)
  end subroutine

  ! case IOPT=3, KDIM=5
  subroutine chg_basis_matrix55_to_tensor4(m55, t4)
    integer, parameter :: KDIM=5
    real(k_real), intent(out) :: t4(3,3,3,3)
    real(k_real), intent(in) :: m55(KDIM,KDIM)

    call chg_basis_matrix_to_tensor4(m55,t4,KDIM)
  end subroutine

  ! case IOPT=3, KDIM=6
  subroutine chg_basis_matrix66_to_tensor4(m66, t4)
    integer, parameter :: KDIM=6
    real(k_real), intent(out) :: t4(3,3,3,3)
    real(k_real), intent(in) :: m66(KDIM,KDIM)

    call chg_basis_matrix_to_tensor4(m66,t4,KDIM)
  end subroutine

  ! I'll leave here the original comment of the change base routine
  ! I split the 4 cases into 4 subroutine, this allows to have more control of
  ! input output, and we can set an intent to the dummy arguments.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !C ************************************************************************
  !C     SUBROUTINE CHG_BASIS    --->   VERSION 19/JUL/01
  !C
  !C     (modif. 06/FEB/98 - same convention as SELFPOLY - C.N.T.)
  !C     (modif. 16/JUN/99 - same convention as Maudlin  - C.N.T.)
  !C     (modif. 10/MAY/01 - KDIM version - R.L.)
  !C
  !C     KDIM=5 or 6, FOR DEVIATORIC or DEV+HYDROST TENSORS, RESPECTIVELY.
  !C     IOPT=0: DEFINES A BASIS OF 6 SECOND ORDER TENSORS B(N).
  !C     IOPT=1: CALCULATES SECOND ORDER TENSOR 'C2' AS AN EXPANSION IN TERMS
  !C             OF VECTOR COMPONENTS CE2(KDIM) AND THE BASIS TENSORS B(KDIM).
  !C     IOPT=2: CALCULATES COMPONENTS OF C2 AS A VECTOR CE2(KDIM).
  !C     IOPT=3: CALCULATES FOURTH ORDER TENSOR 'C4' AS AN EXPANSION IN TERMS
  !C             OF MATRIX COMPONENTS CE4(K,K) AND THE BASIS TENSORS B(KDIM).
  !C     IOPT=4: CALCULATES MATRIX COMPONENTS CE4(K,K) OF TENSOR 'C4'.
  !C **************************************************************************

function chg_basis_vector_to_tensor2_fun(vector, KDIM) result(tensor2)
  integer, intent(in) :: KDIM
  real(k_real), intent(in) :: vector(KDIM)
  real(k_real) :: tensor2(3,3)
  integer :: i,j,n

  ! IOPT = 1
  !C *** CALCULATES CARTESIAN SECOND ORDER TENSOR FROM b-COMPONENTS VECTOR.
  do I=1,3
    do J=1,3
      tensor2(I,J)=0.0
        do N=1,KDIM
          tensor2(I,J)=tensor2(I,J)+vector(N)*B_change_base(I,J,N)
        enddo
      enddo
  enddo
end function

function chg_basis_tensor2_to_vector_fun(tensor2, KDIM) result(vector)
  integer, intent(in) :: KDIM
  real(k_real), intent(in) :: tensor2(3,3)
  real(k_real) :: vector(KDIM)
  integer :: i,j,n

  ! IOPT = 2
  !C *** CALCULATES KDIMx1 b-COMPONENTS VECTOR FROM SECOND ORDER TENSOR.
  do N=1,KDIM
    vector(N)=0.0
    do I=1,3
      do J=1,3
        vector(N)=vector(N)+tensor2(I,J)*B_change_base(I,J,N)
      enddo
    enddo
  enddo
end function

subroutine chg_basis_matrix_to_tensor4(matrix, tensor4, KDIM)
  integer, intent(in) :: KDIM
  real(k_real), intent(out) :: tensor4(3,3,3,3)
  real(k_real), intent(in) :: matrix(KDIM, KDIM)
  integer :: i,j,k,l,m,n

  ! IOPT = 3
  !C *** CALCULATES FOURTH ORDER TENSOR FROM b-COMPONENTS MATRIX.
  do I=1,3
    do J=1,3
      do K=1,3
        do L=1,3
          tensor4(I,J,K,L)=0.0
          do N=1,KDIM
            do M=1,KDIM
              tensor4(I,J,K,L)=tensor4(I,J,K,L)+matrix(N,M)*B_change_base(I,J,N)*B_change_base(K,L,M)
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo

end subroutine

subroutine chg_basis_tensor4_to_matrix(tensor4, matrix, KDIM)
  integer, intent(in) :: KDIM
  real(k_real), intent(in) :: tensor4(3,3,3,3)
  real(k_real), intent(out) :: matrix(KDIM, KDIM)
  integer :: i,j,k,l,m,n

  ! IOPT = 4
  !C *** CALCULATES KDIMxKDIM b-COMPONENTS MATRIX FROM FOURTH ORDER TENSOR.
  do N=1,KDIM
    do M=1,KDIM
      matrix(N,M)=0.0
      do I=1,3
        do J=1,3
          do K=1,3
            do L=1,3
              matrix(N,M)=matrix(N,M)+tensor4(I,J,K,L)*B_change_base(I,J,N)*B_change_base(K,L,M)
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo

end subroutine

end module change_tensor_basis
