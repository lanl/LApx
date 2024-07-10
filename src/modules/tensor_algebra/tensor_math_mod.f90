module tensor_math_mod
  use kinds
  implicit none

#include "macro_debug.fpp"

  real(k_real), parameter, dimension(6,6) :: I66=reshape((/&
    1,0,0,0,0,0,&
    0,1,0,0,0,0,&
    0,0,1,0,0,0,&
    0,0,0,1,0,0,&
    0,0,0,0,1,0,&
    0,0,0,0,0,1/), (/6,6/))

  real(k_real), parameter, dimension(6,6) :: I55_66=reshape((/&
    1,0,0,0,0,0,&
    0,1,0,0,0,0,&
    0,0,1,0,0,0,&
    0,0,0,1,0,0,&
    0,0,0,0,1,0,&
    0,0,0,0,0,0/), (/6,6/))

  real(k_real), parameter, dimension(5,5) :: I55=reshape((/&
    1,0,0,0,0,&
    0,1,0,0,0,&
    0,0,1,0,0,&
    0,0,0,1,0,&
    0,0,0,0,1/), (/5,5/))

  real(k_real), parameter, dimension(3,3) :: I33=reshape((/&
    1,0,0,&
    0,1,0,&
    0,0,1/), (/3,3/))

  interface doubleContraction
    procedure doubleContractionBetweenTwoMatrix
    procedure doubleContractionBetweenTwoVector
    procedure doubleContractionSameMatrix
  end interface

  interface T4ijkl_T2kl
    module procedure T4ijkl_T2klReal
    module procedure T4ijkl_T2klT2Complex
  end interface

  interface computeVMEquivalentStress
    procedure computeVMEquivalentStressTensor
    procedure computeVMEquivalentStressBBasis
  end interface

  interface computeDVMEquivalentStressDstress
    procedure computeDVMEquivalentStressDstressTensor
    procedure computeDVMEquivalentStressDStressBBasis
  end interface

  interface computeDVMEquivalentStressDstressiDStressj
    procedure computeDVMEquivalentStressDstressiDStressjBBasis
  end interface

  interface computePressure
    procedure computePressureTensor
    procedure computePressureBBasis
  end interface

  interface computeDPressureDStress
    procedure computeDPressureDStressTensor
    procedure computeDPressureDStressBBasis
  end interface

  interface normalizeVector
    procedure normalizeVectorOutOfPlace
    procedure normalizeVectorInPlace
  end interface

contains
  
  subroutine getIdentity( I_mtx)
    real(k_real), intent(out), dimension(:,:) :: I_mtx
    integer :: n_row, idx
    n_row = size(I_mtx,1)
    I_mtx = 0._k_real
    do idx=1, n_row
      I_mtx(idx,idx) = 1._k_real
    enddo
  end subroutine

  subroutine normalizeVectorOutOfPlace(vector, unit_vector)
    real(k_real), intent(in), dimension(:) :: vector
    real(k_real), intent(out), dimension(:) :: unit_vector

    unit_vector = vector/vectorNorm(vector)
  end subroutine

  subroutine normalizeVectorInPlace(vector)
    real(k_real), intent(inout), dimension(:) :: vector
    real(k_real) :: vnorm
    vnorm = vectorNorm(vector)
    vector = vector/vnorm
  end subroutine

  subroutine vectorCrossProduct(v1, v2, v_cross)
    real(k_real), intent(in), dimension(3) :: v1, v2
    real(k_real), intent(out), dimension(3) :: v_cross

    v_cross(1) = v1(2)*v2(3) - v1(3)*v2(2)
    v_cross(2) = v1(3)*v2(1) - v1(1)*v2(3)
    v_cross(3) = v1(1)*v2(2) - v1(2)*v2(1)
  end subroutine

  function computeTrace(tensor1) result(trace)
    implicit none
    real(k_real), intent(in), dimension(3,3) :: tensor1
    real(k_real) :: trace
    integer :: i

    trace = 0
    do i= 1,3
      trace=trace+tensor1(i,i)
    enddo

  end function

  function computePressureTensor(tensor1) result(pressure)
    implicit none
    real(k_real), intent(in), dimension(3,3) :: tensor1
    real(k_real) :: pressure
    pressure = computeTrace(tensor1)/3._k_real
  end function

  function computeDPressureDstressTensor(stress33) result(dsh_dsigma_val)
    implicit none
    real(k_real), intent(in), dimension(3,3), target :: stress33
    real(k_real), parameter, dimension(3,3) :: dsh_dsigma = reshape((/&
      1._k_real/3._k_real,0._k_real ,0._k_real ,&
      0._k_real, 1._k_real/3._k_real,0._k_real ,&
      0._k_real, 0._k_real ,1._k_real/3._k_real /), (/3,3/))
    real(k_real) :: dsh_dsigma_val(3,3)
    __DECL_UNUSED_MATRIX_PTR__

    dsh_dsigma_val = dsh_dsigma

    __SUPPRESS_UNUSED_MATRIX_WARNING__(stress33)
  end function

  function getDeviatoricTensor(tensor1) result(dev_tensor)
    implicit none
    real(k_real), intent(in), dimension(3,3) :: tensor1
    real(k_real) :: dev_tensor(3,3), pressure
    integer :: i, j

    pressure = computePressure(tensor1)

    do i= 1,3
      do j= 1,3
        dev_tensor(i,j)=tensor1(i,j)
        if (i==j) dev_tensor(i,i)=dev_tensor(i,i)-pressure
      enddo
    enddo

  end

  function dDevTensor2dTensor2(tensor2) result(ddevtensor_dtensor)
    implicit none
    real(k_real), intent(in) :: tensor2(3,3)
    real(k_real) :: ddevtensor_dtensor(3,3,3,3)
    real(k_real), parameter, dimension(3,3,3,3) :: dT2ij_dT2kl=reshape((/&
      1,0,0, 0,0,0, 0,0,0, &
      0,1,0, 0,0,0, 0,0,0, &
      0,0,1, 0,0,0, 0,0,0, &

      0,0,0, 1,0,0, 0,0,0, &
      0,0,0, 0,1,0, 0,0,0, &
      0,0,0, 0,0,1, 0,0,0, &

      0,0,0, 0,0,0, 1,0,0, &
      0,0,0, 0,0,0, 0,1,0, &
      0,0,0, 0,0,0, 0,0,1 /), (/3,3,3,3/))

    ddevtensor_dtensor = dT2ij_dT2kl - tensor2OuterProduc(I33,computeDPressureDstress(tensor2))

  end function

  function T1ij_T1kl_T2ijkl(tensor1, tensor2) result(scalar)
    implicit none
    real(k_real), intent(in), dimension(3,3,3,3) :: tensor2
    real(k_real), intent(in), dimension(3,3) :: tensor1
    real(k_real) :: scalar
    integer :: i,j,k,l
    scalar  =0._k_real
    do k=1,3
      do l=1,3
        do j=1,3
          do i=1,3
            scalar = scalar + tensor1(i,j)*tensor1(k,l)*tensor2(i,j,k,l)
          enddo
        enddo
      enddo
    enddo
  end function

  function doubleContractionT1ijklT2kl(tensor1, tensor2) result(t2out)
    implicit none
    real(k_real), intent(in), dimension(3,3,3,3) :: tensor1
    real(k_real), intent(in), dimension(3,3) :: tensor2
    real(k_real) :: t2out(3,3)
    integer :: i,j

    do j=1,3
      do i=1,3
        t2out(i,j) = doubleContraction(tensor1(i,j,:,:), tensor2)
      enddo
    enddo
  end function

  function doubleContractionT1ijT2ijkl(tensor1, tensor2) result(t2out)
    implicit none
    real(k_real), intent(in), dimension(3,3) :: tensor1
    real(k_real), intent(in), dimension(3,3,3,3) :: tensor2

    real(k_real) :: t2out(3,3)
    integer :: k,l

    do l=1,3
      do k=1,3
        t2out(k,l) = doubleContraction(tensor1, tensor2(:,:,k,l))
      enddo
    enddo
  end function

  function doubleContractionBetweenTwoMatrix(matrix1, matrix2) result(val)
    implicit none
    real(k_real), intent(in), dimension(:,:) :: matrix1, matrix2
    real(k_real) :: val

    val = sum(matrix1*matrix2)

  end function

  function doubleContractionBetweenTwoVector(vector1, vector2) result(val)
    implicit none
    real(k_real), intent(in), dimension(:) :: vector1, vector2
    real(k_real) :: val

    val = sum(vector1*vector2)
  end function

  function doubleContractionSameMatrix(matrix1) result(val)
    implicit none
    real(k_real), intent(in), dimension(:,:) :: matrix1
    real(k_real) :: val

    val = sum(matrix1*matrix1)
  end function

  function getSymmetricPart(non_symmetric_tensor) result(sym_tensor)
    implicit none
    real(k_real), intent(in), dimension(:,:) :: non_symmetric_tensor
    real(k_real), dimension(size(non_symmetric_tensor,1),size(non_symmetric_tensor,1)) :: sym_tensor

    sym_tensor = 0.5_k_real*(non_symmetric_tensor+transpose(non_symmetric_tensor))
  end function

  function getSkewPart(non_symmetric_tensor) result(sym_tensor)
    implicit none
    real(k_real), intent(in), dimension(:,:) :: non_symmetric_tensor
    real(k_real), dimension(size(non_symmetric_tensor,1),size(non_symmetric_tensor,1)) :: sym_tensor

    sym_tensor = 0.5_k_real*(non_symmetric_tensor-transpose(non_symmetric_tensor))
  end function

  function T4ijkl_T2klReal(t4,t2) result(t2_res)
    implicit none
    real(k_real), intent(in) :: t4(3,3,3,3), t2(3,3)
    real(k_real):: t2_res(3,3)
    integer :: i,j,k,l
    do i=1,3
      do j=1,3
        t2_res(i,j) = 0._k_real
        do k=1,3
          do l=1,3
            t2_res(i,j) = t2_res(i,j)+t4(i,j,k,l)*t2(k,l)
          enddo
        enddo
      enddo
    enddo
  end function

  function T4ijkl_T2klT2Complex(t4,t2) result(t2_res)
    implicit none
    real(k_real), intent(in) :: t4(3,3,3,3)
    COMPLEX(k_real), intent(in) :: t2(3,3)
    COMPLEX(k_real) :: t2_res(3,3)
    integer :: i,j,k,l
    do i=1,3
      do j=1,3
        t2_res(i,j) = 0
        do k=1,3
          do l=1,3
            t2_res(i,j) = t2_res(i,j)+t4(i,j,k,l)*t2(k,l)
          enddo
        enddo
      enddo
    enddo
  end function

  function T4ijkl_T2jl(t4,t2) result(t2_res)
    implicit none
    real(k_real), intent(in) :: t4(3,3,3,3), t2(3,3)
    real(k_real) :: t2_res(3,3)
    integer :: i,j,k,l
    do i=1,3
      do k=1,3
        t2_res(i,k) = 0
        do j=1,3
          do l=1,3
            t2_res(i,k) = t2_res(i,k)+t4(i,j,k,l)*t2(j,l)
          enddo
        enddo
      enddo
    enddo
  end function

  function T2ik_T2jl(t2_1,t2_2) result(t4_res)
    implicit none
    real(k_real), intent(in) :: t2_1(3,3), t2_2(3,3)
    real(k_real) :: t4_res(3,3,3,3)
    integer :: i,j,k,l
    do i=1,3
      do k=1,3
        do j=1,3
          do l=1,3
            t4_res(i,j,k,l) = t2_1(i,k)*t2_2(j,l)
          enddo
        enddo
      enddo
    enddo
  end function


  function vectorNOuterProduct(v1,v2, N) result(mtxNN)
    implicit none
    integer, intent(in) :: N
    real(k_real), intent(in) :: v1(N), v2(N)
    real(k_real) :: mtxNN(N,N)
    integer :: i,j
    do j=1,N
      do i=1,N
        mtxNN(i,j) = v1(i)*v2(j)
      enddo
    enddo
  end function


  function vector6OuterProduct(v1,v2) result(mtx66)
    implicit none
    integer, parameter :: N=6
    real(k_real), intent(in) :: v1(N), v2(N)
    real(k_real) :: mtx66(N,N)
    mtx66 = vectorNOuterProduct(v1,v2,N)
  end function

  function vector5OuterProduct(v1,v2) result(mtx55)
    implicit none
    integer, parameter :: N=5
    real(k_real), intent(in) :: v1(N), v2(N)
    real(k_real) :: mtx55(N,N)
    mtx55 = vectorNOuterProduct(v1,v2,N)
  end function

  function vector3OuterProduct(v1,v2) result(mtx33)
    implicit none
    integer, parameter :: N=3
    real(k_real), intent(in) :: v1(N), v2(N)
    real(k_real) :: mtx33(N,N)
    mtx33 = vectorNOuterProduct(v1,v2,N)
  end function

  function tensor2Norm(t2) result(norm)
    implicit none
    real(k_real), intent(in) :: t2(3,3)
    real(k_real) :: norm
    integer :: i,j

    norm = 0
    do j=1,3
      do i=1,3
        norm = norm + t2(i,j)*t2(i,j)
      enddo
    enddo
    norm = sqrt(norm)
  end function

  function Weightedtensor2Norm(t2) result(norm)
    implicit none
    real(k_real), intent(in) :: t2(3,3)
    real(k_real) :: norm
    integer :: i,j

    norm = 0
    do j=1,3
      do i=1,3
        norm = norm + t2(i,j)*t2(i,j)
      enddo
    enddo
    norm = sqrt(norm/9._k_real)
  end function

  function tensor2OuterProduc(t1, t2) result(t3333)
    implicit none
    real(k_real), intent(in) :: t1(3,3), t2(3,3)
    real(k_real) :: t3333(3,3,3,3)
    integer :: i,j,k,l

    do l=1,3
      do k=1,3
        do j=1,3
          do i=1,3
        t3333(i,j,k,l) = t1(i,j)*t2(k,l)
        enddo
      enddo
    enddo
  enddo

  end function

  function vectorNorm(v) result(norm)
    implicit none
    real(k_real), intent(in) :: v(:)
    real(k_real) :: norm

    norm = sqrt(sum(v*v))
  end function

  function WeightedvectorNorm(v) result(norm)
    implicit none
    real(k_real), intent(in) :: v(:)
    real(k_real) :: norm
    
    norm = sqrt(sum(v*v)/int2real(size(v,1)))
  end function

! some routines used to compute equivalent stress and strains

function computeVMEquivalentStrainNonSymmetric(non_symmetric_tensor) result(eqStrain)
  implicit none
  real(k_real), intent(in), dimension(3,3) :: non_symmetric_tensor
  real(k_real) :: eqStrain
  real(k_real), dimension(3,3) :: symm_tensor

  symm_tensor = getSymmetricPart(non_symmetric_tensor)
  eqStrain = computeVMEquivalentStrain(symm_tensor)
end function

function computeVMEquivalentStressTensor(stress_tensor) result(eqStress)
  implicit none
  real(k_real), intent(in), dimension(3,3) :: stress_tensor
  real(k_real) :: eqStress
  real(k_real), dimension(3,3) :: deviatoric_tensor

  deviatoric_tensor = getDeviatoricTensor(stress_tensor)
  eqStress = sqrt(1.5_k_real*doubleContraction(deviatoric_tensor))
end function

function computeDVMEquivalentStressDstressTensor(stress_tensor) result(deqstress_dstress)
  real(k_real), intent(in), dimension(3,3) :: stress_tensor
  real(k_real), dimension(3,3) :: deqstress_dstress

  real(k_real) :: eqStress, deqstress_dstressdev(3,3), deviatoric_tensor(3,3)

  deviatoric_tensor = getDeviatoricTensor(stress_tensor)
  eqStress = computeVMEquivalentStress(stress_tensor)
  deqstress_dstressdev = 0._k_real

  deqstress_dstressdev = 3._k_real*deviatoric_tensor/eqStress
  deqstress_dstress = doubleContractionT1ijT2ijkl(deqstress_dstressdev, dDevTensor2dTensor2(stress_tensor))

end function

function computeVMEquivalentStressBBasis(stress_vector) result(eqStress)
  use change_tensor_basis, only : chg_basis_vector6_to_tensor2
  implicit none
  real(k_real), intent(in), dimension(6) :: stress_vector
  real(k_real) :: eqStress

  eqStress = sqrt(1.5_k_real*sum(stress_vector(1:5)**2))

end function

function computeDVMEquivalentStressDStressBBasis(stress_vector) result(deqStress_dstress)
  use change_tensor_basis, only : dtensor2_dvector6
  implicit none
  real(k_real), intent(in), dimension(6) :: stress_vector
  real(k_real), dimension(6) :: deqStress_dstress
  real(k_real) :: eqStress

  eqStress = computeVMEquivalentStress(stress_vector)
  deqStress_dstress = 0._k_real
  if (eqStress.ne.0._k_real) &
    deqStress_dstress(1:5) = 1.5_k_real*stress_vector(1:5)/eqStress

end function

function computeDVMEquivalentStressDstressiDStressjBBasis(stress_vector) result(deqStress_dstressi_dstressj)
  use change_tensor_basis, only : dtensor2_dvector6
  implicit none
  real(k_real), intent(in), dimension(6) :: stress_vector
  real(k_real), dimension(6,6) :: deqStress_dstressi_dstressj
  real(k_real) :: eqStress, deqStress_dsi(6)
  integer :: i,j 
  
  eqStress = computeVMEquivalentStress(stress_vector)
  deqStress_dsi = computeDVMEquivalentStressDStress(stress_vector)
  if (eqStress.ne.0._k_real) then
    do j =1,6
      do i =1,6
        deqStress_dstressi_dstressj(i,j) = 0._k_real
        if((i.le.5.).and.(j.le.5)) then
          if (i.eq.j) deqStress_dstressi_dstressj(i,j) = deqStress_dstressi_dstressj(i,j) + eqStress
          deqStress_dstressi_dstressj(i,j) = deqStress_dstressi_dstressj(i,j) - stress_vector(i)*deqStress_dsi(j)
          deqStress_dstressi_dstressj(i,j) = 1.5_k_real*deqStress_dstressi_dstressj(i,j)/(eqStress**2)
        endif
      enddo
    enddo
  endif
end function

function computePressureBBasis(stress_vector) result(pressure)
  use change_tensor_basis, only : chg_basis_vector6_to_tensor2
  implicit none
  real(k_real), intent(in) :: stress_vector(6)
  real(k_real) :: pressure
  real(k_real), parameter :: rsq3 = 1._k_real/sqrt(3._k_real)

  pressure = stress_vector(6)*rsq3
end function

function computeDPressureDStressBBasis(stress_vector) result(dpressure_dstress_val)
  use change_tensor_basis, only : chg_basis_vector6_to_tensor2
  implicit none
  real(k_real), intent(in), target :: stress_vector(6)
  real(k_real), parameter :: rsq3 = 1._k_real/sqrt(3._k_real)
  real(k_real), parameter, dimension(6) :: dpressure_dstress = (/ &
  0._k_real,0._k_real,0._k_real,0._k_real,0._k_real, rsq3/)
  real(k_real) :: dpressure_dstress_val(6)
  __DECL_UNUSED_VECTOR_PTR__

  dpressure_dstress_val = dpressure_dstress

  __SUPPRESS_UNUSED_VECTOR_WARNING__(stress_vector)
end function

function computeVMEquivalentStrain(strain_tensor) result(eqStrain)
  implicit none
  real(k_real), intent(in), dimension(3,3) :: strain_tensor
  real(k_real) :: eqStrain
  real(k_real), dimension(3,3) :: deviatoric_tensor

  deviatoric_tensor = getDeviatoricTensor(strain_tensor)
  eqStrain = sqrt(2._k_real/3._k_real*doubleContraction(deviatoric_tensor))
end function

! matrix operations
function mat66InnerProdct(mat66_1, mat66_2) result(res)
  implicit none
  real(k_real), intent(in), dimension(6,6) :: mat66_1, mat66_2
  real(k_real), dimension(6,6) :: res
  integer :: i,j,k

  do k=1,6
    do i = 1,6
      res(i,j) = 0
      do j = 1,6
        res(i,k) = res(i,k) + mat66_1(i,j)*mat66_2(j,k)
      enddo
    enddo
  enddo

end function

subroutine rotateTensor4(t4, R_mtx, t4_rotated)
  implicit none
  real(k_real), intent(in) :: t4(3,3,3,3), R_mtx(3,3)
  real(k_real), intent(out) :: t4_rotated(3,3,3,3)
  integer :: i,j,k,l,p,q,r,s
  do i = 1,3; do j = 1,3; do k = 1,3; do l = 1,3
    t4_rotated(i,j,k,l) = 0._k_real
    do p = 1,3; do q = 1,3; do r = 1,3; do s = 1,3
      t4_rotated(i,j,k,l) = t4_rotated(i,j,k,l)+R_mtx(i,p)*R_mtx(j,q)*R_mtx(k,r)*R_mtx(l,s)*t4(p, q, r, s)
    enddo; enddo; enddo; enddo
  enddo; enddo; enddo; enddo
end subroutine
end module
