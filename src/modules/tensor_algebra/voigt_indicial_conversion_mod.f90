module voigt_indicial_conversion_mod
  use kinds
implicit none

integer, parameter, dimension(6,2) :: ijv=reshape((/1,2,3,2,1,1,&
                                                    1,2,3,3,3,2/),(/6,2/))
private :: ijv

interface Tensor2ToVectorVoigt
  procedure Tensor2ToVectorVoigtReal
  procedure Tensor2ToVectorVoigtInteger
end interface

interface VectorVoigtToTensor2
  procedure VectorVoigtToTensor2Real
  procedure VectorVoigtToTensor2Integer
end interface

interface Tensor4ToMatrixVoigt
  procedure Tensor4ToMatrixVoigtReal
end interface

interface MatrixVoigtToTensor4
  procedure MatrixVoigtToTensor4Real
end interface

contains


  function Tensor2ToVectorVoigtReal(t2) result(v)
    real(k_real), intent(in) :: t2(3,3)
    real(k_real) :: v(6)
    integer :: i

    do i=1,6
      v(i) = t2(ijv(i,1),ijv(i,2))
    enddo
  end function

  function Tensor2ToVectorVoigtInteger(t2) result(v)
    integer, intent(in) :: t2(3,3)
    integer :: v(6)
    integer :: i

    do i=1,6
      v(i) = t2(ijv(i,1),ijv(i,2))
    enddo
  end function

  function VectorVoigtToTensor2Real(v) result(t2)
    real(k_real), intent(in)  :: v(6)
    real(k_real):: t2(3,3)
    integer :: i

    do i=1,3
      t2(ijv(i,1),ijv(i,2)) =v(i)
    enddo
    do i=4,6
      t2(ijv(i,1),ijv(i,2)) =v(i)
      t2(ijv(i,2),ijv(i,1)) =v(i)
    enddo
  end function

  function VectorVoigtToTensor2Integer(v) result(t2)
    integer, intent(in)  :: v(6)
    integer :: t2(3,3)
    integer :: i

    do i=1,3
      t2(ijv(i,1),ijv(i,2)) =v(i)
    enddo
    do i=4,6
      t2(ijv(i,1),ijv(i,2)) =v(i)
      t2(ijv(i,2),ijv(i,1)) =v(i)
    enddo
  end function

  function Tensor4ToMatrixVoigtReal(t4) result(m66)
    real(k_real), intent(in) :: t4(3,3,3,3)
    real(k_real) :: m66(6,6)
    integer :: i,j

    do i=1,6
      do j=1,6
        m66(i,j) = t4(ijv(i,1),ijv(i,2),ijv(j,1),ijv(j,2))
      enddo
    enddo
  end function

  function MatrixVoigtToTensor4Real(m66) result(t4)
    real(k_real), intent(in) :: m66(6,6)
    real(k_real) :: t4(3,3,3,3)
    integer :: i,j, I1, I2, J1, J2
    !TODO this one we might want to bemake a bit more efficient

    do i=1,6
      I1=IJV(I,1)
      I2=IJV(I,2)
      do j=1,6
        J1=IJV(J,1)
        J2=IJV(J,2)

        t4(i1,i2,j1,j2) = m66(i,j)
        t4(i2,i1,j1,j2) = m66(i,j)
        t4(i1,i2,j2,j1) = m66(i,j)
        t4(i2,i1,j2,j1) = m66(i,j)

      enddo
    enddo

  end function

end module
