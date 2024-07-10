module elasticity_mod
  use kinds
  use change_tensor_basis
  use matrix_inversion
  use tensor_math_mod, only : T4ijkl_T2kl
implicit none
contains
  subroutine computeStressFromElasticStrain(c66,e33,s33)
    real(k_real), intent(in) :: c66(6,6), e33(3,3)
    real(k_real), intent(out) :: s33(3,3)
    real(k_real) :: c3333(3,3,3,3)

    call chg_basis_matrix66_to_tensor4(c66, c3333)

    s33 = T4ijkl_T2kl(c3333, e33)

  end subroutine

  subroutine computeElasticStrainFromStress(c66,s33,e33)
    real(k_real), intent(in) :: c66(6,6), s33(3,3)
    real(k_real), intent(out) :: e33(3,3)
    real(k_real) :: s66(6,6), s3333(3,3,3,3)

    call lapackInverseSymmetric(c66, s66)
    call chg_basis_matrix66_to_tensor4(s66, s3333)

    e33 = T4ijkl_T2kl(s3333, s33)

  end subroutine
end module
