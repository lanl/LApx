module compute_averages_mod
  use kinds
  use mpi_variables_mod
  use mpi_useful_routines_mod, only :

contains

  subroutine UpdateRateFields()
    use global, only: all_mighty_grid_global, &
                      total_strain, total_strain_rate, disgradtot, disgradtot_old, vel_grad, &
                      edotp, ept, ept_old, eelfield_rate, &
                      stiffness66, &
                      grid_stress, eelfield, &
                      etotav, edottav, &
                      tdot, &
                      dis_grad_macro, vel_grad_macro, &
                      dvme, dvmp, dvmtot, edotpav, epav, evme, evmp, evmtot, &
                      rve_equivalent_total_strain, rve_equivalent_total_strain_rate, &
                      rve_stress

    use matrix_inversion, only: matrixInverseSymmetric
    use change_tensor_basis, only: chg_basis_vector6_to_tensor2, chg_basis_matrix66_to_tensor4
    use tensor_math_mod, only: getSymmetricPart, T4ijkl_T2kl, &
                                computeVMEquivalentStrainNonSymmetric, computeVMEquivalentStrain
    use mpi_useful_routines_mod, only : MPIAverageGridMatrix, MPISumMatrix
    implicit none
    integer :: ix, iy, iz
    real(k_real) :: S66_temp(6,6), S3333_temp(3,3,3,3)
    integer ::  x_start_rank, x_end_rank, y_start_rank, y_end_rank, z_start_rank, z_end_rank
  
  
    call all_mighty_grid_global%AMGGetLoopLimitsRank(x_start_rank, x_end_rank, y_start_rank, y_end_rank, z_start_rank, z_end_rank)

    DO iz=z_start_rank, z_end_rank
      DO iy=y_start_rank, y_end_rank
        DO  ix=x_start_rank, x_end_rank
          disgradtot(:,:,ix,iy,iz) = disgradtot_old(:,:,ix,iy,iz) + vel_grad(:,:,ix,iy,iz)*tdot
          ! compute the taol strain  from displacements
          total_strain(:,:,ix,iy,iz) = getSymmetricPart(disgradtot(:,:,ix,iy,iz))
          total_strain_rate(:,:,ix,iy,iz) = getSymmetricPart(vel_grad(:,:,ix,iy,iz))

          ept(:,:,ix,iy,iz) = ept_old(:,:,ix,iy,iz)+edotp(:,:,ix,iy,iz)*tdot

          ! compute elastic strain field
          call matrixInverseSymmetric(stiffness66(:,:,ix,iy,iz), S66_temp)
          call chg_basis_matrix66_to_tensor4(S66_temp, S3333_temp)
          eelfield(:,:,ix,iy,iz) = T4ijkl_T2kl(S3333_temp, grid_stress(:,:,ix,iy,iz))
          eelfield_rate(:,:,ix,iy,iz) = total_strain_rate(:,:,ix,iy,iz) - edotp(:,:,ix,iy,iz)
        end do
      end do
    end do

    ! compute averages
    call MPIAverageGridMatrix(edotp, edotpav)
    call MPIAverageGridMatrix(vel_grad, vel_grad_macro)
    call MPIAverageGridMatrix(grid_stress, rve_stress)
    call MPIAverageGridMatrix(disgradtot, dis_grad_macro)
    call MPIAverageGridMatrix(ept, epav)
    call MPIAverageGridMatrix(total_strain, etotav)
    call MPIAverageGridMatrix(total_strain_rate, edottav)

    ! compute effective RVE quantities
    rve_equivalent_total_strain_rate = computeVMEquivalentStrainNonSymmetric(vel_grad_macro) ! strain rate
    rve_equivalent_total_strain = computeVMEquivalentStrainNonSymmetric(dis_grad_macro) ! strain

    evmp= computeVMEquivalentStrain(epav)
    dvmp= computeVMEquivalentStrain(edotpav)
    evmtot= computeVMEquivalentStrain(etotav)
    dvmtot= computeVMEquivalentStrain(edottav)
    evme= computeVMEquivalentStrain(etotav-epav)
    dvme= computeVMEquivalentStrain(edottav-edotpav)

  end subroutine

  subroutine StressGuessNextOuterLoop()
    use global, only: all_mighty_grid_global, &
                      eelfield, stiffness66, grid_stress
    use tensor_math_mod, only: T4ijkl_T2kl
    use change_tensor_basis, only : chg_basis_matrix66_to_tensor4
    implicit none
    real(k_real) :: Tensor4(3,3,3,3)
    integer :: ix, iy, iz
    integer ::  x_start_rank, x_end_rank, y_start_rank, y_end_rank, z_start_rank, z_end_rank
    
    call all_mighty_grid_global%AMGGetLoopLimitsRank(x_start_rank, x_end_rank, y_start_rank, y_end_rank, z_start_rank, z_end_rank)
    DO iz=z_start_rank, z_end_rank; DO iy=y_start_rank, y_end_rank; DO  ix=x_start_rank,x_end_rank
        call chg_basis_matrix66_to_tensor4(stiffness66(:,:,ix,iy,iz), Tensor4)
        grid_stress(:,:,ix,iy,iz) = T4ijkl_T2kl(Tensor4, eelfield(:,:,ix,iy,iz))
    END DO; END DO; END DO

  end subroutine
end module
