module bc_mod
use kinds
implicit none
  private :: computeInitialGuessMixedBC
  public :: computeInitialGuess
contains

  subroutine computeInitialGuess()
    implicit none

    call computeInitialGuessMixedBC()

  end subroutine

  subroutine computeInitialGuessMixedBC()
    use global, only : vel_grad, grid_stress, grid_stress_old, stiffness66, eelfield, eelfield_old, eelfield_rate, &
                       rve_stress_vm, rve_equivalent_total_strain_rate, rve_stress, &
                       vel_grad_macro, tdot, all_mighty_grid_global
    use elasticity_mod, only : computeStressFromElasticStrain
    use mpi_useful_routines_mod, only : MPIAverageGridMatrix
    use tensor_math_mod, only : computeVMEquivalentStrainNonSymmetric, computeVMEquivalentStress, getSymmetricPart
    use bc_objects_mod, only : i_scauchy_rate, i_vel_grad, i_scauchy_components, i_vel_grad_comp
    use change_tensor_basis, only : chg_basis_matrix66_to_tensor4
    use voigt_indicial_conversion_mod, only : Tensor2ToVectorVoigt
    use matrix_inversion, only : lapackInverseSymmetric
    use solve_linear_system, only : solve_mixed_linear_system
    implicit none
    integer :: ix, iy, iz
    real(k_real) :: S(3,3,3,3), S66(6,6), i_vel_grad_symm(3,3), stress_rate_sol(3,3)
    integer, dimension(6) :: i_vel_grad_comp_6,  i_scauchy_components_6
    integer :: x_start, x_end, y_start, y_end, z_start, z_end

    i_scauchy_components_6 =Tensor2ToVectorVoigt(i_scauchy_components)
    i_vel_grad_comp_6 =Tensor2ToVectorVoigt(i_vel_grad_comp)
    i_vel_grad_symm = getSymmetricPart(i_vel_grad)

    call all_mighty_grid_global%AMGGetLoopLimitsRank(x_start, x_end, y_start, y_end, z_start, z_end)

    do iz=z_start, z_end
      do iy=y_start, y_end
        do ix=x_start, x_end

          call lapackInverseSymmetric(stiffness66(:,:,ix,iy,iz), S66)
          call chg_basis_matrix66_to_tensor4(S66, S)

          call solve_mixed_linear_system(i_vel_grad_comp_6, i_vel_grad_symm, i_scauchy_components_6,i_scauchy_rate, S, &
                eelfield_rate(:,:,ix,iy,iz), stress_rate_sol)

          vel_grad(:,:,ix,iy,iz) = eelfield_rate(:,:,ix,iy,iz)
          eelfield(:,:,ix,iy,iz) = eelfield_old(:,:,ix,iy,iz) + eelfield_rate(:,:,ix,iy,iz)*tdot
          grid_stress(:,:,ix,iy,iz) = grid_stress_old(:,:,ix,iy,iz) + stress_rate_sol*tdot
          end do
       end do
    end do
    call MPIAverageGridMatrix(grid_stress, rve_stress)
    rve_stress_vm = computeVMEquivalentStress(rve_stress)
    call MPIAverageGridMatrix(vel_grad, vel_grad_macro)
    rve_equivalent_total_strain_rate = computeVMEquivalentStrainNonSymmetric(vel_grad_macro)
  end subroutine


end module
