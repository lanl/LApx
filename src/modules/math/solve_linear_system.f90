module solve_linear_system
  use kinds
  implicit none
contains
  subroutine solve_Ax_eq_b(A,b,x)
    real(k_real), dimension(:,:), intent(in) :: A
    real(k_real), dimension(size(A,1)), intent(in) :: b
    real(k_real), dimension(size(A,1)), intent(inout) :: x
    integer, dimension(size(A,1)) :: ipiv   ! pivot indices
    integer :: info, lda, ldb, nrhs, n

    external dgesv

    n = size(A,1)
    nrhs = 1 ! number of right hand sides in b
    lda = n  ! leading dimension of a
    ldb = n  ! leading dimension of b
    x = b
    call dgesv(n, nrhs, a, lda, ipiv, x, ldb, info)

    if (info /= 0) then
       stop 'Sovle linear system failed!'
    end if
  end subroutine

  subroutine solve_mixed_linear_system(ibc_strain,strain,ibc_stress,stress,S, sol_strain, sol_stress)
    ! this subroutine solves a linear system A*x=b where some components of x and b are known
    ! we can think at A is the compliance matrix S, x the strain and b the stress
    use voigt_indicial_conversion_mod, only : Tensor4ToMatrixVoigt, &
                                              Tensor2ToVectorVoigt, VectorVoigtToTensor2
    implicit none
    integer, intent(in), dimension(6) ::  ibc_strain, ibc_stress
    real(k_real), intent(in) :: S(3,3,3,3)
    real(k_real), intent(in), dimension(3,3) :: strain, stress
    real(k_real), intent(out), dimension(3,3) :: sol_strain, sol_stress
    real(k_real), dimension(6) :: strain_6, stress_6, mixed_b, mixed_x, sol_strain_6, sol_stress_6
    real(k_real), dimension(6,6) :: S66, mixed_C
    integer :: i,j
    real(k_real), parameter, dimension(6) :: profac=(/1._k_real,&
                                                      1._k_real,&
                                                      1._k_real,&
                                                      2._k_real,&
                                                      2._k_real,&
                                                      2._k_real/)

    ! convert inputs from tensor to voigt notation
    strain_6 = Tensor2ToVectorVoigt(strain)
    stress_6 = Tensor2ToVectorVoigt(stress)
    S66= Tensor4ToMatrixVoigt(S)

    mixed_b = -strain_6*ibc_strain

    do i=1,6
      mixed_b(i)=-1.d0*ibc_strain(i)*strain_6(i)
      do j=1,6
        mixed_b(i)=mixed_b(i)+S66(i,j)*ibc_stress(j)*stress_6(j)*profac(j)
        mixed_C(i,j)=ibc_stress(j)*(i/j)*(j/i)-ibc_strain(j)*S66(i,j)*profac(j)
      enddo
    enddo

    ! solve the mixed system
    call solve_Ax_eq_b(mixed_C, mixed_b, mixed_x)

    ! compute the solution vector
    do i=1,6
      sol_strain_6(i)=ibc_strain(i)*strain_6(i)+ibc_stress(i)*mixed_x(i)
      sol_stress_6(i)=ibc_stress(i)*stress_6(i)+ibc_strain(i)*mixed_x(i)
    enddo

    ! convert back from voigt to cartesian
    sol_strain = VectorVoigtToTensor2(sol_strain_6)
    sol_stress = VectorVoigtToTensor2(sol_stress_6)

  end

  subroutine solve_mixed_linear_system_EVPBC(average_rve_strain_rate_in, average_rve_stress_in, &
                              imposed_velgrad_components, imposed_stress_components, &
                              imposed_total_strain_rate, imposed_stress_bc, &
                              XMSEC3333, strain_rate_solution_33, stress_solution_33)
    use voigt_indicial_conversion_mod
    use tensor_math_mod
    integer, intent(in) :: imposed_velgrad_components(6), imposed_stress_components(6) !-> imposed stress and displacement components
    real(k_real), dimension(3,3), intent(in) :: average_rve_strain_rate_in, average_rve_stress_in
    real(k_real), dimension(3,3,3,3), intent(in) :: XMSEC3333
    real(k_real), dimension(3,3), intent(in) :: imposed_total_strain_rate , imposed_stress_bc
    real(k_real), dimension(3,3), intent(out) :: stress_solution_33, strain_rate_solution_33

    ! temporary variables
    real(k_real) :: delta_average_constituive_strain_6(6), delta_average_constituive_strain_33(3,3)
    real(k_real) :: A_MATRIX(6,6),B_mixed(6), X_mixed(6), error_strain_6(6)
    real(k_real), dimension(6) :: imposed_stress_bc_6,imposed_strain_rate_6, strain_rate_solution_6, stress_solution_6
    real(k_real) :: XMSEC66(6,6)
    real(k_real), parameter, dimension(6) :: profac=(/1._k_real,1._k_real,1._k_real,2._k_real,2._k_real,2._k_real/) !-> these are voigt factors for strains
    INTEGER ::i,j

  ! transform stress and dispalcments in voigt notation
  imposed_strain_rate_6 = Tensor2ToVectorVoigt(imposed_total_strain_rate)
  imposed_stress_bc_6 = Tensor2ToVectorVoigt(imposed_stress_bc)

  delta_average_constituive_strain_33 = average_rve_strain_rate_in - T4ijkl_T2kl(XMSEC3333, average_rve_stress_in)
  delta_average_constituive_strain_6 = Tensor2ToVectorVoigt(delta_average_constituive_strain_33)
  error_strain_6=imposed_strain_rate_6-delta_average_constituive_strain_6



  ! solving the mix problem
  XMSEC66 = Tensor4ToMatrixVoigt(XMSEC3333)
  ! the known mixed lienar_system term B_mixed
  B_mixed = matmul(XMSEC66, imposed_stress_components*imposed_stress_bc_6*profac) -imposed_velgrad_components*error_strain_6
  ! compute the mix compliance matrix A_MATRIX
  do i=1,6
    do j=1,6
      A_MATRIX(i,j)=imposed_stress_components(j)*I66(i,j)-imposed_velgrad_components(j)*XMSEC66(i,j)*profac(j)
    enddo
  enddo

  call solve_Ax_eq_b(A_MATRIX,B_mixed,X_mixed) !-> solve the linear system

  ! X_mixed contains strain solution for imposed stress components and stress solutions for the imposed strain components

  ! compute the solution
  stress_solution_6=imposed_stress_components*imposed_stress_bc_6+imposed_velgrad_components*X_mixed
  strain_rate_solution_6=imposed_velgrad_components*error_strain_6+imposed_stress_components*X_mixed+delta_average_constituive_strain_6

  strain_rate_solution_33 = VectorVoigtToTensor2(strain_rate_solution_6)
  stress_solution_33 = VectorVoigtToTensor2(stress_solution_6)

  end subroutine solve_mixed_linear_system_EVPBC


end module
