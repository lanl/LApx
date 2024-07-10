module solve_non_linear_system
  use kinds
  implicit none

real(k_real), parameter, private :: nr_tol = 1e-12_k_real
integer, parameter, private  :: n_iter_max = 50

type, abstract :: nr_single_variable
    real(k_real) :: nl_abs_tol = 0
    integer      :: max_nl_iter = 0
    logical      :: initialized = .false.
    real(k_real) :: Residual, dResidual_dx ! resdiual and derivatives
    real(k_real) :: F, dF_dx ! function value and derivatives
    real(k_real) :: x        ! the current x value
    integer      :: iter, &
                    n_params = 0
    real(k_real), pointer, dimension(:) :: dx_dP ! derivative of the function with respect of the parameters
    real(k_real), pointer, dimension(:) :: dF_dP ! derivative of the function with respect of the parameters
    real(k_real), pointer, dimension(:) :: dResidual_dP ! derivative of the residual with respect of the parameters

contains
  procedure(computeFvalAndDFvalDxInterface), deferred :: computeFvalAndDFvalDx
  procedure :: computeResidualAndJacobian =>  computeResidualAndJacobianNRSingleVariable
  procedure :: solve => solveNRSingleVariable
  procedure :: init => initNRSingleVariable
end type

abstract interface
  subroutine computeFvalAndDFvalDxInterface(this)
    use kinds
    import nr_single_variable
    class(nr_single_variable), intent(inout) :: this
  end subroutine
end interface

contains

  subroutine initNRSingleVariable(this, tol_in, max_iter_in, n_params)
    class(nr_single_variable), intent(inout) :: this
    real(k_real), intent(in), optional :: tol_in
    integer, intent(in), optional :: max_iter_in, &
                                     n_params

    !if we provide the tolerance, use it
    if (present(tol_in)) then
      this%nl_abs_tol = tol_in
    else
      this%nl_abs_tol = nr_tol
    endif

    if (this%nl_abs_tol<= 0) error stop "initNRSingleVariable: nl_abs_tol<=0, Abort!"

    !if we provide the max number  of iteration use it
    if (present(max_iter_in)) then
      this%max_nl_iter = max_iter_in
    else
      this%max_nl_iter = n_iter_max
    endif

    if (present(n_params)) then
      if (n_params <= 0) error stop "initNRSingleVariable: number of parameters <= 0, Abort!"
      this%n_params = n_params
      allocate(this%dx_dP(n_params))
      allocate(this%dF_dP(n_params))
      allocate(this%dResidual_dP(n_params))
    endif

    if (this%max_nl_iter<= 0) error stop "initNRSingleVariable: max_nl_iter<=0, Abort!"

    this%initialized = .true.

  end subroutine

  subroutine computeResidualAndJacobianNRSingleVariable(this)
    class(nr_single_variable), intent(inout) :: this

    call this%computeFvalAndDFvalDx()

    associate(R => this%Residual, &
              dR_dx => this%dResidual_dx, &
              F => this%F, &
              dF_dx => this%dF_dx)

    R = F**2
    dR_dx = 2*F*dF_dx

    end associate
  end subroutine

  subroutine solveNRSingleVariable(this, x0, x_sol, dx_dP_out)
    class(nr_single_variable), intent(inout) :: this
    real(k_real), intent(in) :: x0
    real(k_real), intent(out) :: x_sol
    real(k_real), intent(out), dimension(this%n_params), optional :: dx_dP_out
    associate(x => this%x, &
              dx_dP => this%dx_dP, &
              iter => this%iter, &
              R => this%Residual, &
              dR_dx => this%dResidual_dx, &
              F => this%F, &
              dF_dP => this%dF_dP, &
              dF_dx => this%dF_dx, &
              dR_dP => this%dResidual_dP)

    x = x0
    iter = 0
    call this%computeResidualAndJacobian()

    do while ((R.gt.this%nl_abs_tol).and.(iter<=this%max_nl_iter))
      ! update x guess
      x = x - R/dR_dx
      ! compute residual and jacobian
      call this%computeResidualAndJacobian()
      iter = iter + 1
    enddo

    if (R>this%nl_abs_tol) error stop "NR did not converge, abort"
    x_sol = x

    !compute the derivatives of X w.r.t. parameters
    if (this%n_params>0) then
      dR_dP = 2*F*dF_dP
      dx_dP = -dR_dP/dR_dx
    endif

    if (present(dx_dP_out)) then
      dx_dP_out = dx_dP
    endif
    end associate
  end subroutine

end module
