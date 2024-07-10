submodule (hardening_model_mod) voce_hardening_model_smod

contains
  module subroutine initGeneralizedVoceHardening(this, grid_data, n_ss, theta0, theta1, tau0, tau1, m_exponet)
    class(voce_hardening_model), intent(inout) :: this
    class(all_grid_data), intent(in), target :: grid_data
    integer, intent(in) :: n_ss
    real(k_real), intent(in), target :: theta0, theta1, tau0, tau1, m_exponet

    call this%initHardeningModelBase(grid_data, n_ss)

    this%theta0 => theta0
    this%theta1 => theta1
    this%tau0 => tau0
    this%tau1 => tau1
    this%m_exponet => m_exponet

    call this%grid_data%getSSScalarDataPointerByName("gamma_dot", this%gamma_dot_grid)
    call this%grid_data%getSSScalarDataPointerByName("crss", &
      this%strength_old_grid)
    call this%grid_data%getScalarDataPointerByNameOld("total_accumulated_slip", &
        this%accumulated_gamma_old_grid)
  end subroutine

  module subroutine initGeneralizedVoceHardeningFromVector(this, grid_data, param_integer, param_real)
    class(voce_hardening_model), intent(inout) :: this
    class(all_grid_data), intent(in), target :: grid_data
    real(k_real), intent(in), dimension(:) :: param_real
    integer, intent(in), dimension(:) :: param_integer

        associate( theta0 => param_real(1), &
               theta1 => param_real(2), &
               tau0 =>   param_real(3), &
               tau1 =>   param_real(4), &
               m_exponet => param_real(5), &
               n_ss => param_integer(1))

    call this%initGeneralizedVoceHardening(grid_data, n_ss, theta0, theta1, tau0, tau1, m_exponet)
    end associate
  end subroutine

  module subroutine setPointDataGeneralizedVoceHardening(this, ix,iy,iz)
    class(voce_hardening_model), intent(inout) :: this
    integer, intent(in) :: ix,iy,iz

    call this%hardening_model_base%setPointData(ix,iy,iz)
    this%accumulated_gamma_old_ptr => this%accumulated_gamma_old_grid(ix,iy,iz)
    this%gamma_dot_ptr => this%gamma_dot_grid(:,ix,iy,iz)
  end subroutine

  module subroutine computeTauDotGeneralizedVoceHardening(this)
    class(voce_hardening_model), intent(inout) :: this
    real(k_real) :: d_gamma, g_accum, exp_arg

    associate(g_accum_old => this%accumulated_gamma_old_ptr, &
              g_dot => this%gamma_dot_ptr, &
              theta0 => this%theta0, &
              theta1 => this%theta1, &
              tau0 => this%tau0, &
              tau1 => this%tau1, &
              n_ss =>this%n_ss, &
              tau_dot => this%strength_dot)

    d_gamma = sum(abs(g_dot(1:n_ss)))
    g_accum= g_accum_old + d_gamma
    exp_arg = g_accum*theta0/tau1


    tau_dot = d_gamma*(theta1+exp(-exp_arg)*(theta0+theta1*(exp_arg-1)))

    end associate
  end subroutine



end submodule
