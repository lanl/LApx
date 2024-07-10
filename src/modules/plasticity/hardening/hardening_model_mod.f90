module hardening_model_mod
  use kinds
  use all_grid_data_mod, only : all_grid_data

  type :: hardening_model_base
    class(all_grid_data), pointer :: grid_data => null()
    integer :: n_ss
    real(k_real), dimension(:), pointer :: strength_dot
    real(k_real), dimension(:,:,:,:), pointer :: strength_old_grid => null()
    real(k_real), dimension(:), pointer :: strength_old_ptr  => null()

  contains
    procedure :: initHardeningModelBase !-> base function initializing grid_data pointers and allocating the required space
    procedure :: initHardeningVectorData => initHardeningVectorDataBase !-> base function initializing grid_data pointers and allocating the required space
    procedure :: computeStrengthDot => computeStrengthDotHardeningModelBase
    procedure :: getUpdateStrength => getUpdateStrengthModelBase
    procedure :: setPointData => setPointDataHardeningModelBase
  end type hardening_model_base


  type, extends(hardening_model_base) :: voce_hardening_model
    real(k_real), dimension(:,:,:), pointer :: accumulated_gamma_old_grid => null()
    real(k_real), dimension(:,:,:,:), pointer :: gamma_dot_grid => null()

    real(k_real), pointer :: accumulated_gamma_old_ptr
    real(k_real), dimension(:), pointer :: gamma_dot_ptr

    ! model parameters
    real(k_real), pointer :: theta0, theta1, tau0, tau1, m_exponet
  contains
    procedure :: initGeneralizedVoceHardening!-> base function initializing grid_data pointers and allocating the required space
    procedure :: initHardeningVectorData => initGeneralizedVoceHardeningFromVector
    procedure :: computeStrengthDot => computeTauDotGeneralizedVoceHardening
    procedure :: setPointData => setPointDataGeneralizedVoceHardening
  end type voce_hardening_model

  ! interface for generalized_voce_hardening
  interface generalized_voce_hardening
    module subroutine initGeneralizedVoceHardening(this, grid_data, n_ss, theta0, theta1, tau0, tau1, m_exponet)
      class(voce_hardening_model), intent(inout) :: this
      class(all_grid_data), intent(in), target :: grid_data
      integer, intent(in) :: n_ss
      real(k_real), intent(in), target :: theta0, theta1, tau0, tau1, m_exponet
    end subroutine

    module subroutine initGeneralizedVoceHardeningFromVector(this, grid_data, param_integer, param_real)
      class(voce_hardening_model), intent(inout) :: this
      class(all_grid_data), intent(in), target :: grid_data
      real(k_real), intent(in), dimension(:) :: param_real
      integer, intent(in), dimension(:) :: param_integer
    end subroutine

    module subroutine setPointDataGeneralizedVoceHardening(this, ix,iy,iz)
      class(voce_hardening_model), intent(inout) :: this
      integer, intent(in) :: ix,iy,iz
    end subroutine

    module subroutine computeTauDotGeneralizedVoceHardening(this)
      class(voce_hardening_model), intent(inout) :: this
    end subroutine
  end interface generalized_voce_hardening

contains

  subroutine initHardeningVectorDataBase(this, grid_data, param_integer, param_real)
    class(hardening_model_base), intent(inout) :: this
    class(all_grid_data), intent(in), target :: grid_data
    real(k_real), intent(in), dimension(:) :: param_real
    integer, intent(in), dimension(:) :: param_integer

    associate( theta0 => param_real(1), &
               n_ss => param_integer(1))

      call this%initHardeningModelBase(grid_data, n_ss)
    end associate
  end subroutine

  subroutine initHardeningModelBase(this, grid_data, n_ss)
    class(hardening_model_base), intent(inout) :: this
    class(all_grid_data), intent(in), target :: grid_data
    integer, intent(in) :: n_ss

    if (n_ss<1) error stop "initHardeningModelBase n_ss <1. Abort!"
    this%grid_data => grid_data
    this%n_ss = n_ss
    allocate(this%strength_dot(n_ss))
  end subroutine

  subroutine setPointDataHardeningModelBase(this, ix, iy, iz)
    class(hardening_model_base), intent(inout) :: this
    integer, intent(in) :: ix,iy,iz

    this%strength_old_ptr => this%strength_old_grid(:,ix,iy,iz)
  end subroutine

  subroutine computeStrengthDotHardeningModelBase(this)
    class(hardening_model_base), intent(inout) :: this

    error stop "If you end up here it means you did not override computeStrengthDotHardeningModelBase in your class"
    this%strength_dot = 0
  end subroutine

  subroutine getUpdateStrengthModelBase(this, ix, iy, iz, dt, new_ss_strength)
    class(hardening_model_base), intent(inout) :: this
    integer, intent(in) :: ix,iy,iz
    real(k_real), intent(out), dimension(:) :: new_ss_strength
    real(k_real), intent(in) :: dt
    call this%setPointData(ix, iy, iz)

    call this%computeStrengthDot()
    new_ss_strength = this%strength_old_ptr + this%strength_dot*dt

  end subroutine

  subroutine computeDUpdatedStrengthDStressHardeningModelBase(this)
    class(hardening_model_base), intent(inout) :: this
    error stop "computeDUpdatedStrengthDStressHardeningModelBase you forgot to override computeUpdatedStrength"

    ! here to avoid unsued dummy argument warning
    write(*,*) this%n_ss
  end subroutine

end module hardening_model_mod
