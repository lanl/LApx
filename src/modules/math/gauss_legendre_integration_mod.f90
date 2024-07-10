module gauss_legendre_integration_mod
use kinds
implicit none

type :: gauss_legendre_integration
  integer :: n_gauss
  real(k_real) :: average, half_range
  real(k_real), dimension(:), pointer  :: x_i, w_i

  real(k_real) :: dhalf_range_dxmin= -0.5_k_real, dhalf_range_dxmax= 0.5_k_real
  real(k_real) :: daverage_dxmin= 0.5_k_real, daverage_dxmax= 0.5_k_real

contains
  procedure :: init => initGaussLegendreComputeWeightsAndAbscissa
  procedure :: reset => resetGaussLegendreBase
  procedure :: integrate => gaussLegendreIntegrate
  procedure, private :: setIntegrationLimits => setGaussLegendreIntegrationLimits
  procedure :: getXiRealSpace => getGaussLegendreXiRealSpaceBase
end type


! type, extends(gauss_legendre_integration) :: gauss_legendre_gaussian
!   real(k_real) :: mean, std_deviation
!   real(k_real), dimension(:), pointer :: xi_real_space => null(), &
!                                          dxi_real_space_dmean => null(), &
!                                          probability => null(), &
!                                          dprobability_dmean => null()
! contains
!   procedure :: init => initGaussLegendreGaussian
!   procedure :: reset => resetGaussLegendreGaussian
!   procedure :: getXiRealSpaceGaussian
!   procedure :: getDXiRealSpaceDMeanGaussian
!   procedure :: integrate => gaussLegendreIntegrateGaussian
!   procedure :: dIntegralDmean => gaussLegendreDIntegrateGaussianDMean
!   procedure :: computeProbability => computeProbabilityGaussian
!   procedure :: computeDProbabilityDMean => computeDProbabilityGaussianDMean
! end type

! type :: gauss_legendre_gaussian_array
!   type (gauss_legendre_gaussian), pointer :: pt(:) => null()
!   integer :: n_gauss, n_gaussian
! contains
!   procedure :: init => initGaussianArray
!   procedure :: reset => resetGaussianArray
! end type

contains

  ! subroutine initGaussianArray(this, n_gauss, n_gaussian)
  !   implicit none
  !   class(gauss_legendre_gaussian_array), intent(inout) :: this
  !   integer, intent(in) :: n_gauss, n_gaussian
  !   integer :: i
  !   this%n_gauss = n_gauss
  !   this%n_gaussian = n_gaussian
  !   allocate(this%pt(n_gaussian))
  !   do i=1,n_gaussian
  !     call this%pt(i)%init(n_gauss)
  !   enddo
  !
  ! end subroutine
  !
  ! subroutine resetGaussianArray(this)
  !   implicit none
  !   class(gauss_legendre_gaussian_array), intent(inout) :: this
  !   integer :: i
  !
  !   do i=1,this%n_gaussian
  !     call this%pt(i)%reset()
  !   enddo
  !   deallocate(this%pt)
  !   nullify(this%pt)
  ! end subroutine
  !
  ! subroutine resetGaussLegendreGaussian(this)
  !   implicit none
  !   class(gauss_legendre_gaussian), intent(inout) :: this
  !   call this%gauss_legendre_integration%reset()
  !
  !   deallocate(this%xi_real_space, &
  !              this%dxi_real_space_dmean, &
  !              this%probability, &
  !              this%dprobability_dmean )
  !   nullify(this%xi_real_space, &
  !           this%dxi_real_space_dmean, &
  !           this%probability, &
  !           this%dprobability_dmean )
  ! end subroutine

  subroutine resetGaussLegendreBase(this)
    implicit none
    class(gauss_legendre_integration), intent(inout) :: this

    deallocate(this%x_i, &
               this%w_i)
    nullify(this%x_i, &
            this%w_i  )
  end subroutine

subroutine initGaussLegendreComputeWeightsAndAbscissa(this, n_gauss)
  use math_constants, only : PI
  implicit none
  class(gauss_legendre_integration), intent(inout) :: this
  integer, intent(in) :: n_gauss
  real(kind=k_real)            :: x, f, df, dx
  integer                 :: i,  iter, k
  real(kind=k_real), allocatable :: p0(:), p1(:), tmp(:)
  integer, parameter :: max_iter = 100

  this%n_gauss = n_gauss
  allocate(this%x_i(n_gauss), this%w_i(n_gauss))

  associate(n=>this%n_gauss, &
            x_i=>this%x_i, &
            w_i=>this%w_i)

  allocate(p0(1))
  allocate(p1(2))
  p0 = (/1._k_real/)
  p1 = (/1._k_real, 0._k_real/)

  if (n>=2) then
  do k = 2, n
     allocate(tmp(k+1))
     tmp = ((2*k-1)*[p1,0._k_real]-(k-1)*[0._k_real, 0._k_real,p0])/k
     deallocate(p0); allocate(p0(size(p1)))
     p0 = p1
     deallocate(p1); allocate(p1(size(tmp)))
     p1 = tmp;
     deallocate(tmp)
  end do
  endif

  do i = 1, n
    x = cos(pi*(i-0.25_k_real)/(n+0.5_k_real))
    do iter = 1, max_iter
      f = p1(1); df = 0._k_real
      do k = 2, size(p1)
        df = f + x*df
        f  = p1(k) + x * f
      end do
      dx =  f / df
      x = x - dx
      if (abs(dx)<10*epsilon(dx)) exit
    end do
    if (iter==max_iter) error stop " initGaussLegendreComputeWeightsAndAbscissa "
    x_i(i) = x
    w_i(i) = 2._k_real/((1._k_real-x**2)*df**2)
  end do

end associate
end subroutine

! subroutine initGaussLegendreGaussian(this, n_gauss)
!   implicit none
!   class(gauss_legendre_gaussian), intent(inout) :: this
!   integer, intent(in) :: n_gauss
!
!   call this%gauss_legendre_integration%init(n_gauss)
!   allocate(this%xi_real_space(n_gauss), &
!            this%dxi_real_space_dmean(n_gauss), &
!            this%probability(n_gauss), &
!            this%dProbability_dmean(n_gauss) )
! end subroutine

subroutine computeHalfRangeAndAverage(x_min, x_max, half_range, average)
  real(k_real), intent(in) :: x_min, x_max
  real(k_real), intent(out) :: half_range, average
  half_range = (x_max-x_min)*.5
  average = (x_max+x_min)*.5
end subroutine

subroutine GaussLegendreIntegrate(this, f_vals, dfvals_dxi, integral_value, dintegralval_dxi)
  class(gauss_legendre_integration), intent(inout) :: this
  real(k_real), intent(in), dimension(this%n_gauss) :: f_vals, dfvals_dxi
  real(k_real), intent(out):: integral_value, dintegralval_dxi

  integral_value = this%half_range*sum(f_vals*this%w_i)
  dintegralval_dxi = this%half_range*sum(dfvals_dxi*this%w_i)
end subroutine

subroutine setGaussLegendreIntegrationLimits(this, x_min, x_max)
  class(gauss_legendre_integration), intent(inout) :: this
  real(k_real), intent(in):: x_min, x_max

  call computeHalfRangeAndAverage(x_min, x_max, this%half_range, this%average)

end subroutine

subroutine getGaussLegendreXiRealSpaceBase(this, x_min, x_max, xi_real_space)
  class(gauss_legendre_integration), intent(inout) :: this
  real(k_real), intent(out), dimension(:) :: xi_real_space
  real(k_real), intent(in):: x_min, x_max

  call this%setIntegrationLimits(x_min, x_max)
  xi_real_space(1:this%n_gauss) = this%half_range * this%x_i + this%average
end subroutine

! subroutine getXiRealSpaceGaussian(this, mean, std_deviation, xi_real_space)
!   class(gauss_legendre_gaussian), intent(inout) :: this
!   real(k_real), intent(in) :: mean, std_deviation
!   real(k_real), intent(out), dimension(:) :: xi_real_space
!   real(k_real), parameter :: N_std_dev = 5._k_real
!   real(k_real) :: dummy
!   real(k_real) :: x_min, x_max
!
!   this%mean = mean
!   this%std_deviation = std_deviation
!
!   ! I don't understand why we are increasing the std_deviation artificially
!   dummy = N_std_dev*std_deviation
!   x_min = mean - dummy
!   x_max = mean + dummy
!
!   call this%setIntegrationLimits(x_min, x_max)
!   this%xi_real_space = this%half_range * this%x_i + this%average
!   xi_real_space = this%xi_real_space
! end subroutine


! subroutine getDXiRealSpaceDMeanGaussian(this, dxi_real_space_dmean)
!   class(gauss_legendre_gaussian), intent(inout) :: this
!   real(k_real), intent(out), dimension(:), optional :: dxi_real_space_dmean
!   real(k_real), parameter :: dxmin_dmean = 1.
!   real(k_real), parameter :: dxmax_dmean = 1.
!
!   associate( dhalf_range_dxmin=> this%dhalf_range_dxmin, &
!              dhalf_range_dxmax=> this%dhalf_range_dxmax, &
!              daverage_dxmin=> this%daverage_dxmin, &
!              daverage_dxmax=> this%daverage_dxmax, &
!              dx_dmean => this%dxi_real_space_dmean  )
!
!   dx_dmean = &
!     dhalf_range_dxmin*dxmin_dmean * this%x_i + &
!     dhalf_range_dxmax*dxmax_dmean * this%x_i + &
!     daverage_dxmin * dxmin_dmean + &
!     daverage_dxmax * dxmax_dmean
!   if (present(dxi_real_space_dmean)) dxi_real_space_dmean = dx_dmean
!
!   end associate
! end subroutine

! subroutine computeProbabilityGaussian(this)
!   use math_constants, only : PI
!   implicit none
!   class(gauss_legendre_gaussian), intent(inout) :: this
!   associate(probability=>this%probability, &
!             mean=>this%mean, &
!             std_deviation=>this%std_deviation, &
!             x=>this%xi_real_space)
!     probability=EXP( -(x-mean)**2/(2.*std_deviation**2) )/(std_deviation*SQRT(2.*PI))
!   end associate
! end subroutine
!
!
! subroutine computeDProbabilityGaussianDMean(this, dprobability_dmean)
!   use math_constants, only : PI
!   implicit none
!   class(gauss_legendre_gaussian), intent(inout) :: this
!   real(k_real), intent(out), dimension(:), optional :: dprobability_dmean
!
!   call this%getDXiRealSpaceDMeanGaussian()
!   associate(P=>this%probability, &
!             mean=>this%mean, &
!             std_deviation=>this%std_deviation, &
!             x=>this%xi_real_space, &
!             dxi_real_space_dmean => this%dxi_real_space_dmean, &
!             dP_dmean => this%dprobability_dmean)
!
!     dP_dmean = -P*(x-mean)/(std_deviation**2)*(dxi_real_space_dmean-1.)
!     if (present(dprobability_dmean)) dprobability_dmean = dP_dmean
!
!   end associate
!
! end subroutine
!
! subroutine gaussLegendreIntegrateGaussian(this, f_vals, integral_value)
!   class(gauss_legendre_gaussian), intent(inout) :: this
!   real(k_real), intent(in), dimension(:) :: f_vals
!   real(k_real), intent(out):: integral_value
!
!   call this%computeProbability()
!   integral_value = this%half_range*sum(f_vals*this%w_i*this%probability)
! end subroutine
!
! subroutine gaussLegendreDIntegrateGaussianDMean(this, f_vals, dfvals_dxi, dintegral_dmean)
!   class(gauss_legendre_gaussian), intent(inout) :: this
!   real(k_real), intent(in), dimension(:) :: f_vals, dfvals_dxi
!   real(k_real), intent(out):: dintegral_dmean
!
!   call this%computeDProbabilityDMean()
!   dintegral_dmean = this%half_range*sum(this%w_i *( &
!                               dfvals_dxi*this%dxi_real_space_dmean*this%probability + &
!                               f_vals*this%dprobability_dmean) )
! end subroutine

end module
