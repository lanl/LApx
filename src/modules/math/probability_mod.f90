module probability_mod
use kinds

contains
  subroutine gaussianProbabilityFromMeanAndStdDev(mean, std_dev, x, P)
    use math_constants, only : PI
    implicit none
    real(k_real), intent(in) :: mean, std_dev, x
    real(k_real), intent(out) :: P
    associate (exp_arg=> -0.5_k_real * ( (x-mean)/std_dev)**2 )
      if (exp_arg.lt.-210._k_real) then ! preventing underflow
        P = 0._k_real
      else
        P = exp( exp_arg ) /(std_dev*sqrt(2._k_real*PI))
      endif
    end associate
  end subroutine

end module
