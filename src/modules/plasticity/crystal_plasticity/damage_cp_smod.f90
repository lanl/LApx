submodule(cp_base_mod) damage_cp_mod
  implicit none

#include "macro_debug.fpp"
contains

  module subroutine computeLeblondStressCPBase(this, stress6)
    use tensor_math_mod, only : computeVMEquivalentStress
    use print_utils_mod, only : printToScreen
    implicit none
    class(cp_base), intent(inout) :: this
    real(k_real), intent(in) :: stress6(6)
    real(k_real) :: schmid6(6), lambda
    integer :: ss_idx



    associate(n_ss => this%n_ss, &
              flow_direction => this%drss_dstress, &
              ddirection_dstress => this%ddirection_dstress)

    do ss_idx=1,this%n_ss
      schmid6(1:this%n_schmid_components) = this%schmid_ptr(1:this%n_schmid_components,ss_idx)
      if (this%n_schmid_components.eq.5) schmid6(6) = 0._k_real

      call this%NRLeblondStressCPBase(stress6, schmid6(1:5), this%n_exp_damage, this%porosity_ptr, lambda)
      ! we are done solving now we need to update rss flow direction and dflow direction dstress
      call this%LeblondComputeResolvedStressAndDirection(stress6, schmid6, lambda, this%n_exp_damage, this%porosity_ptr, & ! -> these are all inputs
                            this%rss(ss_idx), flow_direction(:,ss_idx), ddirection_dstress(:,:,ss_idx) )  ! -> these are all outputs

    enddo

    end associate
  end subroutine

  module subroutine computeQ(this, stress6, rss, lambda, schmid, f, q4, Q, dQdsi, dQdLambda, dQ2dsidsj, dQ2dsidLambda, dQ2dLambda2)
    use tensor_math_mod, only : computeVMEquivalentStress, computeDVMEquivalentStressDstress, computeDVMEquivalentStressDstressiDStressj
    class(cp_base), intent(in) :: this
    real(k_real), intent(in) :: stress6(6), rss, lambda, schmid(5), f, q4
    real(k_real), intent(out) :: Q
    real(k_real), optional, intent(out) :: dQ2dsidsj(6,6), dQ2dsidLambda(6), dQ2dLambda2, dQdsi(6), dQdLambda
    real(k_real) :: svm, dsvm_dsi(6), dsvm_dsidsj(6,6)
    integer :: i,j
    __DECL_CLASS_UNUSED_THIS__

    __SUPPRESS_CLASS_UNUSED_THIS__
    svm = computeVMEquivalentStress(stress6)
    dsvm_dsi = computeDVMEquivalentStressDstress(stress6)

    associate (kVM => 2._k_real/3._k_real * f * q4)

    Q = (rss/lambda)**2 + kVM * (svm/(lambda*this%qVM))**2
    if (present(dQdsi)) then
      dQdsi(1:5) = 2._k_real*(rss/(lambda**2))*schmid(1:5) +  2._k_real* kVM * svm/(lambda*this%qVM)**2 *dsvm_dsi(1:5)
      dQdsi(6) = 0._k_real
    endif

    if (present(dQdLambda)) then
      dQdLambda = - 2._k_real*rss**2/(lambda**3) -  2._k_real * kVM *svm**2/(this%qVM**2 * lambda**3)
    endif

    if (present(dQ2dsidsj)) then
      dsvm_dsidsj = computeDVMEquivalentStressDstressiDStressj(stress6)
      do j=1,6
        do i=1,6
          if ((i.le.5.).and.(j.le.5)) then
            dQ2dsidsj(i,j) = 2._k_real/(lambda**2) * schmid(i)*schmid(j) + &
            2._k_real * kVM * 1._k_real/(lambda*this%qVM)**2 *dsvm_dsi(i)* dsvm_dsi(j) + &
            2._k_real * kVM * svm/(lambda*this%qVM)**2 * dsvm_dsidsj(i,j)
          else
            dQ2dsidsj(i,j) = 0._k_real
          endif
        enddo
      enddo
    endif

    if (present(dQ2dsidLambda)) then
      dQ2dsidLambda(1:5) = -4._k_real*rss/(lambda**3) * schmid(1:5) -4._k_real* kVM * svm/(lambda**3 *this%qVM**2) *dsvm_dsi(1:5)
      dQ2dsidLambda(6) = 0._k_real
    endif

    if (present(dQ2dLambda2)) then
      dQ2dLambda2 = 6._k_real*rss**2/(lambda**4) + 6._k_real * kVM *svm**2/(this%qVM**2 * lambda**4)
    endif
  end associate 

  end subroutine

  module subroutine computeM(this, lambda, s6, qM, M, dMdsi, dMdLambda, dM2dLambda2, dM2dsidLambda)
    class(cp_base), intent(in) :: this
    real(k_real), intent(in) :: lambda, s6, qM
    real(k_real), intent(out) :: M
    real(k_real), optional, intent(out) :: dMdsi(6), dMdLambda, dM2dLambda2, dM2dsidLambda(6)
    real(k_real), parameter :: s6_min=1e-10_k_real
    real(k_real) :: s6_use

    __DECL_CLASS_UNUSED_THIS__

    __SUPPRESS_CLASS_UNUSED_THIS__

    if (s6.eq.0._k_real) then
#ifdef __DEBUG__
       write(*,*) "voxel ", this%ix, this%iy, this%iz, "has 0. pressure, enforcing very small value "
#endif
       s6_use = s6_min
    else
        s6_use = s6
    endif

    M = qM*s6_use/(Lambda*sqrt(3._k_real))

      if (present(dMdsi)) then
        dMdsi(1:5) = 0; dMdsi(6) = qM/(Lambda*sqrt(3._k_real))
      endif

      if (present(dMdLambda)) then
        dMdLambda = - qM*s6_use/(Lambda**2*sqrt(3._k_real))
      endif

      if (present(dM2dLambda2)) then
        dM2dLambda2 = qM*2._k_real*s6_use/(sqrt(3._k_real)*lambda**3)
      endif

      if (present(dM2dsidLambda)) then
        dM2dsidLambda(1:5)=0._k_real; dM2dsidLambda(6) = - qM/(sqrt(3._k_real)*lambda**2)
      endif
  end subroutine

  module subroutine computeH(this, n, M, h, dhdM, dh2dM2)
    class(cp_base), intent(in) :: this
    real(k_real), intent(in) :: n, M
    real(k_real), intent(out) :: h
    real(k_real), optional, intent(out) :: dhdM, dh2dM2
    real(k_real) :: g, dgdM, dhdg, dh2dg2, dg2dM2
    __DECL_CLASS_UNUSED_THIS__

    __SUPPRESS_CLASS_UNUSED_THIS__
    g = 1._k_real/n *( 1.5_k_real*abs(M) )**( (n+1._k_real)/n )
    h = ( 1._k_real + g/this%qH )**(n)

    if (present(dhdM).or.present(dh2dM2)) then
      dgdM = ((n+1._k_real)/(n**2)) * ((1.5_k_real) * abs(M))**(1._k_real/n) * 1.5_k_real * sign(1._k_real, M)
      dhdg = n/this%qH * ( 1._k_real + g/this%qH )**(n-1._k_real)
    endif

    if (present(dhdM)) then
      dhdM = dhdg*dgdM
    endif

    if (present(dh2dM2)) then
      dh2dg2 = (n**2-n)/(this%qH**2)*(1._k_real+g/this%qH)**(n-2._k_real)
      dg2dM2 = 2.25_k_real*(n+1._k_real)/(n**3)*(1.5_k_real*abs(M))**(1._k_real/n-1._k_real) !2.25=9/4
      dh2dM2 = dh2dg2*dgdM**2 + dhdg*dg2dM2
    endif

  end subroutine

  module subroutine computeLeblondStressResidual(this, rss, lambda, schmid, stress6, n, f, R, dR_dLambda, dR_dsigma, dR2_dsidsj, dR2_dsidLambda, dR2_dLambda2)
    use print_utils_mod, only : printToScreen
    class(cp_base), intent(in) :: this
    real(k_real), intent(in)  :: rss, lambda, schmid(5), stress6(6), n, f
    real(k_real), intent(out)  :: R, dR_dLambda
    real(k_real), optional, intent(out)  :: dR_dsigma(6), dR2_dsidsj(6,6), dR2_dsidLambda(6), dR2_dLambda2
    real(k_real) :: Q, dQdLambda, dQdsi(6), dQ2dsidLambda(6), dQ2dlambda2, dQ2dsidsj(6,6), &
                    M, dMdLambda, dMdsi(6), dM2dLambda2, dm2dsidlambda(6), &
                    H, dh2dM2, dhDM
    integer :: i,j
    real(k_real) :: q1, q2, q3, q4
    q1 = 1._k_real
    q2 = n/(2.083_k_real*n-1._k_real)
    q3 = 1._k_real
    q4 = 1._k_real/(sqrt(3._k_real))


    call this%computeQ(stress6, rss, lambda, schmid, f, q4, Q, dQdLambda=dQdLambda)
    call this%computeM(lambda, stress6(6), q2, M, dMdLambda=dMdLambda)
    call this%computeH(n, M, h, dhDM=dhDM)
    associate( kn =>(n-1._k_real)/(n+1._k_real))

    R = Q + &
        q1*f*(h+kn/h) -1._k_real - q3*kn*f**2

    dR_dLambda = dQdLambda +  q1*f * dhDM*dMdLambda*(1._k_real-kn/h**2)
    if (present(dR_dsigma)) then
    call this%computeQ(stress6, rss, lambda, schmid, f, q4, Q, dQdsi=dQdsi)
    call this%computeM(lambda, stress6(6), q2, M, dMdsi=dMdsi)
    call this%computeH(n, M, h, dhDM=dhDM)
      dR_dsigma = dQdsi  +  q1*f * dhDM*dMdsi*(1._k_real-kn/h**2)
    endif

    if (present(dR2_dsidsj).or.present(dR2_dsidLambda).or.present(dR2_dLambda2)) then
      call this%computeQ(stress6, rss, lambda, schmid,  f, q4, Q, dQdsi=dQdsi, dQ2dsidsj=dQ2dsidsj, dQ2dsidLambda=dQ2dsidLambda, dQ2dLambda2=dQ2dLambda2)
      call this%computeM(lambda, stress6(6), q2, M, dMdsi=dMdsi, dM2dsidLambda=dM2dsidLambda, dM2dLambda2=dM2dLambda2)
      call this%computeH(n, M, h, dhDM=dhDM, dh2dM2=dh2dM2 )
    endif

    if (present(dR2_dsidsj)) then
      do j=1,6
        do i=1,6
          dR2_dsidsj(i,j) =dQ2dsidsj(i,j) +  q1*f * dMdsi(i)*dMdsi(j)*( &
          dh2dM2*(1._k_real-kn/h**2) + 2._k_real * dhdM**2*(kn/h**3) &
          )
      enddo; enddo
    endif

    if (present(dR2_dsidLambda)) then
        dR2_dsidLambda =dQ2dsidLambda + q1*f*( &
          dMdsi * dMdLambda * ( dh2dM2*(1._k_real-kn/h**2) + 2._k_real * dhdM**2*(kn/h**3) )+ &
          dhdM*dM2dsidLambda * (1._k_real-kn/h**2) &
          )
    endif

    if (present(dR2_dLambda2)) then
        dR2_dLambda2 = dQ2dLambda2 + q1*f*( &
          (dh2dM2*dMdLambda**2+dhdM*dM2dLambda2)*(1._k_real-kn/h**2) + &
          2._k_real*(dhdM*dMdLambda)**2*kn/h**3 &
          )
    endif

    end associate
  end subroutine

  module subroutine NRLeblondStressCPBase(this, stress6, schmid5, n, porosity, lambda)
    use print_utils_mod, only :printToScreen
    class(cp_base), intent(in) :: this
    real(k_real), intent(in) :: stress6(6), schmid5(5), n, porosity
    real(k_real), intent(out) :: lambda
    real(k_real) :: R, dlambda, dR_dLambda, NR_error, rss_dev
    real(k_real), parameter :: NR_tol = 1e-6_k_real
    integer, parameter :: max_NR_it = 400
    integer :: NR_iter

    rss_dev = sum(schmid5*stress6(1:5))
    
    lambda = rss_dev
    NR_error = 1._k_real
    NR_iter = 0
    do while((NR_error.ge.NR_tol).and.(NR_iter.le.max_NR_it))
      NR_iter = NR_iter+1
      ! compute residual
      call this%computeLeblondStressResidual(rss_dev, lambda, schmid5, stress6, n, porosity, R, dR_dLambda)
      NR_error = abs(R)
      if (NR_error>NR_tol) then ! prepare for next iteration
        ! compute newotn step
        dlambda = -R/dR_dLambda
        lambda = lambda + dlambda
      endif
    enddo
    if (NR_iter.gt.max_NR_it) error stop "NRLeblondStressCPBase maximum number of iteration reached. Abort"
    if (NR_error.gt.NR_tol) error stop "NRLeblondStressCPBase converged but NR_error>=NR_tol. Abort"
  end subroutine

  module subroutine LeblondComputeResolvedStressAndDirectionCPBase(this, stress6, schmid6, lambda, n, porosity, rss, flow_dir, dflowdir_dsigma )
    class(cp_base), intent(in) :: this
    real(k_real), intent(in) :: stress6(6), schmid6(6), lambda, n, porosity
    real(k_real), intent(out) :: rss, flow_dir(6), dflowdir_dsigma(6,6)

    rss = lambda

    call this%computeLambdaStressDerivatives(stress6, schmid6, lambda, n, porosity, &
                                             flow_dir, dflowdir_dsigma )

  end subroutine

  module subroutine computeLambdaStressDerivatives(this, stress6, schmid6, lambda, n, porosity, dLambda_dsi, dLambda2_dsidsj, rss_dev_out )
    class(cp_base), intent(in) :: this
    real(k_real), intent(in) :: stress6(6), schmid6(6), lambda, n, porosity
    real(k_real), intent(out) :: dLambda_dsi(6), dLambda2_dsidsj(6,6)
    real(k_real), optional, intent(out) :: rss_dev_out
    real(k_real) :: R, rss_dev, dR_dsigma(6), dR_dLambda, dR2_dLambda2, dR2_dsidLambda(6), dR2_dsidsj(6,6)
    integer :: i,j

    rss_dev = sum(stress6(1:5) * schmid6(1:5))

    call this%computeLeblondStressResidual(rss_dev, lambda, schmid6(1:5), stress6, n, porosity, R, dR_dLambda, &
                                          dR_dsigma=dR_dsigma, dR2_dsidsj=dR2_dsidsj, dR2_dsidLambda=dR2_dsidLambda, dR2_dLambda2=dR2_dLambda2)
    dLambda_dsi = - dR_dsigma /dR_dLambda ! this is dLambda/dSigma
    do j=1,6
      do i=1,6
    dLambda2_dsidsj(i,j) = -(dR2_dsidsj(i,j) + &
                             dR2_dsidLambda(i)*dLambda_dsi(j) + dR2_dsidLambda(j)*dLambda_dsi(i) + &
                             dR2_dLambda2 * dLambda_dsi(i)* dLambda_dsi(j))/dR_dLambda ! this is dLambda_dsi_dsj
    enddo; enddo

    if (present(rss_dev_out)) rss_dev_out = rss_dev
  end subroutine
end submodule
