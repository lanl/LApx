module green_function_mod
  use kinds
  use math_constants, only : PI
  use global, only : voxel_size
  implicit none
contains
  subroutine getGreenFunction(i,j,k, nx,ny,nz, xk, kx, ky, kz)
    integer, intent(in) :: i,j,k, nx,ny,nz
    integer, intent(out) :: kx,ky,kz
    real(k_real), intent(out) :: xk(3)

    if(i.le.nx/2) kx=i-1
    if(i.gt.nx/2) kx=i-nx-1

    if(j.le.ny/2) ky=j-1
    if(j.gt.ny/2) ky=j-ny-1

    if(k.le.nz/2) kz=k-1
    if(k.gt.nz/2) kz=k-nz-1

    xk(1)=2._k_real*PI*kx/(voxel_size(1)*nx)
    xk(2)=2._k_real*PI*ky/(voxel_size(2)*ny)
    xk(3)=2._k_real*PI*kz/(voxel_size(3)*nz)

  end subroutine

  subroutine getDDFT(igamma, i ,j, k, nx, ny, nz, xk, D_dft)
    integer, intent(in) :: igamma,i,j,k, nx,ny,nz
    real(k_real), intent(in) :: xk(3)
    real(k_real), intent(out):: D_dft(3,3)
    complex(k_real) :: xkmodc(3)
    complex(k_real), parameter :: ximag = (0.0,1.0)

    integer :: ii, jj
    select case (igamma)
    case (0)
      do ii=1,3
        do jj=1,3
          D_dft(ii,jj) = xk(ii)*xk(jj)
        end do
      end do

     !! AC implementing Miroslav's operator
   case (1)
      D_dft(1,1)=2.0* &
      (cos(2.0*PI*int2real(i-1)/int2real(nx))-1.0)
      D_dft(2,2)=2.0* &
      (cos(2.0*PI*int2real(j-1)/int2real(ny))-1.0)
      D_dft(3,3)=2.0* &
      (cos(2.0*PI*int2real(k-1)/int2real(nz))-1.0)
      D_dft(2,1)=1.0/2.0* &
      (cos(2.0*PI*int2real(i-1)/int2real(nx) + &
           2.0*PI*int2real(j-1)/int2real(ny))- &
       cos(2.0*PI*int2real(i-1)/int2real(nx) - &
           2.0*PI*int2real(j-1)/int2real(ny)))
      D_dft(3,1)=1.0/2.0* &
      (cos(2.0*PI*int2real(i-1)/int2real(nx) + &
           2.0*PI*int2real(k-1)/int2real(nz))- &
       cos(2.0*PI*int2real(i-1)/int2real(nx) - &
           2.0*PI*int2real(k-1)/int2real(nz)))
      D_dft(3,2)=1.0/2.0* &
      (cos(2.0*PI*int2real(j-1)/int2real(ny) + &
           2.0*PI*int2real(k-1)/int2real(nz))- &
       cos(2.0*PI*int2real(j-1)/int2real(ny) - &
           2.0*PI*int2real(k-1)/int2real(nz)))
      D_dft(1,2)=D_dft(2,1)
      D_dft(1,3)=D_dft(3,1)
      D_dft(2,3)=D_dft(3,2)

    case (2)
      do ii=1,3
        xkmodc(ii) = ximag*0.25*tan(xk(ii)/2.0)*&
       (1.0+exp(ximag*xk(1)))*(1.0+exp(ximag*xk(2)))*(1.0+exp(ximag*xk(3)))
      end do
      do ii=1,3
        do jj=1,3
          D_dft(ii,jj) = real(xkmodc(ii)*conjg(xkmodc(jj)), k_real) !-> get only the real part
        end do
      end do

    case default
      error stop "iGamma should be either 0, 1 or 2"

    end select

  end subroutine

end module
