program sin_fft
use,intrinsic :: iso_c_binding
implicit none

include "fftw3.f03"


integer*8:: i1,i2,pfor,pinv
integer:: nx,ny,nn,kx,ky,ksqr,ireal,iimag
real*8::lx,ly,dx,dy,pi,xx,yy
real*8::rk2,scale

  integer,allocatable,dimension(:):: ndim
  real*8,allocatable,dimension(:,:):: omg,psi

nn=128
nx = nn;
ny = nx;

  allocate(ndim(2),omg(nx+2,ny),psi(nx+2,ny))
ndim(1)=nx
ndim(2)=ny
pi = 4.0d0*datan(1.0d0);
lx = 2.0d0*pi;  
ly = lx;
dx = lx/dble(nx);  
dy = ly/dble(ny);
scale = 1.0d0/(dble(nx)*dble(ny));
psi = 0.0d0;
omg = 0.0d0;


!! %-------create plan for forward transform -------------------------
        call dfftw_plan_dft_r2c_2d(pfor,nx, ny, omg, omg,FFTW_MEASURE);
!! %-------create plan for inverse transform -------------------------
        call dfftw_plan_dft_c2r_2d(pinv,nx, ny, psi, psi,FFTW_MEASURE);
    print*,'FFTW Plan Created'


open(unit=36,file='omg.out',status='unknown');
open(unit=37,file='psi.out',status='unknown');
!! %---
  do i2=1,ny
      yy=dble(i2)*dy
    do i1=1,nx
      xx=dble(i1)*dx;
      omg(i1,i2)=sin(5.0d0*xx)
    end do
  end do


do i1=1,nx
write(36,*) (omg(i1,i2),i2=1,ny)
end do
!! Do the forward Fourier transform
   call dfftw_execute(pfor)
    print*,'FFTW forward transform done'

  do i2 = 1,ny
    ky = (i2-1)-(i2/(ny/2+2))*ny
    do i1 = 1,nx/2+1
      kx = i1-1
      ksqr = kx*kx+ky*ky
      ireal = 2*i1-1
      iimag = 2*i1
      rk2 = dble(ksqr)
      if(ksqr==0.or.kx==nx/2)then
      psi(ireal,i2) = 0.0d0;
      psi(iimag,i2) = 0.0d0;
      else
      psi(ireal,i2) = -omg(ireal,i2)/rk2;
      psi(iimag,i2) = -omg(iimag,i2)/rk2;
      endif
    enddo
  enddo
!! Do the inverse Fourier transform
   call dfftw_execute(pinv)
  psi = psi*scale
  psi(nx+1:nx+2,1:ny) = 0.0d0

do i1=1,nx
write(37,*) (psi(i1,i2),i2=1,ny)
end do
close(36)
close(37)

call dfftw_destroy_plan(pfor,pinv)

!print*,ipsi

end program sin_fft
