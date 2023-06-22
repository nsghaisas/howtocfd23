program main
  use mod_part_interp
  use omplib

  implicit none
  integer :: i, i1, i2, n1, n2, Np, Nt, itime, intx, inty, intx_nxt, inty_nxt
  double precision :: dx, dy, length, dt
  real*8, allocatable, dimension(:) :: xp, yp, vxp, vyp, uxp, uyp
  real*8 :: x, y, xi, yj, xip1, yjp1, xip2, yjp2, xim1, yjm1
  real*8, allocatable, dimension(:,:) :: uxreal, uyreal
  real*8, allocatable, dimension(:,:) :: ux, uy
  real*8 :: pi
  character(500) :: fname, fname2

  
  !! Initialisation
  n1 = 128                        !! No. of grids
  n2 = n1
  Np = 10                         !! No. of particles
  Nt = 100                        !! Total no. of iteration
  dt = 1.0d-1
  pi = 4.0d0*datan(1.00d0)
  length = 2.0d0*pi               !! Box size
  dx = length/dble(n1)            !! smallest interval
  dy = dx
  

  allocate(uxreal(n1,n2),uyreal(n1,n2))
  allocate(xp(Np),yp(Np),Vxp(Np),Vyp(Np),uxp(Np),uyp(Np))
  xp=0.;yp=0.
  uxp=0.;uyp=0.
  vxp=0.;vyp=0

  call random_seed()
  !!Initialize particle position
  call random_number(xp(:));
  call random_number(yp(:));
  
  
  write (*,*) 'Initialisation done'

  write (*,*) 'Reading real space data velocity data from datafile'
  open(unit = 31, file = 'uxreal.dat', status ='unknown')
  open(unit = 32, file = 'uyreal.dat', status ='unknown')

  do i2 = 1,n2
     do i1 = 1,n1
        uxreal(i2,i1) = 0.5*(sin(dx*i1)*cos(dy*i2))
        uyreal(i2,i1) = 0.5*(cos(dx*i1)*sin(dy*i2))
      
     enddo
        write(31,*)(uxreal(i2,i1),i1=1,n1)
        write(32,*)(uyreal(i2,i1),i1=1,n1)
  enddo

  do itime = 1,Nt    
    do i = 1,Np
     x = modulo(xp(i)*length,length)
     y = modulo(yp(i)*length,length)

     intx = floor(x/dx)+1
     intx_nxt = (1-intx/n1)*intx+1
     inty = floor(y/dy)+1
     inty_nxt = (1-inty/n1)*inty+1

     !write(*,*)'x= ',x
     !write(*,*)'y= ',y
     !write(*,*)'intx= ',intx
     !write(*,*)'inty= ',inty
     !write(*,*)'intx_nxt= ',intx_nxt
     !write(*,*)'inty_nxt= ',inty_nxt

     xi = dfloat(intx-1)*dx; yj = dfloat(inty-1)*dy;
     xip1 = xi + dx; yjp1 = yj + dy;
     uxp(i) = (xip1 - x)*(yjp1 - y)*uxreal(intx,inty) + &
               (xip1 - x)*(y - yj)*uxreal(intx,inty_nxt) + &
               (x - xi)*(yjp1 - y)*uxreal(intx_nxt,inty) + &
               (x - xi)*(y - yj)*uxreal(intx_nxt,inty_nxt)
     uyp(i) = (xip1 - x)*(yjp1 - y)*uyreal(intx,inty) + &
               (xip1 - x)*(y - yj)*uyreal(intx,inty_nxt) + &
               (x - xi)*(yjp1 - y)*uyreal(intx_nxt,inty) + &
               (x - xi)*(y - yj)*uyreal(intx_nxt,inty_nxt)

     uxp(i) = uxp(i)/(dx*dy)
     uyp(i) = uyp(i)/(dx*dy)     

     vxp(i) = uxp(i)
     vyp(i) = uyp(i)

     xp = xp + vxp*dt
     yp = yp + vyp*dt

    open(unit = 12, file='output/x_Pos_for_controlled.out',status='unknown')
    write(12,*)(xp)
  
    open(unit = 13, file='output/y_Pos_for_controlled.out',status='unknown')
    write(13,*)(yp)
  enddo
enddo
  close(12)
  close(13)

 end program main
