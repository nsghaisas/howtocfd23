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
  real*8 :: pi, rx, ry
  character(500) :: fname, fname2

  
  !! Initialisation
  n1 = 1024                       !! No. of grids
  n2 = n1
  Np = 10                         !! No. of particles
  Nt = 5000                       !! Total no. of iteration
  dt = 1.0d-4
  pi = 4.0d0*datan(1.00d0)
  length = 2.0d0*pi               !! Box size
  dx = length/dble(n1)            !! smallest interval
  dy = dx
  

  allocate(uxreal(n1,n2),uyreal(n1,n2))
  allocate(xp(Np),yp(Np),vxp(Np),vyp(Np),uxp(Np),uyp(Np))
  xp=0.;yp=0.
  uxp=0.;uyp=0.
  vxp=0.;vyp=0

  call random_seed()
  !!Initialize particle position
!  call random_number(xp(:));
!  call random_number(yp(:));
 
 do i=1,Np
   call random_number(rx)
   xp(i) = rx
   call random_number(ry)
   yp(i) = ry
 end do
  
  
  write (*,*) 'Initialisation done'

  write (*,*) 'Reading real space data velocity data from datafile'

  open(unit = 12, file='output/x_Pos.out',status='unknown')
  open(unit = 13, file='output/y_Pos.out',status='unknown')
 
   do itime = 1,Nt

    write(fname,'(I4)')itime
    open(unit = 11, file='energy_analysis/vel_out/vel/ux'//trim(adjustl(fname))//'.dat',form='unformatted', status = 'old')
    read(11)((uxreal(i1,i2),i1=1,n1),i2=1,n2)
    close(11)

    open(unit = 11, file='energy_analysis/vel_out/vel/uy'//trim(adjustl(fname))//'.dat',form='unformatted', status = 'old')
    read(11)((uyreal(i1,i2),i1=1,n1),i2=1,n2)
    close(11)
! enddo

 !write (*,*)  'time evolution starts'
 ! do itime = 1,Nt
   do i = 1,Np
 !    call random_number(xp(i));
 !    call random_number(yp(i));
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

  enddo
    write(12,*)(xp(i), i=1,Np)
    write(13,*)(yp(i), i=1,Np)
enddo
  close(12)
  close(13)

 end program main
