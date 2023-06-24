program main
use nvtx
use openacc
!!use derivative_mod_gpu


implicit none

integer::n2,n1,itime,Ttime
integer::i,j
integer::isend1,isend2,irecv1,irecv2
real*8:: Lx,Ly,pi,dx,dy,Diff,dt
real*8::tini,tfin
real*8,allocatable,dimension(:,:)::T,Trhs
character(100)::fname,ftime
integer,parameter::istr=0

pi=4.d0*datan(1.d0)
Lx=2.e0*pi;
Ly=Lx;
Ttime=100000;


n1=4096;n2=4096
Diff=1.e-2
dt=5.e-5;

dx=Lx/n1;dy=Ly/n2;



allocate(T(0:n1+1,0:n2+1),Trhs(n1,n2))


!! Initialize
T=0;Trhs=0;

!! Initialize the field
open(unit=10,file="init_config.dat",status="unknown")
do j = 1, n2
   do i = 1, n1
      T(i,j) = sin(2.0*pi*i*dx/Lx)*sin(2.0*pi*j*dy/Ly)
      write(10,*) i,j, T(i,j)
   enddo
   write(10,*) 
enddo
close(10)


!!call cpu_time(tini)
!$acc data copy(T,Trhs,n1,n2,dx,dy,Diff,dt) async(istr)

call cpu_time(tini)
do itime=1,Ttime

   call nvtxStartRange("SET PBC")

   !$acc parallel loop private(i) present(T) async(istr)
   do i=1,n1
    T(i,0)=T(i,n2)
    T(i,n2+1)=T(i,1)
   enddo 
   !$acc end parallel

   !$acc parallel loop private(j) present(T) async(istr)
   do j=0,n2+1   
   T(0,j)=T(n1,j)
   T(n1+1,j)=T(1,j)
   enddo
   !$acc end parallel

   call nvtxEndRange


   call nvtxStartRange("Derivative")
   !! Evaluate RHS
   !$acc parallel loop collapse(2) private(i,j) present(T,Trhs,dx,dy,Diff) async(istr) 
   DO j = 1,n2
      DO i = 1, n1
         Trhs(i,j) = Diff*( (T(i+1,j) -2.d0*T(i,j) + T(i-1,j))/(dx*dx)  +&
              ( T(i,j+1) -2.d0*T(i,j) + T(i,j-1))/(dy*dy));
      ENDDO
   ENDDO
   !$acc end parallel   
   call nvtxEndRange

   
   !! Update using Euler scheme
   call nvtxStartRange("Update")
!!   !$acc parallel loop present(T,Trhs,dt) tile(128,8) 
   !$acc parallel loop gang collapse(2)  present(T,Trhs,dt) async(istr)
!!   !$acc kernels present(T,Trhs,dt)
   do j=1,n2
      do i=1,n1
        T(i,j)=T(i,j) +  dt*Trhs(i,j);
        enddo
    enddo
!!   !$acc end kernels
   !$acc end parallel 
   call nvtxEndRange
enddo
!$acc end data
!$acc wait(istr)
call cpu_time(tfin);
write(*,*) itime

write(ftime,'(i8)') itime
open(unit=10,file="diff_sp_"//trim(adjustl(ftime))//".dat",status="unknown")
i=n1/4
do j = 1, n2
!!   do i = 1, n1
      write(10,*) i,j, T(i,j),Trhs(i,j)
!!   enddo
!!   write(10,*) 
enddo
close(10)


write(*,*) "Time elaspsed in seconds:",tfin-tini


END program main
