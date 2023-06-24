program main
use mpi
implicit none

integer::n2,n1,local_n2,local_n1,itime,Ttime
integer::ierr,nprocs,myrank,inext,iprev,i,j,jloc
integer::isend1,isend2,irecv1,irecv2
real*8:: Lx,Ly,pi,dx,dy,Diff,dt
real*8::tini,tfin
real*8,allocatable,dimension(:,:)::T,Trhs
character(100)::fname,ftime

INTEGER istatus(MPI_STATUS_SIZE)
CALL MPI_INIT(ierr)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr) 
CALL MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr) 

pi=4.d0*datan(1.d0)
Lx=2.e0*pi;
Ly=Lx;
Ttime=100000;


n1=2048;n2=2048
Diff=1.e-2
dt=1.e-4;

dx=Lx/n1;dy=Ly/n2;

!! Only slab decomposition  along the Y-direction
local_n2=n2/nprocs


allocate(T(0:n1+1,0:local_n2+1),Trhs(0:n1+1,0:local_n2+1))

!! Initialize
T=0;Trhs=0;

inext = myrank + 1
iprev = myrank - 1
!!IF (myrank == nprocs - 1) inext = MPI_PROC_NULL
!!IF (myrank == 0) iprev = MPI_PROC_NULL

if (myrank == nprocs - 1) inext = 0
if (myrank == 0) iprev = nprocs-1

!! Initialize the field
write(fname,'(i8)') myrank+1
open(unit=10,file="init_config_"//trim(adjustl(fname))//".dat",status="unknown")
do jloc = 1, local_n2
   j = jloc + myrank*local_n2
   do i = 1, n1
      T(i,jloc) = sin(2.0*pi*i*dx/Lx)*sin(2.0*pi*j*dy/Ly)
      write(10,*) i,j, T(i,jloc)
   enddo
   write(10,*) 
enddo
close(10)

tini=mpi_wtime();

do itime=1,Ttime

   !! Set PBC
   CALL MPI_ISEND(T(1,local_n2),n1,MPI_REAL8,inext,1,MPI_COMM_WORLD,isend1,ierr) 
   CALL MPI_ISEND(T(1,1),n1,MPI_REAL8,iprev,1,MPI_COMM_WORLD,isend2,ierr) 
   CALL MPI_IRECV(T(1,0),n1,MPI_REAL8,iprev,1,MPI_COMM_WORLD,irecv1,ierr) 
   CALL MPI_IRECV(T(1,local_n2+1),n1,MPI_REAL8,inext,1,MPI_COMM_WORLD,irecv2,ierr) 
   
   CALL MPI_WAIT(isend1, istatus, ierr)
   CALL MPI_WAIT(isend2, istatus, ierr)
   CALL MPI_WAIT(irecv1, istatus, ierr)
   CALL MPI_WAIT(irecv2, istatus, ierr)
   
   T(0,0:local_n2+1)=T(n1,0:local_n2+1)
   T(n1+1,0:local_n2+1)=T(1,0:local_n2+1)


   !! Evaluate RHS
   DO j = 1,local_n2
      DO i = 1, n1
         Trhs(i,j) = Diff*( (T(i+1,j) -2.d0*T(i,j) + T(i-1,j))/(dx*dx)  +&
              (T(i,j-1) + T(i,j+1) -2.d0*T(i,j))/(dy*dy));
      ENDDO
   ENDDO
   
   !! Update using Euler scheme
   T=T + dt*Trhs
enddo

tfin=mpi_wtime();



write(fname,'(i8)') myrank+1
write(ftime,'(i8)') itime
open(unit=10,file="diff_"//trim(adjustl(ftime))//"_"//trim(adjustl(fname))//".dat",status="unknown")
do jloc = 1, local_n2
   j = jloc + myrank*local_n2
   do i = 1, n1
      write(10,*) i,j, T(i,jloc)
   enddo
   write(10,*) 
enddo
close(10)


if(myrank==0) write(*,*) "Time elaspsed in seconds:",tfin-tini


CALL MPI_FINALIZE(ierr) 
END program main
