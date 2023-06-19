program berger_2d
  implicit none

  integer, parameter :: rp = selected_real_kind(8)
  integer :: i, j, itr, n, itr_max
  real(rp) :: dx, dy, dt, x, maxT, nu, lambda
  real(rp), allocatable, dimension(:,:) :: u, unew, v, vnew, up, vp
  
  ! initialize parameters
  dx = 0.01_rp
  dy = dx
  dt = 0.001_rp
  maxT = 1.0_rp
  n = 1000
  itr_max = floor(maxT / dt)
  nu = 0.01_rp
  lambda = dt / dx

  allocate(u(n+1,n+1),v(n+1,n+1),unew(n+1,n+1),vnew(n+1,n+1),up(n+1,n+1),vp(n+1,n+1))
 
  ! initialize arrays to zero
  u = 0.0_rp
  unew = 0.0_rp
  v = 0.0_rp
  vnew = 0.0_rp
  up = 0.0_rp
  vp = 0.0_rp
 
  ! intial condition
  do i = 100, 401
    do j = 1, n+1
      u(i,j) = 1.0_rp
      v(i,j) = 1.0_rp
    end do
  end do

  open(unit=1, file='heaviside_2d_f90.dat')
  do i = 1, n+1
    do j = 1, n+1
      write(1,*) (i-1)*dx, (j-1)*dy, u(i,j), v(i,j)
    end do
  end do
  close(1)
  
  ! scheme
  do itr = 1, itr_max
    do i = 2, n
      do j = 2, n

         !..........

      end do
    end do
    do i = 2, n
      do j = 2, n
          !..........
      end do
    end do
    do i = 2, n
      do j = 2, n
          !............
      end do
    end do
  end do

  ! writing output
  open(unit=2, file='burger_2d_f90.dat', form='formatted')
  do i = 1, n+1
    do j = 1, n+1
      write(2,*) (i-1)*dx, (j-1)*dy, u(i,j), v(i,j), sqrt(u(i,j)**2 + v(i,j)**2)
    end do
  end do
  close(2)

  deallocate(u,v,unew,vnew,up,vp)
end program 
