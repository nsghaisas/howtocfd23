program 1d_berger
  implicit none
  
  integer, parameter :: rp = selected_real_kind(8)
  integer :: i, itr, n, scheme, itr_max
  real(rp) :: dx, dt, x, lambda, maxT, nu
  real(rp), dimension(:), allocatable :: f, fnew
  real(rp), dimension(:), allocatable :: fp
  character(len=10) :: scheme_name 
  
  ! initialize parameters
  dx = 0.01_rp
  dt = 0.000001_rp
  nu = 0.01_rp
  maxT = 1.0_rp
  n = 1000
  itr_max = floor(maxT / dt)
  lambda = dt / dx

  ! allocate vectors
  allocate(f(n+1), fnew(n+1))
  allocate(fp(n+1)) ! for mac-cormack scheme

  ! initialize vectors to zero
  f = 0.0_rp
  fnew = 0.0_rp
  fp = 0.0_rp


  do i = 1, n+1
    if (i > 100 .and. i < 300) then
      f(i) = 1.0_rp
    else
      f(i) = 0.0_rp
    end if
  end do

  !write(*,*) "Choose your favorite scheme"
  !write(*,*) "1=FTCS", "; 2=FTBS", "; 3=Lax", "; 4=Lax-Wendroff"
  !write(*,*) "5=Beam-warming", "; 6=Mac-Cormak"
  !read(*,*) scheme

  !if (scheme == 1) then ! FTCS
  !  do itr = 1, itr_max
  !    do i = 2, n
  !      fnew(i) = f(i) - 0.5 * lambda * f(i) * (f(i+1) - f(i-1))
  !    end do
  !    do i = 2, n
  !      f(i) = fnew(i)
  !    end do
  !    f(n+1) = f(2)
  !    f(1) = f(n)
  !  end do
  !end if
  !if (scheme == 2) then ! FTBS
  !  do itr = 1, itr_max
  !    do i = 2, n
  !      fnew(i) = f(i) - lambda * f(i) * (f(i) - f(i-1))
  !    end do
  !    do i = 2, n
  !      f(i) = fnew(i)
  !    end do
  !    f(2) = f(n+1)
  !    f(1) = f(n)
  !  end do
  !end if
  !if (scheme == 3) then ! Lax-Friedrichs
  !  do itr = 1, itr_max
  !    do i = 2, n
  !      fnew(i) = 0.5 * (f(i+1) + f(i-1) - lambda * f(i) * (f(i+1) - f(i-1)))
  !    end do
  !    do i = 2, n
  !      f(i) = fnew(i)
  !    end do
  !    f(n+1) = f(2)
  !    f(1) = f(n)
  !  end do
  !end if
  !if (scheme == 4) then ! Lax-Wendroff
  !  do itr = 1, itr_max
  !    do i = 2, n
  !      fnew(i) = f(i) - 0.5 * lambda * f(i) * (f(i+1) - f(i-1)) + &
  !                0.5 * lambda * lambda * f(i) * f(i) * (f(i+1) - 2*f(i) + f(i-1))
  !    end do
  !    do i = 2, n
  !      f(i) = fnew(i)
  !    end do
  !    f(n+1) = f(2)
  !    f(1) = f(n)
  !  end do
  !end if
  !if (scheme == 5) then ! Beam-Warming
  !  do itr = 1, itr_max
  !    do i = 2, n
  !      fnew(i) = (1 - 1.5 * lambda * f(i) + 0.5 * lambda * lambda * f(i) * f(i)) * f(i) + &
  !                lambda * f(i) * (2 - lambda * f(i) ) * f(i-1) + &
  !                0.5 * lambda * f(i) * (lambda * f(i) - 1) * f(i-2)
  !    end do
  !    do i = 2, n
  !      f(i) = fnew(i)
  !    end do
  !    f(n+1) = f(2)
  !    f(1) = f(n)
  !  end do
  !end if
  
  ! Mac-Cormak
    do itr = 1, itr_max
      !Write Mac-Cormak scheme here
    end do
  
  write(scheme_name, '(I10)') scheme
  open(unit=20, file="burger_1d_f90.dat", form="formatted")
  do i = 1, n+1
    write(20,*) (i-1)*dx, f(i)
  end do
  close(20)
  deallocate(f, fnew, fp)
end 
