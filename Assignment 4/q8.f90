program ring_rk4
  implicit none
  integer, parameter :: n=50, nsteps=2000
  real(8), parameter :: dt=0.02d0
  real(8) :: y(n), v(n), t
  real(8) :: k1y(n), k2y(n), k3y(n), k4y(n)
  real(8) :: k1v(n), k2v(n), k3v(n), k4v(n)
  integer :: i, ip, im, step

  ! Initialize positions and velocities
  y = 0.0d0
  v = 0.0d0
  y(1) = 0.8d0
  y(26) = 0.8d0

  do step = 1, nsteps
     t = step * dt

     ! Calculate k1
     do i=1,n
        ip = mod(i, n) + 1
        im = i - 1; if (im < 1) im = n
        k1y(i) = v(i)
        k1v(i) = y(ip) + y(im) - 2.0d0*y(i)
     end do

     ! Calculate k2
     do i=1,n
        ip = mod(i, n) + 1
        im = i - 1; if (im < 1) im = n
        k2y(i) = v(i) + 0.5d0*dt*k1v(i)
        k2v(i) = ( (y(ip) + 0.5d0*dt*k1y(ip)) + (y(im) + 0.5d0*dt*k1y(im)) - 2.0d0*(y(i) + 0.5d0*dt*k1y(i)) )
     end do

     ! Calculate k3
     do i=1,n
        ip = mod(i, n) + 1
        im = i - 1; if (im < 1) im = n
        k3y(i) = v(i) + 0.5d0*dt*k2v(i)
        k3v(i) = ( (y(ip) + 0.5d0*dt*k2y(ip)) + (y(im) + 0.5d0*dt*k2y(im)) - 2.0d0*(y(i) + 0.5d0*dt*k2y(i)) )
     end do

     ! Calculate k4
     do i=1,n
        ip = mod(i, n) + 1
        im = i - 1; if (im < 1) im = n
        k4y(i) = v(i) + dt*k3v(i)
        k4v(i) = ( (y(ip) + dt*k3y(ip)) + (y(im) + dt*k3y(im)) - 2.0d0*(y(i) + dt*k3y(i)) )
     end do

     ! Update y and v
     do i=1,n
        y(i) = y(i) + dt*(k1y(i) + 2.0d0*k2y(i) + 2.0d0*k3y(i) + k4y(i)) / 6.0d0
        v(i) = v(i) + dt*(k1v(i) + 2.0d0*k2v(i) + 2.0d0*k3v(i) + k4v(i)) / 6.0d0
     end do
  end do

  print *, 'Position of particle 1 after 2000 iterations:', y(1)

end program ring_rk4

