
! Fortran program to generate random numbers with exponential and Gaussian distributions
program random_distributions
  implicit none

  integer, parameter :: n = 10000000     ! Number of random numbers to generate
  real :: uniform_rand, exp_rand, gauss_rand
  real, dimension(n) :: exp_numbers, gauss_numbers
  integer :: i
  open(10, file='exponential_data.dat', status='replace')
  open(20, file='gaussian_data.dat', status='replace')

  ! Seed the random number generator
  call random_seed()

  ! Generate exponential random numbers (f(x) = 2 * exp(-2x))
  do i = 1, n
    call random_number(uniform_rand)            ! Generate a uniform random number in [0, 1)
    exp_numbers(i) = -0.5 * log(1.0 - uniform_rand)  ! Transform to exponential distribution
    write(10, *) exp_numbers(i)
  end do

  ! Generate Gaussian random numbers (mean = 0, standard deviation = 2)
  do i = 1, n, 2
    call random_number(uniform_rand)            ! First uniform random number
    call random_number(gauss_rand)              ! Second uniform random number
    
    ! Box-Muller transformation
    exp_rand = sqrt(-2.0 * log(uniform_rand)) * cos(2.0 * 3.14159 * gauss_rand)
    gauss_rand = sqrt(-2.0 * log(uniform_rand)) * sin(2.0 * 3.14159 * gauss_rand)
    
    gauss_numbers(i) = 2.0 * exp_rand           ! Scale by standard deviation (SD = 2)
    if (i < n) gauss_numbers(i+1) = 2.0 * gauss_rand
    write(20, *) gauss_numbers(i)
    if (i < n) write(20, *) gauss_numbers(i+1)
  end do

  close(10)
  close(20)

  print *, "Data written to 'exponential_data.dat' and 'gaussian_data.dat'."

end program random_distributions


