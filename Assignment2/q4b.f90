program gaussian_random
    implicit none
    integer :: i, n
    real :: u1, u2, z1, sigma

    n = 10          ! Number of random numbers to generate
    sigma = 2.0     ! Standard deviation

    call random_seed()  ! Initialize random number generator

    print *, "Gaussian random numbers (mean = 0, SD = 2):"
    do i = 1, n
        call random_number(u1)  ! Generate two uniform random numbers U1 and U2
        call random_number(u2)
        
        ! Apply Box-Muller transform
        z1 = sqrt(-2.0 * log(u1)) * cos(2.0 * 3.141592653589793 * u2)
        z1 = z1 * sigma  ! Scale by the standard deviation
        
        print *, z1
    end do
end program gaussian_random

