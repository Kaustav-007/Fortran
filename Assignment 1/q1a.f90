program print_random_numbers
    implicit none
    integer :: i
    real*8 :: rand_num

    ! Seed the random number generator
    call random_seed()

    ! Print 10 random numbers between 0 and 1
    do i = 1, 10
        call random_number(rand_num)
        print *, rand_num
    end do

end program print_random_numbers

