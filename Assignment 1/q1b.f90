program write_random_numbers
    implicit none
    integer :: i
    real :: random_num

    open(unit=10, file='test_ran.dat', status='replace')
    
    call random_seed()  ! Initialize the random number generator

    do i = 1, 10
        call random_number(random_num)
        write(10, '(F6.4)') random_num
    end do

    close(10)
end program write_random_numbers

