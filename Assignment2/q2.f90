program write_random_numbers
    implicit none
    integer :: i,m
    double precision :: random_num
    m=100000

    open(unit=10, file='ran_uni.dat', status='unknown')
    


    do i = 1, m
        call random_number(random_num)
        write(10, *) random_num
    end do

    close(10)
end program write_random_numbers

