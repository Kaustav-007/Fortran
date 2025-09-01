program write_random_numbers
    implicit none
    integer :: i
    integer,parameter::m=10000
    double precision :: random_num(m)
    
    open(unit=5, file='ran_2b.dat', status='unknown')
    call random_number(random_num)
    


    do i = 1, m,2
        write(5, *) random_num(i),"",random_num(i+1)
    end do

    close(5)
end program write_random_numbers

