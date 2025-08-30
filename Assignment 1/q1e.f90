program q1e
        implicit none
        integer :: i
        real :: rand_num(10),sum, average !rand_num(10) is an 1d array of 10 elements 

        !now generating 10 random numbers and writing it in test_ran.dat
        open(10, file = 'test_ran.dat', position = "append")
        call random_seed() !seed initializes a new query for random number generator
        do i=1,10
        call random_number(rand_num(i))
        end do

        !calculating the average
        sum = 0d0 !it indicates 0 as a double precision with 0 decimal points to the right of decimal
        do i=1,10
        sum = sum + rand_num(i) !we need the sum variable as it updates with each iteration of the do loop
        end do
        average = sum/10d0
        
        !writing random numbers to the opened dat file
        write(10,*) "Newly generated random numbers"
        do i=1,10
        write(10,*) rand_num(i)
        end do
        write(10,*) 'NOW calculating the average of 10 random numbers'
        write(10,*) 'Average=', average
        close(10)
        print*, 'Average=',average


end program q1e
