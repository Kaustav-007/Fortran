program questiond
        implicit none
        real :: rand_num(10)
        integer :: m
        integer :: seed

        open (unit =10, file="test_run.txt", action="write")
        
        
        call random_seed()
        
        do m=1,10
           call random_number(rand_num(m))                                                                                
        end do
        
        open (unit =10, file="test_run.txt", action="write")
        do m=1,10
             write (10,*) rand_num(m)
        end do

        write (10,*) 'Changing seed and generating 10 new random numbers' 
        close (10)
end program questiond
