program calculate_absolute_differences
  implicit none
  integer :: i, size
  real(8) :: average, difference
  real(8), dimension(:), allocatable :: random_numbers
  real(8), parameter :: reference_value = 0.50d0

  ! Array sizes to calculate averages for
  integer, dimension(3) :: sizes = [100, 10000, 1000000]

  ! Loop over each size
  do i = 1, 3
     size = sizes(i)
     
     ! Allocate the array for random numbers
     allocate(random_numbers(size))

     ! Generate random numbers between -1 and 1
     call random_seed()  ! Initialize the random number generator
     call random_number(random_numbers)
     random_numbers = 2.0 * random_numbers - 1.0  ! Transform to [-1, 1]

     ! Calculate the average
     average = sum(random_numbers) / size

     ! Calculate the absolute difference from 0.50d0
     difference = abs(reference_value - average)

     ! Print the result
     print *, 'Size:', size
     print *, 'Average:', average
     print *, 'Absolute difference from 0.50:', difference

     ! Deallocate the array
     deallocate(random_numbers)
  end do

end program calculate_absolute_differences

