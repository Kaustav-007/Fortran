program random_sum_distribution
  implicit none

  integer, parameter :: n = 100000 ! Number of random numbers per trial
  integer, parameter :: m = 100000 ! Number of trials
  real :: upper, lower, u_max, l_max, ranz, rand_num
  real, dimension(m) :: sums
  integer :: i, j, k, cnt

  open(2, file="dist_l.dat", status="unknown")

  ! Generate sums of random numbers for each trial
  do i = 1, m
     sums(i) = 0.0
     do j = 1, n
        call random_number(rand_num)
        if (rand_num < 0.5) then
            sums(i) = sums(i) - 1.0
        else
            sums(i) = sums(i) + 1.0
        end if
     end do
  end do

  ! Define histogram parameters
  u_max = 600.0
  l_max = -600.0
  ranz = 10.0
  cnt = int((u_max - l_max) / ranz)
  lower = l_max

  ! Calculate histogram and write results to file
  do k = 1, cnt
     upper = lower + ranz
     cnt = 0
     do i = 1, m
        if (sums(i) >= lower .and. sums(i) < upper) then
           cnt = cnt + 1
        end if
     end do
     write(2, '(F8.3, 1X, F8.5)') (lower + upper) / 2.0, real(cnt) / (m * ranz)
     lower = upper
  end do

  close(2)

end program random_sum_distribution

