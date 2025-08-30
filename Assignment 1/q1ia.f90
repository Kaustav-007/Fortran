program random_sum_distribution
  implicit none

  integer, parameter :: n = 10000  ! Number of random numbers per trial
  integer, parameter :: m = 10000  ! Number of trials
  real :: upper, lower, u_max, l_max, ranz, rand_num
  real :: sums(m)
  integer :: i, j, k, l, cnt

  open(1, file="dist_rand_walk_bin_0.5.dat", status="unknown")

  ! Initialize sums array to zero
  sums = 0.0

  ! Generate sums of random numbers for each trial
  do i = 1, m
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
  ranz = 0.5
  l = int((u_max - l_max) / ranz)

  lower = l_max

  ! Calculate histogram and write results to file
  do k = 1, l
     upper = lower + ranz
     cnt = 0
     do i = 1, m
        if (sums(i) >= lower .and. sums(i) < upper) then
           cnt = cnt + 1
        end if
     end do
     write(1, '(F8.3, 1X, I8, 1X, F8.5)') (lower + upper) / 2.0, cnt, real(cnt) / m
     lower = upper
  end do

  close(10)

end program random_sum_distribution
