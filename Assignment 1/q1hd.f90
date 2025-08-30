program random_sum_distribution
  implicit none
  integer, parameter :: n = 10000
  integer, parameter :: m = 10000
  real :: random_value, upper, lower, u_max, l_max, ranz,r_1,a,b,r
  real :: sums(m)
  integer :: i, j, k, l, cnt
  a=-1.0
  b=1.0
  open(3,file="dist_range(-1 TO +1)",status="unknown")
  ! Generate sums of 10000 random numbers between 0 and 1 for 10000 trials
  do i = 1, m
     sums(i) = 0.0d0
     do j = 1, n
        call random_number(r)
        r_1=a+(b-a)*r
        sums(i) = sums(i) + r_1
     enddo
  enddo

 u_max = 250 ;  l_max =-250
 ranz = 0.5; l = (u_max-l_max)/ranz


lower  = l_max
 do k = 1, l
        upper = lower+ranz
        cnt = 0
        do i = 1, m
                if(sums(i)>=lower .and. sums(i)<upper) cnt = cnt + 1
        enddo
        write(3,*) (lower+upper)/2, cnt, real(cnt)/m
        lower = upper
 enddo
  close(10)

  END
