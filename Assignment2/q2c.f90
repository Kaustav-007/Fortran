program correlation
implicit none
integer*8::i,j,n,k
real*8::p,avg,sigma,sum1,sq_sum,sq_avg,c_k
real*8,allocatable,dimension(:)::q

!avg=mean,sigma=st.dev,c_k=autocorrelation
sum1=0.0d0
sq_sum=0.0d0
open(unit=8,file='c_k.dat',status='unknown')
open(unit=9,file='scatter.dat',status='unknown')
write(*,*) "How many random numbers do you want to use?"
read(*,*)n

do i=1,n
  call random_number(p)
  sum1=sum1+p
  sq_sum=sq_sum+(p*p)
end do
avg=sum1/real(n)
sq_avg=sq_sum/real(n)
sigma=sqrt(sq_avg-(avg*avg))


write(*,*) "The mean is=",avg
write(*,*) "The standard deviation is=",sigma

!The section computes the autocorrelation function
allocate(q(n))
q=0.0d0
call random_number(q)
do j=1,3000
  sum1=0.0d0
  k=j-1
  do i=1,n-k
    sum1=sum1+(q(i)*q(i+k))
  end do
  sum1=sum1/real(n-k)
  c_k=(sum1-(avg*avg))/(sigma*sigma)
  write(8,*) k,c_k
end do
   do i=1,n-1
      write(9,*) q(i),q(i+1)
   end do
end program correlation
