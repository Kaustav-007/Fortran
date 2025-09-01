program trapezoida_integration
implicit none

integer::i,n
real*8::a,b,h, func, x, fa, fb, trap_sum,error,pi
pi=22.0d0/7.0d0
open(3,file='trap_sin.dat',status='unknown')

write(*,*) "Give the value of n"
read(*,*) n

write(*,*) "give the lower limit of the integral"
read(*,*) a

write(*,*) "give the upper limit of the integral"
read(*,*) b

do
 n=n*5
 h=(b-a)/n

 fa=func(a)/2.0d0
 fb=func(b)/2.0d0

trap_sum=0.0d0

do i=1,n-1
   x=a+h*i
   trap_sum=trap_sum+func(x)
   end do
   
trap_sum=(trap_sum+fa+fb)*h
error=(pi-trap_sum)

write(3,*) h,"","",trap_sum,"","",error
if(n.ge.100000000)exit
end do

end program

real*8 function func(x)
implicit none
real*8::x
func=sin(x)
end function
	

