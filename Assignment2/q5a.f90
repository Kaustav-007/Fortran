program brute_force_integration
!6-d integral
integer::n,i,j
real*8::x(6),y,func,int_mc,var,sigma,p(6)
real*8::length,volume

length=5.0d0
volume=(2.0d0*length)**6

open(unit=1,file='q5a_brute.dat')
write(*,*)"Enter the value of n"
read(*,*) n

int_mc=0.0d0
var=0.0d0
sigma=0.0d0

7	do i=1,n
	  call random_number(p)
	  do j=1,6
	    x(j)=-length+2.0d0*length*p(j)
	  end do
	  int_mc=int_mc+func(x)
	  sigma=sigma+func(x)*func(x)
	end do

	int_mc=int_mc/real(n)
	sigma=sigma/real(n)
	var=sigma-int_mc*int_mc
	int_mc=volume*int_mc

	sigma=volume*sqrt(var/real(n))
	write(1,*) n,"",int_mc,"",sigma

n=n*10
if(n.lt.100000000) goto 7
 
 end program brute_force_integration
 
real*8 function func(x)
implicit none
real*8::x(6),xx,yy,xy
real*8::a,b
a=1.0d0
b=0.5d0

xx=x(1)*x(1)+x(2)*x(2)+x(3)*x(3)
yy=x(4)*x(4)+x(5)*x(5)+x(6)*x(6)
xy=(x(1)-x(4))**2+(x(2)-x(5))**2+(x(3)-x(6))**2
func=exp(-a*xx-a*yy-b*xy)

end function
