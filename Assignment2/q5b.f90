program multi_dim_imp_mc

implicit none
integer::n,i,j
real*8::x(6),y,func,int_mc,var,sigma,p(2)
real*8::length,volume,sqrt2,gauss_dev

length=5.0d0
volume=acos(-1.0d0)**3
sqrt2=1.0d0/sqrt(2.0d0)

open(unit=7,file='q5b_impsamp.dat',status='unknown')

write(*,*) "Enter the value of n"
read(*,*) n

8	int_mc=0.0d0
	var=0.0d0
	sigma=0.0d0
	do i=1,n
	  do j=1,6
	    call random_number(p)
	    x(j)=gauss_dev(p)*sqrt2
	  end do
	  int_mc=int_mc+func(x)
	  sigma=sigma+func(x)*func(x)
	end do
	int_mc=int_mc/real(n)
	sigma=sigma/real(n)
	var=sigma-int_mc*int_mc


	int_mc=int_mc*volume
	sigma=volume*sqrt(var/real(n))

write(7,*) n,"",int_mc,"",sigma

	n=n*10
if (n.lt.1000000000) goto 8

end program multi_dim_imp_mc

 
real*8 function func(x)
implicit none
real*8::x(6),xy
real*8::a,b
a=0.5d0



xy=(x(1)-x(4))**2+(x(2)-x(5))**2+(x(3)-x(6))**2
func=exp(-a*xy)

end function
real*8 function gauss_dev(x)
implicit none
real*8::fact,sqr,p,x1,x2,x(2)

7   	call random_number(p)
	x1=2.0d0*p-1.0d0
	call random_number(p)
	x2=2.0d0*p-1.0d0

	sqr=x1*x1+x2*x2
if(sqr.ge.1.0d0.or.sqr.eq.0.0d0) goto 7

fact=sqrt(-2.0d0*log(sqr)/sqr)
gauss_dev=x2*fact
end function

