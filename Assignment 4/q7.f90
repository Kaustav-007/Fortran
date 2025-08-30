program SOSHO
implicit none

real(8)::x0,v0,h,f0,x1,v1,f1,f1v,x2,v2,f2,f2v,x3,v3,f3,f3v,E,f0v
integer::i,niter
x0=0
v0=2.00001
niter=5000
h=0.01
open(39,file='q7_SHO.dat')
do i =1,niter
f0=v0;f0v=-sin(x0)
x1=x0+f0*(h/2);v1=v0+f0v*(h/2);f1=v1;f1v=-sin(x1)
x2=x0+f1*(h/2);v2=v0+f1v*(h/2);f2=v2;f2v=-sin(x2)
x3=x0+f2*h;v3=v0+f2v*h;f3=v3;f3v=-sin(x3)


x0=x0+(h/6)*(f0+2*f1+2*f2+f3)
v0=v0+(h/6)*(f0v+2*f1v+2*f2v+f3v)

E=v0*(v0/2)-cos(x0)
write(39,*)h*dfloat(i),x0,v0,E
end do
print*,x0
end program
