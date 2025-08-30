program euler
implicit none
real(8)::dx,x,y,y1,y2,fxy,f0,f1,fe,niter,ya,ye,yme,yie,diff_e,diff_me,diff_ie
integer::i,n

ya=48.078

!Euler method!
open(30, file='euler1.dat')
print*,'Enter the step size'
read*,dx
x=0
y=0
niter=1.55/dx
do i =1,int(niter)+1
fxy=1+y*y
y=y+fxy*dx
write(30,*)dfloat(i)*dx,y
end do
ye=y
print*,ye
diff_e=ya-ye
write(*,*)'The value of the difference y_A-y_E at x=1.550 is', diff_e
close(30)


!Modified Euler!

open(31, file='euler_mod1.dat')
x=0
y=0
do i=1, int(niter)+1
f0 = 1 + y*y
y1 = y + (f0)*dx/2
f1 = 1 + y1*y1
y2 = y + f1*dx
write(31,*)dfloat(i)*dx, y2
y = y2
end do
yme = y
print*, yme
diff_me = ya - yme

write(*,*)'The value the difference y_A - y_ME at x=1.550 is', diff_me

close(31)
!Improved Euler
open(32, file='euler_mod2.dat')
x=0
y=0
do i=1,int(niter)+1
f0=1+y*y
y1=y+(f0)*dx
fe=1+y1*y1
fxy = (f0 + fe)/2
y = y + fxy*dx
write(32,*)dfloat(i)*dx, y
end do

yie = y
print*, yie
diff_ie = ya - yie
write(*,*)'The value the difference y_A - y_IE at x=1.550 is', diff_ie
close(32)

end program euler


