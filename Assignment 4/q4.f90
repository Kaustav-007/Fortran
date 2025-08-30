program RK4
    implicit none
    real*8 :: x0, y0, ya, dx,x1,f0,f1,f2,f3,x2,x3, y1, y2, y3,yf, diff_rk4
    integer :: i, niter
    ya=48.078
    open(30, file='RK4.dat')
    niter = 155
    ya=tan(1.550d0)
    x0=0
    y0=0
    dx = 1.550d0/niter
    print*, 'Step Size =', dx
    do i =1,niter
        f0=1+y0*y0

        x1=x0+(dx/2)
        y1=y0+(dx/2)*f0
        f1=1+y1*y1

        x2=x0+(dx/2)
        y2=y0+(dx/2)*f1
        f2=1+y2*y2
        
        x3=x0+dx
        y3=y0+dx*f2
        f3=1+y3*y3
        
        y0=y0+dx*(f0+2*f1+2*f2+f3)/6
        x0=x0+dx
        if(i==155) then
            yf=y0
        end if
        diff_rk4=ya-y0
    write(30,*)dfloat(i)*dx,y0, diff_rk4, abs(diff_rk4)/ya
    end do
    
    print*,y0
    y0=yf
    diff_rk4=ya-yf
    write(*,*)'The value of the difference y_A-y_rk4 at x=1.550 is', diff_rk4

    close(30)
    end program RK4
