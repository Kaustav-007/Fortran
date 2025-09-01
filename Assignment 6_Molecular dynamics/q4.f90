module parameter
    implicit none
    integer :: zzz=4280145, zzzz=77777, i, j, k, c ,p,divide
    integer, parameter :: npart = 1200, niter=300000
    integer, parameter :: lx=20, ly=20, lz=20
    real*8, parameter :: mass=1.0d0, kb_T=1.0d0, rc=2.50d0, rs=4.50d0, sigma=1.0d0, eps=4.0d0, delta_t=0.0025d0, pi=acos(-1.0d0)
    real*8, parameter :: sigma6=sigma**6, sigma12=sigma**12
    real*8, parameter :: fc=eps*((12.0d0*sigma12/(rc**13))-(6.0d0*sigma6/(rc**7)))
    real*8, parameter :: ufc=fc*rc + eps*(((sigma/rc)**12)-((sigma/rc)**6))
    real*8 :: avg_vx, avg_vy, avg_vz, PE, KE, gx, gy, gz
    real*8 :: x, y, z, x1, y1, z1, x2, y2, z2, dx, dy, dz,fs, r, lj, lj_force
    real*8 :: invr, ir2, ir6, theoryKE, scalef
    real*8, parameter :: llx=dfloat(lx), lly=dfloat(ly), llz=dfloat(lz)
    real*8, parameter :: llxby2=llx/2.0d0, llyby2=lly/2.0d0, llzby2=llz/2.0d0
    real*8, dimension(3*npart) :: pos, vel, force, old_force, acc, old_acc
    integer ,allocatable :: no_neighb(:) , neigh_list (:,:)
    integer :: bin, nbins
    real(8) ::  dr , rmax, vol, rho, shell_vol
    real(8), allocatable :: gr(:)
 end module parameter
 
 !==========MAIN PROGRAM==========
 program md
    use parameter
    implicit none
    allocate (no_neighb(npart))
    allocate (neigh_list(npart,npart))
    allocate(gr(nbins))
    !integer :: time
    
    call position_init
    call velocity_init
    call neighbour_list
    call calc_force
    
    open(14,file='force.dat',status='unknown',form='formatted')
    do i=1,npart
       write(14,*) force(3*i-2), force(3*i-1), force(3*i)
    end do
    close(14)  

    open(30,file='grerror.dat',status='unknown',form='formatted')
    open(10,file='neghbourlist.dat',status='unknown',form='formatted')
    open(11,file='PE+KE+TE(WITH THERMOSTAT ).dat',status='unknown',form='formatted')
    open(20, file='q5_paircorrelation.dat', status='unknown', form='formatted')
    
    divide=0
    gr = 0.0d0
    do k=1,niter

       call update_pos
       if(mod(k,40)==0) then
        call neighbour_list
        end if
       call calc_force
       call update_vel
       call thermostat
       
       if (k>50000) then
         if(mod(k,100)==0) then
            divide=divide+1
         
         call compute_gr
         end if 
      end if
       

      
     if(mod(k,1000)==0) then
      !  write(11,*) dfloat(k),PE/dfloat(npart),kE/dfloat(npart),(KE+PE)/dfloat(npart)
      write(*,*) k
      end if
       
    end do
    do bin = 1, nbins
      r = (bin - 0.5d0) * dr
      write(20,*) r, gr(bin)/dfloat(divide)
  end do
  close(20)
deallocate(gr)
    close(11)
           
 end program md
 
 !==========INITIATION OF POSITION==========  
 
 subroutine position_init
    use parameter
    use omp_lib
    implicit none
    integer :: n
    integer,allocatable :: seed(:)
 
    call random_seed(size=n)
    allocate(seed(n))
    seed = zzz  
    call random_seed(put=seed)
 
    do i=1,npart
  62    call random_number(x)
       x= 0.5d0 + (x*(lx-1))
       call random_number(y)
       y= 0.5d0 + (y*(ly-1))
       call random_number(z)
       z= 0.5d0 + (z*(lz-1))
       do j= 1, i-1
          if (sqrt((x - pos(3*j-2))**2 + (y - pos(3*j-1))**2 + (z - pos(3*j))**2) < sigma) then
             goto 62
          end if
       end do
      
       pos(3*i-2) = x
       pos(3*i-1) = y
       pos(3*i) = z
    end do      
    open(1,file='initial_positions.dat',status='unknown')
    do i=1,npart
       write(1,*) pos(3*i-2), pos(3*i-1), pos(3*i)  
    end do
    close(1)
    deallocate(seed)
 end subroutine position_init
 
 !==========INITIATION OF VELOCITY==========
 
 subroutine velocity_init
    use parameter
    use omp_lib
    implicit none
    real*8 :: vel_const, vx, vy, vz
    integer :: n
    integer,allocatable :: seed(:)
 
    call random_seed(size=n)
    allocate(seed(n))
    seed = zzzz  
    call random_seed(put=seed)
    
    vel_const = dsqrt(12.0d0*kb_T/mass)
    
    do i= 1, npart
       call random_number(vx)
       call random_number(vy)
       call random_number(vz)
       vel(3*i-2) = vel_const*(vx-0.5d0)
       vel(3*i-1) = vel_const*(vy-0.5d0)
       vel(3*i) = vel_const*(vz-0.5d0)
    end do  
    
    avg_vx = 0.0d0; avg_vy = 0.0d0; avg_vz = 0.0d0
    do i = 1, npart
       avg_vx = vel(3*i-2)+avg_vx; avg_vy = vel(3*i-1)+avg_vy; avg_vz = vel(3*i)+avg_vz
    end do
    avg_vx = avg_vx/dfloat(npart); avg_vy = avg_vy/dfloat(npart); avg_vz = avg_vz/dfloat(npart)
    
    do i= 1, npart
       vel(3*i-2) = avg_vx - vel(3*i-2)
       vel(3*i-1) = avg_vy - vel(3*i-1)
       vel(3*i) = avg_vz - vel(3*i)
       KE = KE + (vel(3*i-2)*vel(3*i-2)+vel(3*i-1)*vel(3*i-1)+vel(3*i)*vel(3*i))
    end do
    KE = 0.50d0*mass*KE
    open(2,file='initial_velocities.dat',status='unknown')
    do i=1,npart
       write(2,*) vel(3*i-2), vel(3*i-1), vel(3*i)  
    end do
    close(2)
    deallocate(seed)
 end subroutine velocity_init
 
 !==========CALCULATION OF FORCE==========
 
 subroutine calc_force
    use parameter
    use omp_lib
    implicit none
    PE=0.0d0
    force=0.0d0
    lj_force=0.0d0
    
    do i=1,npart-1
       x1=pos(3*i-2); y1=pos(3*i-1); z1=pos(3*i)
       do j=1,no_neighb(i)
        p=neigh_list(i,j)
          x2=pos(3*p-2); y2=pos(3*p-1); z2=pos(3*p)
          dx=x1-x2; dy=y1-y2; dz=z1-z2
          
          if(abs(dx).ge.llxby2) dx=(llx-abs(dx))*((-1.0d0*dx)/abs(dx))
          if(abs(dy).ge.llyby2) dy=(lly-abs(dy))*((-1.0d0*dy)/abs(dy))
          if(abs(dz).ge.llzby2) dz=(llz-abs(dz))*((-1.0d0*dz)/abs(dz))
          r=dsqrt(dx*dx+dy*dy+dz*dz)
          
          if(r<=rc) then
             lj=eps*(((sigma/r)**12)-((sigma/r)**6))-ufc+fc*r
             PE=PE+lj
             lj_force=eps*((12.0d0*sigma12/(r**13))-(6.0d0*sigma6/(r**7)))-fc
             force(3*i-2)=force(3*i-2)+(lj_force*dx/r)
             force(3*i-1)=force(3*i-1)+(lj_force*dy/r)
             force(3*i)=force(3*i)+(lj_force*dz/r)
             force(3*p-2)=force(3*p-2)-(lj_force*dx/r)
             force(3*p-1)=force(3*p-1)-(lj_force*dy/r)
             force(3*p)=force(3*p)-(lj_force*dz/r)
          end if
       end do
    end do
    acc=force/mass  
 end subroutine calc_force
 
 !==========POSITION UPDATION==========
 
 subroutine update_pos
    use parameter
    use omp_lib
    implicit none
    real*8 :: dt2by2
    old_force=0.0d0
    
    dt2by2 = 0.50d0*delta_t**2    
   ! if(mod(k,100)==0) then
    !write(10,*) npart 
    !write(10,*)
    !end if
    do i=1,npart
       pos(3*i-2)=pos(3*i-2) + (vel(3*i-2)*delta_t) + (dt2by2*force(3*i-2))
       pos(3*i-1)=pos(3*i-1) + (vel(3*i-1)*delta_t) + (dt2by2*force(3*i-1))
       pos(3*i)=pos(3*i) + (vel(3*i)*delta_t) + (dt2by2*force(3*i))
      
       pos(3*i-2)=modulo(pos(3*i-2),llx)
       pos(3*i-1)=modulo(pos(3*i-1),lly)
       pos(3*i)=modulo(pos(3*i),llz)
      
       old_force(3*i-2)=force(3*i-2)
       old_force(3*i-1)=force(3*i-1)
       old_force(3*i)=force(3*i)  
      ! if(mod(k,100)==0) then
     !  write(10,*) pos(3*i-2),pos(3*i-1),pos(3*i)
       !end if
    end do
    old_acc=acc  
 end subroutine update_pos
 
 !==========VELOCITY UPDATION==========
 
 subroutine update_vel
    use parameter
    use omp_lib
    implicit none
    KE=0.0d0
    fs = 0.5d0*delta_t/mass
    
    do i=1,npart
       vel(3*i-2) = vel(3*i-2) +(fs*(old_force(3*i-2)+force(3*i-2)))
       vel(3*i-1) = vel(3*i-1) +(fs*(old_force(3*i-1)+force(3*i-1)))
       vel(3*i) = vel(3*i) +(fs*(old_force(3*i)+force(3*i)))
       KE = KE + (vel(3*i-2)*vel(3*i-2)+vel(3*i-1)*vel(3*i-1)+vel(3*i)*vel(3*i))
    end do
    KE = 0.50d0*mass*KE  
 end subroutine update_vel  
 
 !==========VELOCITY RESCALING(THERMOSTAT)=========
 
 subroutine thermostat
    use parameter
    use omp_lib
    implicit none
 
    if(mod(k,100)==0) then
       theoryKE=1.5d0*dble(npart)*kb_T
       scalef=dsqrt(theoryKE/KE)
       vel = vel*scalef
      
       KE = 0.0d0
       do i=1,npart
          KE = KE + (vel(3*i-2)*vel(3*i-2)+vel(3*i-1)*vel(3*i-1)+vel(3*i)*vel(3*i))
       end do  
       KE = 0.50d0*mass*KE
    end if
 end subroutine thermostat 

 subroutine neighbour_list
    use parameter
    use omp_lib
    implicit none

    
    neigh_list=0 ; no_neighb=0

  
    do i=1,npart-1
        x1=pos(3*i-2); y1=pos(3*i-1); z1=pos(3*i)
        do j=i+1,npart
           x2=pos(3*j-2); y2=pos(3*j-1); z2=pos(3*j)
           dx=x1-x2; dy=y1-y2; dz=z1-z2
           
           if(abs(dx).ge.llxby2) dx=(llx-abs(dx))*((-1.0d0*dx)/abs(dx))
           if(abs(dy).ge.llyby2) dy=(lly-abs(dy))*((-1.0d0*dy)/abs(dy))
           if(abs(dz).ge.llzby2) dz=(llz-abs(dz))*((-1.0d0*dz)/abs(dz))
           r=dsqrt(dx*dx+dy*dy+dz*dz)
          
           if (r<rs) then
           no_neighb(i)=no_neighb(i)+1
           neigh_list(i,no_neighb(i))=j
           end if
        end do
    end do
    !if (20000<k .and. k<30000) then
     ! do i=1 ,npart                           !maximum neighbour list cheking
      !write(10,*) k , maxval(no_neighb)
   !end do
  ! end if

   !do j=1 , no_neighb(i)
     
    !write(10,*) no_neighb(i)
    !end do
   


end subroutine neighbour_list


subroutine compute_gr
    use parameter
    use omp_lib
    implicit none

    
    integer, allocatable :: hist(:)

    
    dr = 0.1d0
    rmax = 0.5d0 * min(llx, lly, llz)
    nbins = int(rmax / dr)

    
    allocate(hist(nbins))
   
    hist = 0

    ! Loop over all unique particle pairs
    do i = 1, npart - 1
        x1 = pos(3*i - 2); y1 = pos(3*i - 1); z1 = pos(3*i)
        do j = i + 1, npart
            x2 = pos(3*j - 2); y2 = pos(3*j - 1); z2 = pos(3*j)

            dx = x1 - x2
            dy = y1 - y2
            dz = z1 - z2

            ! Periodic Boundary Conditions (Minimum Image Convention)
            if(abs(dx).ge.llxby2) dx=(llx-abs(dx))*((-1.0d0*dx)/abs(dx))
            if(abs(dy).ge.llyby2) dy=(lly-abs(dy))*((-1.0d0*dy)/abs(dy))
            if(abs(dz).ge.llzby2) dz=(llz-abs(dz))*((-1.0d0*dz)/abs(dz))
            r=dsqrt(dx*dx+dy*dy+dz*dz)
            

            if (r < rmax) then
                bin = int(r / dr) + 1      ! one add kela because fortran cha array one pasun start hoto
                if (bin <= nbins) hist(bin) = hist(bin) + 2     ! donhi side che particle consider kele ahet
            end if
        end do
    end do

    vol = llx * lly * llz
    rho = dfloat(npart) / vol

    ! Normalize to compute g(r)
    do bin = 1, nbins
        r = (bin - 0.5d0) * dr
        shell_vol = 4.0d0 / 3.0d0 * pi * ((r + dr)**3.0d0 - (r )**3.0d0)
        gr(bin) = gr(bin)+(dfloat(hist(bin)) / (rho * shell_vol * dfloat(npart)))
      !write(*,*)  gr(16)
    end do

    ! Output g(r) to file
   
    

    deallocate( hist)
end subroutine compute_gr
