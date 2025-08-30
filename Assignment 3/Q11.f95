program ising_2d
    implicit none
    integer::i,j,k,L,p,a,b,c,d,g,f,niter,time,mm,nn,oo,N,temp
real*8 ::r,E,M,mag,Ei,dE,u,Ef,h
real*8 ::T,J_ising=1.0
real*8:: final_m4,av_m4,av_m2,final_m2_2,result
integer,dimension(:,:,:),allocatable::spin
integer::n_equilib,n_stat

print*,'enter the number of lattice points in one dimension'
read*,L
print*,'enter the number of iterations'
read*,niter                                !niter mhanje number of independent microstate that do you want to genrate


allocate(spin(L,L,L))
E=0.0
M=0.0
N=L*L*L

n_equilib=10000
n_stat=10
p=0
do i=1,L
    do j=1,L
        do k=1,L
        call random_number(r)
       ! spin(k,j,i)=1
       if(r<0.5) then 
            spin(k,j,i)=-1
        else 
            spin(k,j,i)=+1
        end if
    end do
    end do
end do
!now we calculating initial magnetization and energy

do i=1,L
    do j=1,L
        do k=1,L
    
        ! we identify the neighbour
    
        a=i+1;b=i-1;c=j+1;d=j-1;g=k+1;f=k-1   
    
        if(i==L) a=1
        if(i==1) b=L            !periodic boundary condition
        if(j==1) d=L
        if(j==L) c=1
        if(k==1) f=L
        if(k==L) g=1
        
    
    
    M=M+spin(i,j,k)
    
    
    
    E=E-J_ising*dfloat((spin(k,j,i))*(spin(a,j,k)+spin(b,j,k)+spin(i,c,k)+spin(i,d,k)+spin(i,j,g)+spin(i,j,f)))
    end do
    end do
    end do
    
    mag=M/(dfloat(N))
    E=E*0.5             !extra counting mule 1/2 ne divide kel
    
    print*,'intial energy E,E per spin=',E,E/dfloat(N)
    print*,'intial magnetisation M, M per spin=',M,mag    

! initiallastion  complete zale atta apan thermal reservior shi system la equilibrate karu

open(10,file='q11next_for_L7.dat')

do temp= 470,380,-2

    T= (dfloat(temp))/100.0d0
    final_m4=0.0d0 ; final_m2_2 =0.0d0 ;av_m2=0.0d0 ;av_m4=0.0d0 ; result=0.0d0

do time=1,niter
     do mm=1,L
        do nn=1,L
            do oo=1,L
            call random_number(r); i=int(r*float(L))+1       !choose randomly lattice in x direction 
            call random_number(r); j=int(r*float(L))+1      !choose random lattice in y direction
            call random_number(r); k=int(r*float(L))+1
            a=i+1;b=i-1;c=j+1;d=j-1;g=k+1;f=k-1  

            if(i==L) a=1
            if(i==1) b=L            !periodic boundary condition
            if(j==1) d=L
            if(j==L) c=1
            if(k==1) f=L
            if(k==L) g=1
           !metrapholis alogorithm sathi energy calculate karat ahe 
            !intial energy
            Ei=-J_ising*dfloat((spin(i,j,k))*(spin(a,j,k)+spin(b,j,k)+spin(i,c,k)+spin(i,d,k)+spin(i,j,g)+spin(i,j,f)))
            !trial flip dili
            spin(i,j,k)=-spin(i,j,k)
            Ef=-J_ising*dfloat((spin(i,j,k))*(spin(a,j,k)+spin(b,j,k)+spin(i,c,k)+spin(i,d,k)+spin(i,j,g)+spin(i,j,f)))

            dE=Ef-Ei

            if(dE<=0.0)then
                ! if energy less than zero asel tar accepting this microstate with probability one
                E=E+dE   ! energy and magnetisation update kele   
                M=M+(2.0*dfloat(spin(i,j,k)))   !ek spin flip kelyavar 2 cha change hoto magnetiasion madhe (eg 1,-1,1)
            else
                u=exp(-dE/(T))
                call random_number(h)
                if(h<u) then                       
                    E=E+dE
                    M=M+(2.0*dfloat(spin(i,j,k)))
                else
                    spin(i,j,k)=-spin(i,j,k)    
                end if
            end if
        end do
    end do
end do
    if (time>n_equilib)then
       ! if (mod(time,n_stat).eq.0) then
          av_m4=av_m4+(M*M*M*M)
            av_m2=av_m2+(M*M) 
       ! end if
    end if
end do

final_m2_2= 3*(av_m2/(dfloat(niter-n_equilib)))**2
final_m4=(av_m4/(dfloat(niter-n_equilib)))
result=1-(final_m4/final_m2_2)

write(10,*) T,result
end do
close(10)
deallocate(spin)
print*,"your file is ready"
end program ising_2d
