program ising_2d
    implicit none
    integer::i,j,k,L,p,a,b,c,d,g,f,niter,time,mm,nn,oo,N
real  ::r,q,E,M,mag,Ei,dE,u,Ef,h
real::T=4.9,J_ising=1.0
integer,dimension(:,:,:),allocatable::spin
integer::seed
character(len=30)::charac_a,charac_b
seed=44859
charac_b='store_config'
print*,'enter the number of lattice points in one dimension'
read*,L
print*,'enter the number of iterations'
read*,niter                                !niter mhanje number of independent microstate that do you want to genrate


allocate(spin(L,L,L))
E=0.0
M=0.0
N=L*L*L

call random_seed
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
!Energy calculate karat ahe atta sum over Si and Sj ithe spin(j,i) che tychya neigbour walya shi interaction consider karat ahet


E=E-J_ising*float((spin(k,j,i))*(spin(a,j,k)+spin(b,j,k)+spin(i,c,k)+spin(i,d,k)+spin(i,j,g)+spin(i,j,f)))
end do
end do
end do

mag=M/(float(N))
E=E*0.5             !extra counting mule 1/2 ne divide kel

print*,'intial energy E,E per spin=',E,E/float(N)
print*,'intial magnetisation M, M per spin=',M,mag
 
! initiallastion  complete zale atta apan thermal reservior shi system la equilibrate karu

open(13,file='q3_ising.dat')

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
            Ei=-J_ising*float((spin(i,j,k))*(spin(a,j,k)+spin(b,j,k)+spin(i,c,k)+spin(i,d,k)+spin(i,j,g)+spin(i,j,f)))
            !trial flip dili
            spin(i,j,k)=-spin(i,j,k)
            Ef=-J_ising*float((spin(i,j,k))*(spin(a,j,k)+spin(b,j,k)+spin(i,c,k)+spin(i,d,k)+spin(i,j,g)+spin(i,j,f)))

            dE=Ef-Ei

            if(dE<=0.0)then
                ! if energy less than zero asel tar accepting this microstate with probability one
                E=E+dE   ! energy and magnetisation update kele   
                M=M+(2.0*float(spin(i,j,k)))   !ek spin flip kelyavar 2 cha change hoto magnetiasion madhe (eg 1,-1,1)
            else
                u=exp(-dE/(T))
                call random_number(h)
                if(h<u) then                       ! he metrapolise algorithm che convention ahe
                    E=E+dE
                    M=M+(2.0*float(spin(i,j,k)))
                else
                    spin(i,j,k)=-spin(i,j,k)    !trail flip accept nahi honar and energy and magnetastion update nahi honar
                end if
            end if
        end do
    end do
end do
    write(13,*)time,M/float(N)
end do
end program ising_2d
