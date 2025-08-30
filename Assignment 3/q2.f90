!Program to implememnt 3D ising model in a spin lattice

implicit none

integer:: i,j,k,L,p,a,b,c,d,e,f,niter,time,nn,mm,oo,N

real::r,q,En,M,Ef,dE,u,Ei,h,Mag

real:: T=4.90,J_ising=1.0

integer, dimension(:,:,:),allocatable::spin

character(len=30):: charac_a,charac_b

charac_b='store_config'

print*,'Enter the number of lattice points in one dimension'
read*,L

print*, 'Enter the number of iterations'
read*,niter
allocate(spin(L,L,L))

En=0.0 ! initialize the energy value which is instantenous

M=0.0 ! initialize the value of magnetisationwhich is instanteneous
N=L*L*L

open(11, file='three_D ising.dat')

p=0
do i=1,L

   do j=1,L

      do k=1,L
call random_number(r)
spin(k,j,i)=+1
! if (r<0.5) then
!             spin(k,j,i)=-1
! else
!     spin(k,j,i)=+1
!         end if

write(10,fmt='(4g10.8)') float(i),float(j),float(k),float(p),float(spin(k,j,i))
      end do
   end do
end do
close(11)

! calculate the initial energy and magnetisation
do i =1,L
   do j =1,L
      do k=1,L
      a=i+1;b=i-1;c=j+1;d=j-1;e=k+1;f=k-1
     
      if(i==L)a=1
      if(i==1)b=L
      if(j==L)c=1
      if(j==1)d=L
      if(k==L)e=1
      if(k==1)f=L
     
      M=M+spin(i,j,k)
      En=En-J_ising*float((spin(k,j,i))*(spin(a,j,k)+spin(b,j,k)+spin(i,c,k)+spin(i,d,k)+spin(i,j,e)+spin(i,j,f)))
     
      end do
   end do
     
         
end do
Mag=M/(float(N))
En=En*0.5d0

print*,'initial energy E,E per spin=',En,En/float(N)
print*,'initial magnetisation M, M per spin=',M,Mag

!Evolve it to reach Equilibrium

open(11,file='ising_T2_N40_init_random.dat')

do time=1,niter! loop over no of MCS
   do mm=1,L
      do nn=1,L
         do oo=1,L
            call random_number(r);  i=int(r*float(L))+1
            call random_number(r);  j=int(r*float(L))+1
            call random_number(r);  k=int(r*float(L))+1
            a=i+1; b=i-1; c=j+1; d=j-1; e=k+1; f=k-1
            if(i==L)a=1; if(i==1)b=L;   if(j==L)c=1; if(j==1)d=L; if(k==L)e=1; if(k==1)f=L
            Ei=-J_ising*float((spin(k,j,i))*(spin(a,j,k)+spin(b,j,k)+spin(i,c,k)+spin(i,d,k)+spin(i,j,e)+spin(i,j,f)))! Energy before the spin flip
            spin(i,j,k)=-spin(i,j,k) ! Trial flip
            Ef=-J_ising*float((spin(k,j,i))*(spin(a,j,k)+spin(b,j,k)+spin(i,c,k)+spin(i,d,k)+spin(i,j,e)+spin(i,j,f))) ! Energy after the flip
           
            dE=Ef-Ei !Difference in the energies
            if(dE<=0.0) then
               En=En+dE
               M=M+(2.0*float(spin(i,j,k)))
            else
                u=exp(-dE/(T))
                call random_number(h)
                if(h<u) then
                   En=En+dE
                   M=M+(2.0*float(spin(i,j,k)))
                else
                    spin(i,j,k)=spin(i,j,k)
                end if
            end if
               
         end do
      end do
   end  do
   write(11,*)time,M/float(N),En/float(N)
end do


end program
