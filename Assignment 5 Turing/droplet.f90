PROGRAM LAPLACE
	IMPLICIT NONE
	INTEGER, PARAMETER 	:: lx=60, ly=60
	REAL*8 				:: old_tempA(lx,ly), tempA(lx,ly)
	REAL*8 				:: old_tempB(lx,ly), tempB(lx,ly)
	INTEGER 			:: i, j, ii, jj, kk
	REAL*8 				:: bound_temp, increment_temp, prefactor, densityA, densityB
	INTEGER				:: xp, xn, yp, yn !Dummy
	
	INTEGER 			:: test, counter, niter
		
	! PARAMETERS For Simulation
	REAL*8 	:: dx=1.0d0, dy=1.0d0, dt=0.0020d0
	REAL*8 	:: toler=0.0001d0
	REAL*8 	:: Diff_a=1.0d0, Diff_b=100.0d0 !Diffusion Constants
	REAL*8 	:: alpha=0.0050d0, beta=10.0d0 !FitzHugh-Nagumo Parameters
	INTEGER :: n_snapshot=1000, tot_niter=40000
	
	OPEN(21,file='initialize.txt',status='unknown')
	OPEN(22,file='initializeB.txt',status='unknown')

	! WRITE(*,*) 'old_temp(1,ly)', old_tempA(1,ly)
	CALL RANDOM_NUMBER(old_tempA)
	CALL RANDOM_NUMBER(old_tempB)
	
	! old_tempA=0.0d0 ; old_tempA(lx/2,ly/2)=1.0d0
	old_tempA = old_tempA - 0.50d0 ; old_tempB = old_tempB - 0.50d0
	tempA = old_tempA ; tempB = old_tempB
	
	! Boundary Conditions at Zeroth Iteration
	niter = 0
	DO ii = 1,lx
		DO jj = 1,ly
			WRITE(21,*) niter*dt, ii, jj, old_tempA(ii,jj)
		END DO
	END DO
	
	DO ii = 1,lx
		DO jj = 1, ly
			WRITE(22,*) niter*dt, ii, jj, old_tempB(ii,jj)
		END DO
	END DO
	WRITE(21,*) ; WRITE(22,*)
	
	! dx=0.05d0 ; dy = 0.05d0
	test=0 ; counter=0 ; prefactor = (0.50d0*dx*dx*dy*dy)/(dx*dx+dy*dy)
	
	DO niter = 1,tot_niter
		! counter = counter + 1
		test = 0
		
		DO jj = 1,ly ! Update Step : Going across a Row
			xp = jj+1 ; xn = jj-1
			
			DO ii = 1,lx ! Going down a Column
			
				IF(jj==1) xn = ly ; if(jj==ly) xp = 1
				yp = ii+1 ; yn = ii - 1 ; if(ii==1) yn = lx ; if(ii==lx) yp = 1
				
				tempA(ii,jj) = old_tempA(ii,jj) + dt*diff_A*(old_tempA(yn,jj) + old_tempA(yp,jj) + &
				& old_tempA(ii,xp) + old_tempA(ii,xn) - 4.0d0*old_tempA(ii,jj))
				
				tempA(ii,jj) = tempA(ii,jj) + dt*(old_tempA(ii,jj) - (old_tempA(ii,jj))**3 + alpha - old_tempB(ii,jj))
				
				tempB(ii,jj) = old_tempB(ii,jj) + dt*diff_B*(old_tempB(yn,jj) + old_tempB(yp,jj) + &
				& old_tempB(ii,xp) + old_tempB(ii,xn) - 4.0d0*old_tempB(ii,jj))
				
				tempB(ii,jj) = tempB(ii,jj) + dt*beta*(old_tempA(ii,jj) - old_tempB(ii,jj))
			END DO
			! tempB=0.0d0
		END DO
		
		old_tempA = tempA ! After Test Condition
		old_tempB = tempB ! After Test Condition
		
		! DO jj = 2,ly-1 ! Checking for Convergence at Each Lattice Site
			! DO ii = 2,lx-1
				! IF((ABS(temp(jj,ii) - old_temp(jj,ii))) > toler) test = 1
			! END DO
		! END DO
		
		! IF(test == 0) EXIT ! EXIT Condition
		
		! Write the Evolving Equation RESULT : Solution to the Equation with PBC
		
		IF(MOD(niter, n_snapshot) == 0) THEN
		
			densityA = 0.0d0 ; densityB =0.0d0
			
			DO ii = 1,lx
				DO jj = 1,ly
					WRITE(21,*) niter*dt, ii, jj, tempA(ii,jj)
					WRITE(22,*) niter*dt, ii, jj, tempB(ii,jj)
					densityA = densityA + tempA(ii,jj)
					densityB = densityB + tempB(ii,jj)
				END DO
			END DO
			WRITE(*,*) niter, densityA, densityB
			WRITE(21,*) ; WRITE(22,*)
		ENDIF
	END DO ! do niter = 1, tot_niter
	
	CLOSE(21)
	CLOSE(22)	
END PROGRAM LAPLACE
