PROGRAM GaussSeidelBVP
    IMPLICIT NONE
    
    INTEGER, PARAMETER :: NMAX = 100        ! number of sub-intervals (0..100 => 101 points)
    REAL*8, PARAMETER :: dx = 0.01D0        ! grid spacing
    REAL*8, PARAMETER :: tol = 1.0D-4       ! convergence tolerance
    INTEGER, PARAMETER :: maxIter = 100000  ! safety limit on iterations
  
    REAL*8 :: x(0:NMAX), y(0:NMAX)
    REAL*8 :: diff, oldVal
    REAL*8 :: a, b, c
    INTEGER :: i, iter
  
    !-------------------------------------------------------------------------
    ! 1. Initialize the grid points and boundary values
    !-------------------------------------------------------------------------
    DO i = 0, NMAX
       x(i) = dx * DBLE(i)
    END DO
    
    ! Boundary conditions:
    !   y(0)   = 0
    !   y(100) = 2
    DO i = 0, NMAX
       y(i) = 0.0D0
    END DO
    y(0)   = 0.0D0
    y(NMAX) = 2.0D0
  
    !-------------------------------------------------------------------------
    ! 2. Constants for the finite-difference equation
    !    9750*y_{i+1} + 10250*y_{i-1} - 19990*y_i = 10*x_i
    !-------------------------------------------------------------------------
    a = 9750.0D0
    b = 10250.0D0
    c = -19990.0D0   ! We will move c to the other side when solving for y_i.
  
    !-------------------------------------------------------------------------
    ! 3. Gauss-Seidel Iteration
    !-------------------------------------------------------------------------
    DO iter = 1, maxIter
  
       diff = 0.0D0   ! track maximum change in y(i) this iteration
  
       ! Sweep over interior points
       DO i = 1, NMAX-1
          oldVal = y(i)
          
          ! y_i = ( a*y_{i+1} + b*y_{i-1} - 10*x_i ) / ( -c )
          ! but -c = 19990, so:
          y(i) = ( a*y(i+1) + b*y(i-1) - 10.0D0*x(i) ) / (-c)
  
          ! Track the largest update
          diff = MAX( diff, ABS(y(i) - oldVal) )
       END DO
  
       ! Check convergence
       IF (diff < tol) THEN
          EXIT
       END IF
  
    END DO
  
    !-------------------------------------------------------------------------
    ! 4. Report the result at x=0.80
    !-------------------------------------------------------------------------
    WRITE(*,'(A,I5)') 'Number of Gauss-Seidel iterations used = ', iter
    WRITE(*,'(A,F12.6)') 'y at x=0.80 is approximately = ', y(80)
  
  END PROGRAM GaussSeidelBVP
