program laplace_q3
    implicit none
  
    !-----------------------------------------------------------------------
    ! 1. Parameters and Declarations
    !-----------------------------------------------------------------------
    integer, parameter :: lx = 68, ly = 68    ! Changed from 34 to 68
    integer, parameter :: maxIter = 10000
    real(8), parameter :: tol = 1.0d-4        ! Convergence limit = 0.0001
  
    real(8) :: T(1:lx, 1:ly)
    real(8) :: diff, oldVal
    integer :: i, j, iter
  
    !-----------------------------------------------------------------------
    ! 2. Initialize Arrays and Set Boundary Conditions
    !-----------------------------------------------------------------------
    ! Set the interior initial guess to zero.
    T = 0.0d0
  
    ! Left boundary (x=1): For all j, T(1,j) = 3.7.
    do j = 1, ly
       T(1,j) = 3.7d0
    end do
  
    ! Right boundary (x=lx): For all j, T(68,j) = 0.4.
    do j = 1, ly
       T(lx,j) = 0.4d0
    end do
  
    ! Bottom boundary (y=1): T decreases linearly from 3.7 to 0.4 as x changes from 1 to 68.
    ! Since there are 67 intervals in x, the decrement per step is:
    !   delta = (3.7 - 0.4) / (lx - 1) = 3.3/67 ≈ 0.049253731
    do i = 1, lx
       T(i,1) = 3.7d0 - (3.3d0/67.0d0)*(i-1)
    end do
  
    ! Top boundary (y=ly): Similarly, T decreases linearly from 3.7 to 0.4.
    do i = 1, lx
       T(i,ly) = 3.7d0 - (3.3d0/67.0d0)*(i-1)
    end do
  
    !-----------------------------------------------------------------------
    ! 3. Gauss–Seidel Iteration (In-Place Update)
    !-----------------------------------------------------------------------
    ! The interior nodes (i=2:lx-1, j=2:ly-1) are updated by the 5-point formula.
    write(*,*) "Starting Gauss-Seidel iteration on a ", lx, "x", ly, " grid..."
    
    diff = 1.0d0
    iter = 0
  
    do while (diff > tol .and. iter < maxIter)
       diff = 0.0d0
  
       ! Update interior points
       do j = 2, ly - 1
          do i = 2, lx - 1
             oldVal = T(i,j)
             T(i,j) = 0.25d0 * ( T(i+1,j) + T(i-1,j) + T(i,j+1) + T(i,j-1) )
             diff = max(diff, dabs(T(i,j) - oldVal))
          end do
       end do
  
       iter = iter + 1
    end do
  
    if (diff > tol) then
       write(*,*) "Warning: did not fully converge after ", maxIter, " iterations."
    else
       write(*,*) "Converged after ", iter, " iterations, final diff = ", diff
    end if
  
    !-----------------------------------------------------------------------
    ! 4. Output the Temperature at (40,40)
    !-----------------------------------------------------------------------
    write(*,*) "The temperature at (40,40) is ", T(40,40)    ! Changed from (20,20) to (40,40)
  
end program laplace_q3
