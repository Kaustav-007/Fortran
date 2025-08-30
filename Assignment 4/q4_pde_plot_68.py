program laplace_sor
    implicit none
 
    !------------------------------------------------------------------------
    ! Parameters matching the C++ code
    !------------------------------------------------------------------------
    integer, parameter :: N = 68
    real*8, parameter :: CONV_LIMIT = 1.0d-5      ! 0.00001
    integer, parameter :: MAX_ITER = 50000
    real*8, parameter :: OMEGA = 1.5d0            ! Relaxation factor
 
    !------------------------------------------------------------------------
    ! Variables
    !------------------------------------------------------------------------
    integer :: i, j, iterations
    real*8  :: dx, dy
    real*8  :: A, B, C, D
    real*8  :: max_diff, shift
    real*8  :: gauss_seidel
    real*8  :: corner11, corner1N, cornerN1, cornerNN
    real*8  :: T(1:N, 1:N), T_old(1:N, 1:N)
 
    !------------------------------------------------------------------------
    ! Initialization
    !------------------------------------------------------------------------
    dx = 1.0d0
    dy = 1.0d0
 
    ! Neumann boundary condition parameters
    A = -70.0d0  ! dT/dx at x = 1
    B = -40.0d0  ! dT/dx at x = 34
    C =  20.0d0  ! dT/dy at y = 1
    D = -10.0d0  ! dT/dy at y = 34
 
    ! Initialize the temperature arrays
    T(:,:)     = 0.0d0
    T_old(:,:) = 0.0d0
 
    ! Initial condition: T(1,1) = 2000
    T(1,1) = 2000.0d0
 
    max_diff   = 1.0d0
    iterations = 0
 
    !------------------------------------------------------------------------
    ! Main SOR Iteration Loop
    !------------------------------------------------------------------------
    do while (max_diff > CONV_LIMIT .and. iterations < MAX_ITER)
       iterations = iterations + 1
 
       ! Save the old temperature
       T_old(:,:) = T(:,:)
 
       !---------------------------------------------------------------------
       ! 1) Update Interior Points with SOR
       !---------------------------------------------------------------------
       do i = 2, N-1
          do j = 2, N-1
             gauss_seidel = 0.25d0 * ( T(i+1,j) + T(i-1,j) + T(i,j+1) + T(i,j-1) )
             T(i,j) = T(i,j) + OMEGA * (gauss_seidel - T(i,j))
          end do
       end do
 
       !---------------------------------------------------------------------
       ! 2) Update Boundaries (excluding corners)
       !---------------------------------------------------------------------
       ! Left boundary (x=1), skip corners
       do j = 2, N-1
          gauss_seidel = 0.25d0 * ( &
               2.0d0*T(2,j) - 2.0d0*dx*A + T(1,j+1) + T(1,j-1) )
          T(1,j) = T(1,j) + OMEGA * (gauss_seidel - T(1,j))
       end do
 
       ! Right boundary (x=N), skip corners
       do j = 2, N-1
          gauss_seidel = 0.25d0 * ( &
               2.0d0*T(N-1,j) + 2.0d0*dx*B + T(N,j+1) + T(N,j-1) )
          T(N,j) = T(N,j) + OMEGA * (gauss_seidel - T(N,j))
       end do
 
       ! Bottom boundary (y=1), skip corners
       do i = 2, N-1
          gauss_seidel = 0.25d0 * ( &
               T(i+1,1) + T(i-1,1) + 2.0d0*T(i,2) - 2.0d0*dy*C )
          T(i,1) = T(i,1) + OMEGA * (gauss_seidel - T(i,1))
       end do
 
       ! Top boundary (y=N), skip corners
       do i = 2, N-1
          gauss_seidel = 0.25d0 * ( &
               T(i+1,N) + T(i-1,N) + 2.0d0*T(i,N-1) + 2.0d0*dy*D )
          T(i,N) = T(i,N) + OMEGA * (gauss_seidel - T(i,N))
       end do
 
       !---------------------------------------------------------------------
       ! 3) Update Corners
       !---------------------------------------------------------------------
       corner11 = 0.5d0 * ( T(1,2) - dy*C + T(2,1) - dx*A )
       T(1,1) = T(1,1) + OMEGA * (corner11 - T(1,1))
 
       corner1N = 0.5d0 * ( T(1,N-1) + dy*D + T(2,N) - dx*A )
       T(1,N) = T(1,N) + OMEGA * (corner1N - T(1,N))
 
       cornerN1 = 0.5d0 * ( T(N-1,1) + dx*B + T(N,2) - dy*C )
       T(N,1) = T(N,1) + OMEGA * (cornerN1 - T(N,1))
 
       cornerNN = 0.5d0 * ( T(N-1,N) + dx*B + T(N,N-1) + dy*D )
       T(N,N) = T(N,N) + OMEGA * (cornerNN - T(N,N))
       !---------------------------------------------------------------------
       ! 4) Normalize to keep T(1,1) = 2000
       !---------------------------------------------------------------------
       shift = 2000.0d0 - T(1,1)
       do i = 1, N
          do j = 1, N
             T(i,j) = T(i,j) + shift
          end do
       end do
       T(1,1) = 2000.0d0
 
       !---------------------------------------------------------------------
       ! 5) Compute max_diff for convergence check
       !---------------------------------------------------------------------
       max_diff = 0.0d0
       do i = 1, N
          do j = 1, N
             if (abs(T(i,j) - T_old(i,j)) > max_diff) then
                max_diff = abs(T(i,j) - T_old(i,j))
             end if
          end do
       end do
 
       ! Optional: print status every 5000 iterations
       if (mod(iterations, 5000) == 0) then
          write(*,'(A,I0,A,ES12.4)') "Iteration ", iterations, &
                ", max_diff = ", max_diff
          write(*,'(A,ES12.4)') "T(20,20) = ", T(20,20)
       end if
 
    end do  ! End of main SOR loop
 
    !------------------------------------------------------------------------
    ! Check whether the loop ended due to convergence or max_iterations
    !------------------------------------------------------------------------
    if (iterations >= MAX_ITER) then
       write(*,*) "Warning: Did not converge within ", MAX_ITER, " iterations!"
    else
       write(*,*) "Converged in ", iterations, " iterations"
    end if
 
    ! Print final temperature at (10,10)
    write(*,'(A,ES12.4)') "Temperature at (20,20): ", T(20,20)
 
    !------------------------------------------------------------------------
    ! Save the final temperature profile to a file
    !------------------------------------------------------------------------
    open(unit=10, file='temperature_profile.txt', status='unknown', action='write')
    write(10,*) "Temperature Profile (34x34 grid):"
    do i = 1, N
       do j = 1, N
          ! Write space-separated values in a single row
          write(10,'(F10.4,1x)', advance='no') T(i,j)
       end do
       write(10,*)  ! New line after each row
    end do
    close(10)
 
    
 
 end program laplace_sor
