module md_simulation
    implicit none
    integer, parameter :: dp = kind(1.0d0)
    integer, parameter :: N = 2139
    integer, parameter :: box_size = 20
    real(dp), parameter :: dt = 0.005_dp
    real(dp), parameter :: eps = 1.0_dp
    real(dp), parameter :: sigma = 1.0_dp
    real(dp), parameter :: rc = 2.5_dp
    real(dp), parameter :: rc2 = rc * rc
    real(dp), parameter :: rs = 4.5_dp
    real(dp), parameter :: rs2 = rs * rs
    real(dp), parameter :: kBT = 1.0_dp
    integer, parameter :: nl_update_freq = 40
    integer, parameter :: equil_steps = 50000
    integer, parameter :: prod_steps = 250000
    integer, parameter :: sample_freq = 100
    real(dp), parameter :: v_max = 5.0_dp
    real(dp), parameter :: dv = 0.05_dp
    integer, parameter :: n_vbins = int(v_max / dv)
  
  contains
  
    subroutine initialize_particles(pos, vel)
      real(dp), intent(out) :: pos(N,3), vel(N,3)
      integer :: i, j, k, idx, n_per_side
      real(dp) :: spacing, v_scale
      idx = 0
      n_per_side = nint((N / 4.0_dp)**(1.0_dp/3.0_dp))
      spacing = box_size / real(n_per_side, dp)
      do i = 0, n_per_side-1
        do j = 0, n_per_side-1
          do k = 0, n_per_side-1
            if (idx < N) then
              idx = idx + 1
              pos(idx,:) = [real(i, dp), real(j, dp), real(k, dp)] * spacing
            end if
            if (idx < N) then
              idx = idx + 1
              pos(idx,:) = [real(i, dp)+0.5_dp, real(j, dp), real(k, dp)+0.5_dp] * spacing
            end if
            if (idx < N) then
              idx = idx + 1
              pos(idx,:) = [real(i, dp), real(j, dp)+0.5_dp, real(k, dp)+0.5_dp] * spacing
            end if
            if (idx < N) then
              idx = idx + 1
              pos(idx,:) = [real(i, dp)+0.5_dp, real(j, dp)+0.5_dp, real(k, dp)] * spacing
            end if
          end do
        end do
      end do
      pos = modulo(pos, real(box_size, dp))
      call random_number(vel)
      vel = vel - 0.5_dp
      do i = 1, 3
        vel(:,i) = vel(:,i) - sum(vel(:,i)) / real(N, dp)
      end do
      v_scale = sqrt(3.0_dp * kBT * N / sum(vel**2))
      vel = vel * v_scale
    end subroutine initialize_particles
  
    subroutine build_neighbor_list(pos, neighbor_list, num_neighbors)
      real(dp), intent(in) :: pos(N,3)
      integer, intent(out) :: neighbor_list(N,N), num_neighbors(N)
      real(dp) :: rij(3), dist2
      integer :: i, j
      neighbor_list = 0
      num_neighbors = 0
      do i = 1, N-1
        do j = i+1, N
          rij = pos(i,:) - pos(j,:)
          rij = rij - box_size * nint(rij / box_size)
          dist2 = sum(rij**2)
          if (dist2 < rs2) then
            num_neighbors(i) = num_neighbors(i) + 1
            neighbor_list(i, num_neighbors(i)) = j
            num_neighbors(j) = num_neighbors(j) + 1
            neighbor_list(j, num_neighbors(j)) = i
          end if
        end do
      end do
    end subroutine build_neighbor_list
  
    subroutine compute_forces(pos, neighbor_list, num_neighbors, forces, pe)
      real(dp), intent(in) :: pos(N,3)
      integer, intent(in) :: neighbor_list(N,N), num_neighbors(N)
      real(dp), intent(out) :: forces(N,3), pe
      real(dp) :: rij(3), dist2, inv_dist2, inv_dist6, f_scalar
      integer :: i, j, k
      forces = 0.0_dp
      pe = 0.0_dp
      do i = 1, N
        do k = 1, num_neighbors(i)
          j = neighbor_list(i,k)
          if (j > i) then
            rij = pos(i,:) - pos(j,:)
            rij = rij - box_size * nint(rij / box_size)
            dist2 = sum(rij**2)
            if (dist2 < rc2 .and. dist2 > 0.0_dp) then
              inv_dist2 = 1.0_dp / dist2
              inv_dist6 = inv_dist2**3
              f_scalar = 48.0_dp * eps * inv_dist6 * (inv_dist6 - 0.5_dp) * inv_dist2
              forces(i,:) = forces(i,:) + f_scalar * rij
              forces(j,:) = forces(j,:) - f_scalar * rij
              pe = pe + 4.0_dp * eps * (inv_dist6 * (inv_dist6 - 1.0_dp))
            end if
          end if
        end do
      end do
    end subroutine compute_forces
  
    subroutine apply_thermostat(vel)
      real(dp), intent(inout) :: vel(N,3)
      real(dp) :: ke, temp, scale
      ke = 0.5_dp * sum(vel**2)
      temp = 2.0_dp * ke / (3.0_dp * real(N, dp))
      scale = sqrt(kBT / temp)
      vel = vel * scale
    end subroutine apply_thermostat
  
    subroutine velocity_verlet(pos, vel, forces, neighbor_list, num_neighbors, pe, ke, step)
      real(dp), intent(inout) :: pos(N,3), vel(N,3), forces(N,3)
      integer, intent(inout) :: neighbor_list(N,N), num_neighbors(N)
      real(dp), intent(out) :: pe, ke
      integer, intent(in) :: step
      real(dp) :: new_forces(N,3)
      pos = pos + vel * dt + 0.5_dp * forces * dt**2
      pos = modulo(pos, real(box_size, dp))
      if (mod(step, nl_update_freq) == 0) then
        call build_neighbor_list(pos, neighbor_list, num_neighbors)
      end if
      call compute_forces(pos, neighbor_list, num_neighbors, new_forces, pe)
      vel = vel + 0.5_dp * (forces + new_forces) * dt
      forces = new_forces
      ke = 0.5_dp * sum(vel**2)
      if (isnan(ke)) then
        print *, "Error: NaN in KE at step ", step
        stop
      end if
    end subroutine velocity_verlet
  
    subroutine compute_speed_distribution(vel, v_hist, v_count)
      real(dp), intent(in) :: vel(N,3)
      real(dp), intent(inout) :: v_hist(n_vbins)
      integer, intent(inout) :: v_count
      real(dp) :: speed
      integer :: i, bin
      do i = 1, N
        speed = sqrt(sum(vel(i,:)**2))
        if (speed < v_max) then
          bin = int(speed / dv) + 1
          if (bin <= n_vbins) then
            v_hist(bin) = v_hist(bin) + 1.0_dp
          end if
        end if
      end do
      v_count = v_count + 1
    end subroutine compute_speed_distribution
  
    subroutine finalize_speed_distribution(v_hist, v_count)
      real(dp), intent(in) :: v_hist(n_vbins)
      integer, intent(in) :: v_count
      real(dp) :: v, f_v, norm
      integer :: i
      real(dp) :: v_dist(n_vbins)
      open(10, file='speed_dist.dat', status='replace')
      norm = real(N, dp) * real(v_count, dp) * dv
      do i = 1, n_vbins
        v = (i - 0.5_dp) * dv
        v_dist(i) = v_hist(i) / norm
        f_v = sqrt(2.0_dp / 3.14159265359_dp) * v**2 * exp(-0.5_dp * v**2)
        write(10, *) v, v_dist(i), f_v
      end do
      close(10)
    end subroutine finalize_speed_distribution
  end module md_simulation
  
  program main
    use md_simulation
    implicit none
    real(dp) :: pos(N,3), vel(N,3), forces(N,3)
    integer :: neighbor_list(N,N), num_neighbors(N)
    real(dp) :: pe, ke, v_hist(n_vbins)
    integer :: step, v_count
    call initialize_particles(pos, vel)
    neighbor_list = 0
    num_neighbors = 0
    call build_neighbor_list(pos, neighbor_list, num_neighbors)
    call compute_forces(pos, neighbor_list, num_neighbors, forces, pe)
    ke = 0.5_dp * sum(vel**2)
    v_hist = 0.0_dp
    v_count = 0
    do step = 1, equil_steps
      call velocity_verlet(pos, vel, forces, neighbor_list, num_neighbors, pe, ke, step)
      call apply_thermostat(vel)
      if (mod(step, 1000) == 0) then
        print *, "Equil Step:", step, "PE:", pe, "KE:", ke, "Total E:", pe + ke
      end if
    end do
    do step = 1, prod_steps
      call velocity_verlet(pos, vel, forces, neighbor_list, num_neighbors, pe, ke, step)
      if (mod(step, sample_freq) == 0) then
        call compute_speed_distribution(vel, v_hist, v_count)
      end if
      if (mod(step, 1000) == 0) then
        print *, "Prod Step:", step, "PE:", pe, "KE:", ke, "Total E:", pe + ke
      end if
    end do
    call finalize_speed_distribution(v_hist, v_count)
  end program main