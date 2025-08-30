import numpy as np
import matplotlib.pyplot as plt

# Grid size and constants
N = 34
CONV_LIMIT = 1.0e-5   # 0.00001
MAX_ITER = 50000
OMEGA = 1.5

# Physical parameters
dx = 1.0
dy = 1.0
A = -70.0  # dT/dx at x=0   (Fortran's x=1)
B = -40.0  # dT/dx at x=N-1 (Fortran's x=34)
C =  20.0  # dT/dy at y=0   (Fortran's y=1)
D = -10.0  # dT/dy at y=N-1 (Fortran's y=34)

# Allocate temperature arrays
T     = np.zeros((N, N), dtype=float)
T_old = np.zeros((N, N), dtype=float)

# Initial condition: T(0,0) = 2000  (matches Fortran's T(1,1))
T[0, 0] = 2000.0

max_diff   = 1.0
iterations = 0

# Main SOR loop
while max_diff > CONV_LIMIT and iterations < MAX_ITER:
    iterations += 1
    T_old[:, :] = T[:, :]

    #----------------------------------------------------------
    # 1) Update interior points with SOR (i=1..N-2, j=1..N-2)
    #----------------------------------------------------------
    for i in range(1, N-1):
        for j in range(1, N-1):
            gauss_seidel = 0.25 * (T[i+1, j] + T[i-1, j] + T[i, j+1] + T[i, j-1])
            T[i, j] += OMEGA * (gauss_seidel - T[i, j])

    #----------------------------------------------------------
    # 2) Update boundaries (excluding corners)
    #    Left boundary (x=0), Right boundary (x=N-1),
    #    Bottom boundary (y=0), Top boundary (y=N-1)
    #----------------------------------------------------------

    # Left boundary (x=0), skip corners
    for j in range(1, N-1):
        gauss_seidel = 0.25 * (
            2.0 * T[1, j] - 2.0 * dx * A
            + T[0, j+1] + T[0, j-1]
        )
        T[0, j] += OMEGA * (gauss_seidel - T[0, j])

    # Right boundary (x=N-1), skip corners
    for j in range(1, N-1):
        gauss_seidel = 0.25 * (
            2.0 * T[N-2, j] + 2.0 * dx * B
            + T[N-1, j+1] + T[N-1, j-1]
        )
        T[N-1, j] += OMEGA * (gauss_seidel - T[N-1, j])

    # Bottom boundary (y=0), skip corners
    for i in range(1, N-1):
        gauss_seidel = 0.25 * (
            T[i+1, 0] + T[i-1, 0]
            + 2.0 * T[i, 1] - 2.0 * dy * C
        )
        T[i, 0] += OMEGA * (gauss_seidel - T[i, 0])

    # Top boundary (y=N-1), skip corners
    for i in range(1, N-1):
        gauss_seidel = 0.25 * (
            T[i+1, N-1] + T[i-1, N-1]
            + 2.0 * T[i, N-2] + 2.0 * dy * D
        )
        T[i, N-1] += OMEGA * (gauss_seidel - T[i, N-1])

    #----------------------------------------------------------
    # 3) Update corners
    #----------------------------------------------------------
    corner11 = 0.5 * (T[0, 1] - dy*C + T[1, 0] - dx*A)
    T[0, 0] += OMEGA * (corner11 - T[0, 0])

    corner1N = 0.5 * (T[0, N-2] + dy*D + T[1, N-1] - dx*A)
    T[0, N-1] += OMEGA * (corner1N - T[0, N-1])

    cornerN1 = 0.5 * (T[N-2, 0] + dx*B + T[N-1, 1] - dy*C)
    T[N-1, 0] += OMEGA * (cornerN1 - T[N-1, 0])

    cornerNN = 0.5 * (T[N-2, N-1] + dx*B + T[N-1, N-2] + dy*D)
    T[N-1, N-1] += OMEGA * (cornerNN - T[N-1, N-1])

    #----------------------------------------------------------
    # 4) Normalize so T(0,0) stays at 2000
    #----------------------------------------------------------
    shift = 2000.0 - T[0, 0]
    T += shift
    T[0, 0] = 2000.0

    #----------------------------------------------------------
    # 5) Compute max_diff for convergence
    #----------------------------------------------------------
    diff_array = np.abs(T - T_old)
    max_diff = np.max(diff_array)

    # Optional progress print every 5000 iterations
    if iterations % 5000 == 0:
        print(f"Iteration {iterations}, max_diff = {max_diff:.4e}")
        print(f"T(10,10) = {T[9,9]:.4f}")  # Fortran (10,10) -> Python (9,9)

#--------------------------------------------------------------------
# End of SOR loop
#--------------------------------------------------------------------
if iterations >= MAX_ITER:
    print(f"Warning: Did not converge within {MAX_ITER} iterations!")
else:
    print(f"Converged in {iterations} iterations")

print(f"Temperature at (10,10): {T[9,9]:.4f}")

#--------------------------------------------------------------------
# Plot the final temperature distribution
#--------------------------------------------------------------------
plt.figure(figsize=(7,6))
plt.imshow(T.T, origin='lower', cmap='hot', interpolation='nearest')
plt.colorbar(label="Temperature")
plt.title("Temperature Profile from SOR")
plt.xlabel("X coordinate")
plt.ylabel("Y coordinate")
plt.show()
