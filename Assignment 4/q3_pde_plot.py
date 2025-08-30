import numpy as np
import matplotlib.pyplot as plt

#-----------------------------------------------------------
# 1. Parameters and Grid Setup
#-----------------------------------------------------------
lx = 34  # number of points in x-direction (x = 1 ... 34)
ly = 34  # number of points in y-direction (y = 1 ... 34)
tol = 1e-4       # Convergence tolerance
maxIter = 10000  # Maximum number of iterations

# Create a 2D array for the temperature field
T = np.zeros((lx, ly))

#-----------------------------------------------------------
# 2. Apply Boundary Conditions
#-----------------------------------------------------------
# Left boundary: x = 1, T = 3.7 for all y
T[0, :] = 3.7

# Right boundary: x = 34, T = 0.4 for all y
T[-1, :] = 0.4

# Bottom boundary: y = 1, T decreases linearly from 3.7 at x=1 to 0.4 at x=34.
# There are lx points, so the decrement per step is (3.7-0.4)/(lx-1) = 3.3/33 = 0.1.
for i in range(lx):
    T[i, 0] = 3.7 - i * 0.1

# Top boundary: y = 34, similarly T decreases linearly from 3.7 to 0.4.
for i in range(lx):
    T[i, -1] = 3.7 - i * 0.1

# Initialize the interior (optional initial guess):
# Here, we use the average of the bottom and top boundaries for interior points.
for i in range(1, lx-1):
    for j in range(1, ly-1):
        T[i, j] = 0.5 * (T[i, 0] + T[i, -1])

#-----------------------------------------------------------
# 3. Gauss–Seidel Iteration
#-----------------------------------------------------------
diff = 1.0
it = 0

while diff > tol and it < maxIter:
    diff = 0.0
    # Update interior points only
    for i in range(1, lx-1):
        for j in range(1, ly-1):
            oldVal = T[i, j]
            T[i, j] = 0.25 * (T[i+1, j] + T[i-1, j] + T[i, j+1] + T[i, j-1])
            diff = max(diff, abs(T[i, j] - oldVal))
    it += 1

if it == maxIter:
    print("Warning: Maximum iterations reached without full convergence.")
else:
    print(f"Converged after {it} iterations, final diff = {diff:.2e}")

# Temperature at (20,20) -> index [19,19]
temp_20_20 = T[19, 19]
print(f"The temperature at (20,20) is {temp_20_20:.4f}")

#-----------------------------------------------------------
# 4. Plotting the Temperature Profile
#-----------------------------------------------------------
# Create grid arrays corresponding to x and y coordinates (from 1 to 34)
x = np.linspace(1, 34, lx)
y = np.linspace(1, 34, ly)
X, Y = np.meshgrid(x, y, indexing='ij')

plt.figure(figsize=(8,6))
contour = plt.contourf(X, Y, T, levels=20, cmap='viridis')
plt.colorbar(contour, label='Temperature')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Temperature Profile on a 34x34 Plate')

# Mark the point (20,20)
plt.scatter(20, 20, color='red', marker='o', label=f'(20,20): {temp_20_20:.2f}')
plt.legend()
plt.show()
