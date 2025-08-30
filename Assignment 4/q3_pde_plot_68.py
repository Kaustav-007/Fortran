import numpy as np
import matplotlib.pyplot as plt

# Parameters and Grid Setup
lx = 68  # x-points (1 to 68)
ly = 68  # y-points (1 to 68)
tol = 1e-4  # Convergence tolerance
maxIter = 10000  # Max iterations

# Initialize temperature grid
T = np.zeros((lx, ly))

# Boundary Conditions
T[0, :] = 3.7  # x = 1
T[-1, :] = 0.4  # x = 68
step = (3.7 - 0.4) / (lx - 1)  # 3.3 / 67
for i in range(lx):
    T[i, 0] = 3.7 - i * step  # y = 1
    T[i, -1] = 3.7 - i * step  # y = 68

# Initial guess for interior
for i in range(1, lx-1):
    for j in range(1, ly-1):
        T[i, j] = 0.5 * (T[i, 0] + T[i, -1])

# Gauss-Seidel Iteration
diff = 1.0
it = 0
while diff > tol and it < maxIter:
    diff = 0.0
    for i in range(1, lx-1):
        for j in range(1, ly-1):
            oldVal = T[i, j]
            T[i, j] = 0.25 * (T[i+1, j] + T[i-1, j] + T[i, j+1] + T[i, j-1])
            diff = max(diff, abs(T[i, j] - oldVal))
    it += 1

# Output results
if it == maxIter:
    print("Warning: Max iterations reached.")
else:
    print(f"Converged after {it} iterations, diff = {diff:.2e}")
temp_40_40 = T[39, 39]
print(f"Temperature at (40,40) is {temp_40_40:.4f}")

# Plotting
x = np.linspace(1, 68, lx)
y = np.linspace(1, 68, ly)
X, Y = np.meshgrid(x, y, indexing='ij')
plt.figure(figsize=(8,6))
contour = plt.contourf(X, Y, T, levels=20, cmap='viridis')
plt.colorbar(contour, label='Temperature')
plt.xlabel('x', fontsize=20)
plt.ylabel('y', fontsize=20)
plt.title('Temperature Profile on a 68x68 Plate', fontsize=24)
plt.scatter(40, 40, color='red', marker='o', label=f'(40,40): {temp_40_40:.2f}')
plt.legend()
plt.show()
