import numpy as np
import math
import matplotlib.pyplot as plt

# -------------------------------------------------------
# 1. Setup the problem
# -------------------------------------------------------
N = 100               # Number of sub-intervals (0..100 => 101 points)
dx = 1.0 / N          # Grid spacing, h=0.01
x = np.linspace(0, 1, N+1)  # x_0=0, x_N=1 in steps of 0.01

# Boundary conditions
y = np.zeros(N+1, dtype=np.float64)
y[0]   = 0.0   # y(0)
y[N]   = 2.0   # y(1)

# -------------------------------------------------------
# 2. Finite-difference constants
#    Derived from:
#    y''(x_i) - 5 y'(x_i) + 10 y(x_i) = 10 x_i
# -------------------------------------------------------
# Using central differences:
#    y''(x_i) ~ (y_{i+1} - 2y_i + y_{i-1}) / dx^2
#    y'(x_i)  ~ (y_{i+1} - y_{i-1}) / (2 dx)
#
# After rearranging, for i=1..N-1 we get an equation of the form:
#   9750*y_{i+1} + 10250*y_{i-1} - 19990*y_i = 10*x_i
#
a = 9750.0
b = 10250.0
c = -19990.0   # The diagonal coefficient in the rearranged equation

# -------------------------------------------------------
# 3. Gauss-Seidel iteration
# -------------------------------------------------------
tol = 1.0e-4       # Convergence tolerance
maxIter = 100000   # Safety limit on number of iterations

for iteration in range(maxIter):
    diff = 0.0
    # Sweep through interior points
    for i in range(1, N):
        oldVal = y[i]
        # y_i = [a*y_{i+1} + b*y_{i-1} - 10*x_i] / (-c)
        y[i] = (a*y[i+1] + b*y[i-1] - 10.0*x[i]) / (-c)
        # Track largest update
        diff = max(diff, abs(y[i] - oldVal))
    # Check for convergence
    if diff < tol:
        print(f"Converged in {iteration} iterations with max update = {diff:.2e}")
        break

# -------------------------------------------------------
# 4. Exact (closed-form) solution for comparison
# -------------------------------------------------------
# ODE: y'' - 5y' + 10y = 10x
# 1) Solve homogeneous eqn: y_h'' - 5y_h' + 10y_h = 0
#    Characteristic eqn: r^2 - 5r + 10 = 0 => r = (5 +/- i sqrt(15))/2
#    => y_h(x) = e^{(5/2)x}[C1 cos((sqrt(15)/2)x) + C2 sin((sqrt(15)/2)x)]
#
# 2) Find a particular solution: guess y_p = A x + B
#    => A=1, B=0.5
#    => y_p(x) = x + 0.5
#
# => General solution: y(x) = y_h(x) + y_p(x)
# Apply boundary conditions:
#    y(0)=0 => e^0 [C1 cos(0) + C2 sin(0)] + 0.5 = 0 => C1+0.5=0 => C1=-0.5
#    y(1)=2 => e^(5/2)*[ -0.5 cos(sqrt(15)/2) + C2 sin(sqrt(15)/2 )] + 1.5=2
#             => -0.5 cos(...) + C2 sin(...) = 0.5 / e^(5/2)
#             => C2 = [0.5/e^(2.5) + 0.5 cos(...)] / sin(...)
#

C1 = -0.5
cosVal = math.cos(math.sqrt(15)/2)
sinVal = math.sin(math.sqrt(15)/2)
lhs = 0.5 / math.exp(2.5) + 0.5*cosVal
C2 = lhs / sinVal

def exact_solution(x):
    return (math.exp(2.5*x) * 
            (C1*math.cos((math.sqrt(15)/2)*x) + 
             C2*math.sin((math.sqrt(15)/2)*x))
           ) + (x + 0.5)

y_exact = [exact_solution(xi) for xi in x]

# -------------------------------------------------------
# 5. Plot the numerical vs exact solution
# -------------------------------------------------------
plt.figure(figsize=(7,5))
plt.plot(x, y, 'b.-', label='Gauss-Seidel FD')
plt.plot(x, y_exact, 'r--', label='Exact solution')
plt.xlabel('x')
plt.ylabel('y(x)')
plt.title('Solution of y\'\' - 5y\' + 10y = 10x with y(0)=0, y(1)=2')
plt.grid(True)
plt.legend()
plt.show()

# Print the value near x=0.80 for reference
i_80 = 80  # index for x=0.80
print(f"At x = 0.80, numerical y = {y[i_80]:.6f}, exact y = {y_exact[i_80]:.6f}")
