import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Parameters
N = 50  # number of particles
dt = 0.02  # time step
steps = 2000  # total steps, t = 40
radius = 5.0  # ring radius

# Initial conditions
state = np.zeros(2 * N)  # [y1, v1, y2, v2, ..., y50, v50]
state[0] = 0.8  # y1(0) = 0.8
state[50] = 0.8  # y26(0) = 0.8 (index 50 = 2*26-2)

# Derivative function
def derivatives(state):
    dydt = np.zeros(2 * N)
    for i in range(N):
        idx_y = 2 * i
        idx_v = 2 * i + 1
        dydt[idx_y] = state[idx_v]  # dy/dt = v
        i_next = (i + 1) % N
        i_prev = (i - 1) % N
        dydt[idx_v] = (state[2 * i_next] + state[2 * i_prev] - 2 * state[idx_y])  # dv/dt = accel
    return dydt

# RK4 step
def rk4_step(state, dt):
    k1 = derivatives(state)
    k2 = derivatives(state + (dt / 2) * k1)
    k3 = derivatives(state + (dt / 2) * k2)
    k4 = derivatives(state + dt * k3)
    return state + (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4)

# Simulation
positions = [state[0::2].copy()]  # Store y_i for each step
current_state = state.copy()
for _ in range(steps):
    current_state = rk4_step(current_state, dt)
    positions.append(current_state[0::2].copy())

# Result at t = 40
y1_final = current_state[0]
print(f"Position of particle 1 at t = 40: y_1 = {y1_final:.6f}")

# Ring coordinates
theta = np.linspace(0, 2 * np.pi, N, endpoint=False)
x = radius * np.cos(theta)
z = radius * np.sin(theta)

# Setup figure
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')
ax.set_xlim(-6, 6)
ax.set_ylim(-1.5, 1.5)  # y-axis range reflects oscillation amplitude
ax.set_zlim(-6, 6)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
scatter = ax.scatter(x, positions[0], z, c='b', s=50)

# Update function for animation
def update(frame):
    y = positions[frame]
    scatter._offsets3d = (x, y, z)
    ax.set_title(f"Time: {frame * dt:.2f} s")
    return scatter,

# Create animation
ani = FuncAnimation(fig, update, frames=range(0, steps, 10), interval=50, blit=True)

# Display animation
plt.show()

# Optional: Save animation (uncomment and ensure ffmpeg is installed)
ani.save('ring_oscillation.mp4', writer='ffmpeg', dpi=300)
