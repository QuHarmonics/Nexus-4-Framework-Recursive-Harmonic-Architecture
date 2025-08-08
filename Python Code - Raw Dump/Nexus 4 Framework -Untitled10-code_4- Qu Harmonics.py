# Animated ψ-Vector Convergence Toward 0.35
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Phase-space attractor
psi_target = 0.35

# Simulated Δ-trajectory (from decoded digit behavior)
initial_points = [0.10, 0.20, 0.50, 0.75, 0.60, 0.90, 0.25, 0.40]

# Convergence dynamics
def update_points(points, alpha=0.15):
    return [p - alpha * (p - psi_target) for p in points]

# Set up animation
fig, ax = plt.subplots(figsize=(8, 6))
x_base = np.linspace(0, 1, len(initial_points))
scat = ax.scatter(x_base, initial_points, color='purple')

ax.axhline(y=psi_target, color='green', linestyle='--', label='ψ = 0.35')
ax.set_ylim(0, 1)
ax.set_title("ψ-Vector Convergence to Attractor at 0.35")
ax.set_ylabel("ψ-State")
ax.legend()

state = [initial_points[:]]

def animate(i):
    state[0] = update_points(state[0])
    scat.set_offsets(np.c_[x_base, state[0]])
    return scat,

ani = animation.FuncAnimation(fig, animate, frames=30, interval=300, blit=True)
plt.show()