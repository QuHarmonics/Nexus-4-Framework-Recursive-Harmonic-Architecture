import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from IPython.display import HTML

# Parameters
steps = 55
radius_step = 0.15
harmonic_base = 7.75
angle_step = 2 * np.pi / harmonic_base

# Build spiral R → L (future → past)
theta = np.array([(steps - i) * angle_step for i in range(steps)])
radii = np.array([(steps - i) * radius_step for i in range(steps)])
x = radii * np.cos(theta)
y = radii * np.sin(theta)

# Setup plot
fig, ax = plt.subplots(figsize=(6, 6))
spiral, = ax.plot([], [], 'o-', lw=2, color='crimson')
collapse_point, = ax.plot([], [], 'ko', markersize=8)
ax.set_xlim(-2, 4)
ax.set_ylim(-2, 4)
ax.set_aspect('equal')
ax.set_title("Recursive R→L Spiral of Observation Collapse")
ax.axis('off')

# Init function
def init():
    spiral.set_data([], [])
    collapse_point.set_data([], [])
    return spiral, collapse_point

# Frame update function
def update(frame):
    spiral.set_data(x[:frame], y[:frame])
    if frame == steps:
        collapse_point.set_data([x[-1]], [y[-1]])  # must pass as list/tuple
    else:
        collapse_point.set_data([], [])
    return spiral, collapse_point

# Animate
ani = animation.FuncAnimation(fig, update, frames=steps + 1,
                              init_func=init, blit=True, interval=100)

# Show animation
HTML(ani.to_jshtml())
