# Fix: Handle case when no centroids should be plotted yet
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
from IPython.display import HTML

# Spiral parameters
total_harmonic_states = 55
harmonic_base = 7.75
angle_step = 2 * np.pi / harmonic_base
radius_step = 0.15

# Triangle centroid points (from earlier harmonic sequence)
centroid_points = [
    (1.6667, 0.0),         # π-ray
    (2.20833, 0.48412),    # 2-3-4
    (1.5, 0.44096),        # 2-2-3
    (2.6667, 1.0)          # 3-4-5
]

# Generate spiral points
theta = np.array([i * angle_step for i in range(total_harmonic_states)])
radii = np.array([i * radius_step for i in range(total_harmonic_states)])
x = radii * np.cos(theta)
y = radii * np.sin(theta)

# Setup plot
fig, ax = plt.subplots(figsize=(6, 6))
line, = ax.plot([], [], 'o-', lw=2, color='darkorange')
dots, = ax.plot([], [], 'ko')  # centroid dots
ax.set_xlim(-2, 4)
ax.set_ylim(-2, 4)
ax.set_aspect('equal')
ax.set_title("Animated Recursive Harmonic Spiral of Light")
ax.axis('off')

# Initialization function
def init():
    line.set_data([], [])
    dots.set_data([], [])
    return line, dots

# Animation update function
def update(frame):
    xdata = x[:frame]
    ydata = y[:frame]
    line.set_data(xdata, ydata)

    # Conditionally plot centroids as spiral expands
    visible_points = [(cx, cy) for i, (cx, cy) in enumerate(centroid_points) if frame > (i + 1) * 10]
    if visible_points:
        cx, cy = zip(*visible_points)
        dots.set_data(cx, cy)
    else:
        dots.set_data([], [])
    return line, dots

# Create and display animation
ani = animation.FuncAnimation(fig, update, frames=total_harmonic_states + 1, init_func=init,
                              blit=True, interval=100)
HTML(ani.to_jshtml())
