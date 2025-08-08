import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
from IPython.display import HTML
import matplotlib as mpl
mpl.rcParams['animation.embed_limit'] = 150  # Set limit to 50 MB instead of 20 MB

# Set up figure
fig, ax = plt.subplots(figsize=(10, 10))

# Parameters
num_nodes = 8
angles = np.linspace(0, 2 * np.pi, num_nodes, endpoint=False)
base_radii = 4 + 0.5 * np.sin(5 * angles)
frozen_indices = [2, 7, 13]

# Node positions
x_nodes = base_radii * np.cos(angles)
y_nodes = base_radii * np.sin(angles)

# Plinko probes
num_probes = 32
probes_x = np.random.uniform(-5, 5, num_probes)
probes_y = np.random.uniform(-5, 5, num_probes)
probes_vx = np.zeros(num_probes)
probes_vy = np.zeros(num_probes)

# Breathing field update
def update(frame):
    ax.clear()
    breathing_radii = base_radii * (1 + 0.05 * np.sin(frame * 0.2))
    x_nodes = breathing_radii * np.cos(angles)
    y_nodes = breathing_radii * np.sin(angles)

    ax.scatter(x_nodes, y_nodes, s=300, edgecolors='black', facecolors='none', linewidths=1.5)
    circle = plt.Circle((0, 0), 5.0, color='lightblue', fill=False, linestyle='dotted', linewidth=1.5)
    ax.add_artist(circle)

    for i in range(num_nodes):
        for j in range(i + 1, num_nodes):
            distance = np.linalg.norm([x_nodes[i] - x_nodes[j], y_nodes[i] - y_nodes[j]])
            if 1.5 < distance < 4.5:
                ax.plot([x_nodes[i], x_nodes[j]], [y_nodes[i], y_nodes[j]], alpha=0.2, linestyle='dashed')

    for i in range(num_probes):
        distances = np.sqrt((x_nodes - probes_x[i])**2 + (y_nodes - probes_y[i])**2)
        nearest_idx = np.argmin(distances)
        dx = x_nodes[nearest_idx] - probes_x[i]
        dy = y_nodes[nearest_idx] - probes_y[i]
        probes_vx[i] += 0.02 * dx
        probes_vy[i] += 0.02 * dy

        probes_x[i] += probes_vx[i]
        probes_y[i] += probes_vy[i]

    ax.scatter(probes_x, probes_y, s=50, color='magenta', alpha=0.7)
    ax.set_aspect('equal')
    ax.set_xlim(-7, 7)
    ax.set_ylim(-7, 7)
    ax.axis('off')
    ax.set_title(" Plinko Breathing Exploration ")

# Create animation
ani = animation.FuncAnimation(fig, update, frames=512, interval=256)

# THIS LINE is critical for Jupyter:
HTML(ani.to_jshtml())
