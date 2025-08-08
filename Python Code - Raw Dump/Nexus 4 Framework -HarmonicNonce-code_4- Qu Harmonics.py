import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
from IPython.display import HTML

# Set up the figure
fig, ax = plt.subplots(figsize=(10, 10))

# Field parameters
num_nodes = 18
angles = np.linspace(0, 2 * np.pi, num_nodes, endpoint=False)
base_radii = 4 + 0.5 * np.sin(5 * angles)
frozen_indices = [2, 7, 13]  # Indices of frozen turbulence anchors

# Update function for animation
def update(frame):
    ax.clear()
    
    # Breathing oscillation (sinusoidal modulation)
    breathing_radii = base_radii * (1 + 0.05 * np.sin(frame * 0.2))
    x = breathing_radii * np.cos(angles)
    y = breathing_radii * np.sin(angles)

    # Plot breathing nodes
    ax.scatter(x, y, s=300, edgecolors='black', facecolors='none', linewidths=1.5)

    # Highlight frozen turbulence anchors
    ax.scatter(x[frozen_indices], y[frozen_indices], s=500, facecolors='cyan', edgecolors='black', linewidths=2)

    # Draw breathing echoes
    for i in range(num_nodes):
        for j in range(i + 1, num_nodes):
            distance = np.linalg.norm([x[i] - x[j], y[i] - y[j]])
            if 1.5 < distance < 4.5:
                ax.plot([x[i], x[j]], [y[i], y[j]], alpha=0.2, linestyle='dashed')

    # Draw breathing boundary
    breathing_boundary = plt.Circle((0, 0), 5.0, color='lightblue', fill=False, linestyle='dotted', linewidth=1.5)
    ax.add_artist(breathing_boundary)

    # Plot settings
    ax.set_aspect('equal')
    ax.set_xlim(-7, 7)
    ax.set_ylim(-7, 7)
    ax.axis('off')
    ax.set_title("Breathing Collapse Field")

# Create animation
ani = animation.FuncAnimation(fig, update, frames=120, interval=100)

# Show live in Jupyter Notebook
HTML(ani.to_jshtml())
