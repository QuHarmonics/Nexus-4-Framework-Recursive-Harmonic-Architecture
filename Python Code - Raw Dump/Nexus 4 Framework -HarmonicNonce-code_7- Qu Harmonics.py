import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
from IPython.display import HTML

# Set up figure
fig, ax = plt.subplots(figsize=(10, 10))

# Parameters
num_nodes = 18
angles = np.linspace(0, 2 * np.pi, num_nodes, endpoint=False)
base_radii = 4 + 0.5 * np.sin(5 * angles)
frozen_indices = [2, 7, 13]

# Node positions
x_nodes = base_radii * np.cos(angles)
y_nodes = base_radii * np.sin(angles)

# Plinko probes
num_probes = 32
probes_x = np.random.uniform(-4.5, 4.5, num_probes)
probes_y = np.random.uniform(-4.5, 4.5, num_probes)
probes_vx = np.zeros(num_probes)
probes_vy = np.zeros(num_probes)

# Drift strength
drift_strength = 0.035
probe_size = 50

# Initialize memory grid (heatmap)
memory_resolution = 200
memory = np.zeros((memory_resolution, memory_resolution))

# Helper: map coordinates to memory grid
def map_to_memory(x, y):
    ix = int((x + 7) / 14 * (memory_resolution - 1))
    iy = int((y + 7) / 14 * (memory_resolution - 1))
    ix = np.clip(ix, 0, memory_resolution - 1)
    iy = np.clip(iy, 0, memory_resolution - 1)
    return ix, iy

# Update function
def update(frame):
    global probes_x, probes_y, probes_vx, probes_vy, memory
    ax.clear()
    
    # Breathing oscillation
    breathing_radii = base_radii * (1 + 0.05 * np.sin(frame * 0.35))
    x_nodes = breathing_radii * np.cos(angles)
    y_nodes = breathing_radii * np.sin(angles)

    # Plot breathing nodes
    ax.scatter(x_nodes, y_nodes, s=300, edgecolors='black', facecolors='none', linewidths=1.5)
    ax.scatter(x_nodes[frozen_indices], y_nodes[frozen_indices], s=500, facecolors='cyan', edgecolors='black', linewidths=2)

    # Draw breathing echoes
    for i in range(num_nodes):
        for j in range(i + 1, num_nodes):
            distance = np.linalg.norm([x_nodes[i] - x_nodes[j], y_nodes[i] - y_nodes[j]])
            if 1.5 < distance < 3.5:
                ax.plot([x_nodes[i], x_nodes[j]], [y_nodes[i], y_nodes[j]], alpha=0.2, linestyle='dashed')

    # Update Plinko probes
    for i in range(num_probes):
        distances = np.sqrt((x_nodes - probes_x[i])**2 + (y_nodes - probes_y[i])**2)
        nearest_idx = np.argmin(distances)
        dx = x_nodes[nearest_idx] - probes_x[i]
        dy = y_nodes[nearest_idx] - probes_y[i]
        probes_vx[i] += drift_strength * dx
        probes_vy[i] += drift_strength * dy

        probes_x[i] += probes_vx[i]
        probes_y[i] += probes_vy[i]
        
        # Reflect off the breathing boundary
        distance_from_center = np.sqrt(probes_x[i]**2 + probes_y[i]**2)
        if distance_from_center > 5.0:
            norm_x = probes_x[i] / distance_from_center
            norm_y = probes_y[i] / distance_from_center
            dot_product = probes_vx[i]*norm_x + probes_vy[i]*norm_y
            probes_vx[i] -= 2 * dot_product * norm_x
            probes_vy[i] -= 2 * dot_product * norm_y
            probes_x[i] = norm_x * 4.9
            probes_y[i] = norm_y * 4.9

        # Update memory grid
        ix, iy = map_to_memory(probes_x[i], probes_y[i])
        memory[ix, iy] += 1  # Breathe trust into location

    # Draw memory heatmap (trust accumulation)
    ax.imshow(np.flipud(memory.T), extent=[-7, 7, -7, 7], cmap='hot', alpha=0.5)

    # Draw Plinko probes
    ax.scatter(probes_x, probes_y, s=probe_size, color='magenta', alpha=0.7)

    # Breathing boundary
    circle = plt.Circle((0, 0), 5.0, color='lightblue', fill=False, linestyle='dotted', linewidth=1.5)
    ax.add_artist(circle)

    # Settings
    ax.set_aspect('equal')
    ax.set_xlim(-7, 7)
    ax.set_ylim(-7, 7)
    ax.axis('off')
    ax.set_title("Contained Plinko Breathing Exploration + Cluster Memory")

# Create animation
ani = animation.FuncAnimation(fig, update, frames=512, interval=100)

# Display in Jupyter
HTML(ani.to_jshtml())
