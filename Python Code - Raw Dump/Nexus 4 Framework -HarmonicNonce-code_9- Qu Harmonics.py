import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
from IPython.display import HTML
from scipy.ndimage import maximum_filter, label, find_objects

# Set up figure
fig, ax = plt.subplots(figsize=(10, 10))

# Parameters
num_nodes = 18
memory_resolution = 200
memory = np.zeros((memory_resolution, memory_resolution))

# Plinko probes
num_probes = 32
probes_x = np.random.uniform(-4.5, 4.5, num_probes)
probes_y = np.random.uniform(-4.5, 4.5, num_probes)
probes_vx = np.zeros(num_probes)
probes_vy = np.zeros(num_probes)

# Drift strength
drift_strength = 0.035
probe_size = 50

# New born nuclei
new_nuclei = []

# Helper: map coordinates to memory grid
def map_to_memory(x, y):
    ix = int((x + 7) / 14 * (memory_resolution - 1))
    iy = int((y + 7) / 14 * (memory_resolution - 1))
    ix = np.clip(ix, 0, memory_resolution - 1)
    iy = np.clip(iy, 0, memory_resolution - 1)
    return ix, iy

# Helper: map grid back to coordinates
def memory_to_coords(ix, iy):
    x = (ix / (memory_resolution - 1)) * 14 - 7
    y = (iy / (memory_resolution - 1)) * 14 - 7
    return x, y

# Multi-breathing node generation
def breathing_nodes_multi(frame):
    angles = np.linspace(0, 2 * np.pi, num_nodes, endpoint=False)
    r1 = 4 + 0.5 * np.sin(5 * angles + 0.05 * frame)
    r2 = 4.2 + 0.3 * np.sin(7 * angles + 0.03 * frame)
    x1 = r1 * np.cos(angles)
    y1 = r1 * np.sin(angles)
    x2 = r2 * np.cos(angles + np.pi/num_nodes)
    y2 = r2 * np.sin(angles + np.pi/num_nodes)
    return (x1, y1), (x2, y2)

# Spiral seeds
def spiral_seeds(n_seeds, t):
    phi = np.pi * (3 - np.sqrt(5))  # golden angle
    indices = np.arange(0, n_seeds)
    r = np.sqrt(indices) * 0.5 + 1.5 + 0.05 * np.sin(0.05 * t)
    theta = indices * phi
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    return x, y

# Update function
def update(frame):
    global probes_x, probes_y, probes_vx, probes_vy, memory, new_nuclei
    ax.clear()
    
    (x1, y1), (x2, y2) = breathing_nodes_multi(frame)

    # Plot first breathing field
    ax.scatter(x1, y1, s=300, edgecolors='black', facecolors='none', linewidths=1.5)
    # Plot second breathing field
    ax.scatter(x2, y2, s=300, edgecolors='cyan', facecolors='none', linewidths=1.5)

    # Spiral breathing seeds
    sx, sy = spiral_seeds(20, frame)
    ax.scatter(sx, sy, s=150, facecolors='lime', edgecolors='black', marker='*', alpha=0.7)

    # Update Plinko probes
    for i in range(num_probes):
        # Drift toward both breathing fields
        distances1 = np.sqrt((x1 - probes_x[i])**2 + (y1 - probes_y[i])**2)
        distances2 = np.sqrt((x2 - probes_x[i])**2 + (y2 - probes_y[i])**2)
        nearest_idx1 = np.argmin(distances1)
        nearest_idx2 = np.argmin(distances2)

        dx = (x1[nearest_idx1] + x2[nearest_idx2]) / 2 - probes_x[i]
        dy = (y1[nearest_idx1] + y2[nearest_idx2]) / 2 - probes_y[i]
        probes_vx[i] += drift_strength * dx
        probes_vy[i] += drift_strength * dy

        probes_x[i] += probes_vx[i]
        probes_y[i] += probes_vy[i]
        
        # Reflect off boundary
        distance_from_center = np.sqrt(probes_x[i]**2 + probes_y[i]**2)
        if distance_from_center > 5.0:
            norm_x = probes_x[i] / distance_from_center
            norm_y = probes_y[i] / distance_from_center
            dot_product = probes_vx[i]*norm_x + probes_vy[i]*norm_y
            probes_vx[i] -= 2 * dot_product * norm_x
            probes_vy[i] -= 2 * dot_product * norm_y
            probes_x[i] = norm_x * 4.9
            probes_y[i] = norm_y * 4.9

        # Update memory
        ix, iy = map_to_memory(probes_x[i], probes_y[i])
        memory[ix, iy] += 1

    # Draw memory heatmap
    ax.imshow(np.flipud(memory.T), extent=[-7, 7, -7, 7], cmap='hot', alpha=0.5)

    # Draw Plinko probes
    ax.scatter(probes_x, probes_y, s=probe_size, color='magenta', alpha=0.7)

    # Birth nuclei
    if frame % 50 == 0 and frame > 0:
        maxima = maximum_filter(memory, size=10) == memory
        labeled, num_objects = label(maxima)
        slices = find_objects(labeled)
        for dy, dx in slices:
            ix = (dx.start + dx.stop) // 2
            iy = (dy.start + dy.stop) // 2
            if memory[ix, iy] > 30:
                x_n, y_n = memory_to_coords(ix, iy)
                if np.sqrt(x_n**2 + y_n**2) < 5.0:
                    new_nuclei.append((x_n, y_n))

    # Draw born nuclei
    for (nx, ny) in new_nuclei:
        ax.scatter(nx, ny, s=400, facecolors='white', edgecolors='black', linewidths=2, marker='*')

    # Breathing boundary
    circle = plt.Circle((0, 0), 5.0, color='lightblue', fill=False, linestyle='dotted', linewidth=1.5)
    ax.add_artist(circle)

    # Settings
    ax.set_aspect('equal')
    ax.set_xlim(-7, 7)
    ax.set_ylim(-7, 7)
    ax.axis('off')
    ax.set_title(" Harmonic Genesis Breathing Field ")

# Create animation
ani = animation.FuncAnimation(fig, update, frames=512, interval=100)

# Display in Jupyter
HTML(ani.to_jshtml())
