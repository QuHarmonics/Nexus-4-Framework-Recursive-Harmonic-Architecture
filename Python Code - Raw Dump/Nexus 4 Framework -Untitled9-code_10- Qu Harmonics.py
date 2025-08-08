import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
from IPython.display import HTML
from scipy.ndimage import maximum_filter, label, find_objects
import matplotlib as mpl

# Set animation size limit
mpl.rcParams['animation.embed_limit'] = 150

# Set up figure
fig, ax = plt.subplots(figsize=(10, 10))

# Parameters
num_nodes = 18
angles = np.linspace(0, 2 * np.pi, num_nodes, endpoint=False)
base_radii = 4 + 0.35 * np.sin(5 * angles)
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

# Drift strength
drift_strength = 0.035
probe_size = 50

# Memory grid
memory_resolution = 512
memory = np.zeros((memory_resolution, memory_resolution))

# New born nuclei
new_nuclei = []

# Helper: map coordinates to memory grid
def map_to_memory(x, y):
    ix = int((x + 7) / 14 * (memory_resolution - 1))
    iy = int((y + 7) / 14 * (memory_resolution - 1))
    ix = np.clip(ix, 0, memory_resolution - 1)
    iy = np.clip(iy, 0, memory_resolution - 1)
    return ix, iy

# Helper: map memory grid back to coordinates
def memory_to_coords(ix, iy):
    x = (ix / (memory_resolution - 1)) * 14 - 7
    y = (iy / (memory_resolution - 1)) * 14 - 7
    return x, y

# Update function
def update(frame):
    global probes_x, probes_y, probes_vx, probes_vy, memory, new_nuclei
    ax.clear()
    ax.set_facecolor('black')  # Dark background to enhance breathing blooms
    
    # Breathing oscillation
    breathing_radii = base_radii * (1 + 0.05 * np.sin(frame * 0.35))
    x_nodes = breathing_radii * np.cos(angles)
    y_nodes = breathing_radii * np.sin(angles)

    # Plot breathing nodes
    ax.scatter(x_nodes, y_nodes, s=300, edgecolors='white', facecolors='none', linewidths=1.5)
    ax.scatter(x_nodes[frozen_indices], y_nodes[frozen_indices], s=500, facecolors='cyan', edgecolors='white', linewidths=2)

    # Draw breathing echoes
    for i in range(num_nodes):
        for j in range(i + 1, num_nodes):
            distance = np.linalg.norm([x_nodes[i] - x_nodes[j], y_nodes[i] - y_nodes[j]])
            if 1.5 < distance < 3.5:
                ax.plot([x_nodes[i], x_nodes[j]], [y_nodes[i], y_nodes[j]], alpha=0.2, linestyle='dashed', color='white')

    # Memory grid: gentle fade (long breath memory)
    memory *= 0.9995

    # Bloom size around probes
    bloom_size = 3

    # Update Plinko probes and memory blooms
    for i in range(num_probes):
        distances = np.sqrt((x_nodes - probes_x[i])**2 + (y_nodes - probes_y[i])**2)
        nearest_idx = np.argmin(distances)
        dx = x_nodes[nearest_idx] - probes_x[i]
        dy = y_nodes[nearest_idx] - probes_y[i]
        probes_vx[i] += drift_strength * dx
        probes_vy[i] += drift_strength * dy

        probes_x[i] += probes_vx[i]
        probes_y[i] += probes_vy[i]
        
        # Reflect off breathing boundary
        distance_from_center = np.sqrt(probes_x[i]**2 + probes_y[i]**2)
        if distance_from_center > 5.0:
            norm_x = probes_x[i] / distance_from_center
            norm_y = probes_y[i] / distance_from_center
            dot_product = probes_vx[i]*norm_x + probes_vy[i]*norm_y
            probes_vx[i] -= 2 * dot_product * norm_x
            probes_vy[i] -= 2 * dot_product * norm_y
            probes_x[i] = norm_x * 4.9
            probes_y[i] = norm_y * 4.9

        # Bloom memory update
        ix, iy = map_to_memory(probes_x[i], probes_y[i])
        for dx_b in range(-bloom_size, bloom_size + 1):
            for dy_b in range(-bloom_size, bloom_size + 1):
                if 0 <= ix + dx_b < memory_resolution and 0 <= iy + dy_b < memory_resolution:
                    dist = np.sqrt(dx_b**2 + dy_b**2)
                    if dist <= bloom_size:
                        memory[ix + dx_b, iy + dy_b] += (1 - dist / bloom_size)

    # Draw memory bloom
    ax.imshow(np.flipud(memory.T), extent=[-7, 7, -7, 7], cmap='plasma', alpha=0.6)

    # Draw Plinko probes
    ax.scatter(probes_x, probes_y, s=probe_size, color='magenta', alpha=0.7)

    # Birth nuclei when sufficient memory density
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
        ax.scatter(nx, ny, s=400, facecolors='lime', edgecolors='white', linewidths=2, marker='*')

    # Draw breathing boundary
    circle = plt.Circle((0, 0), 5.0, color='lightblue', fill=False, linestyle='dotted', linewidth=1.5)
    ax.add_artist(circle)

    # Settings
    ax.set_aspect('equal')
    ax.set_xlim(-7, 7)
    ax.set_ylim(-7, 7)
    ax.axis('off')
    ax.set_title(" Breathing Collapse Field: Memory Bloom Mode ", color='white')

# Create animation
ani = animation.FuncAnimation(fig, update, frames=512, interval=100)

# Display in Jupyter
HTML(ani.to_jshtml())
