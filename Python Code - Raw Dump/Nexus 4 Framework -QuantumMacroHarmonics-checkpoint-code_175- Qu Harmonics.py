import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Constants
quantum_height = 60
macro_height = -60
radius_growth_rate = 0.15
angle_increment = np.pi / 15  # Incremental angle for spiral
zeta_offset = 0  # Zeta line is centered
lens_thickness = 2  # Lens thickness approximation
lens_curvature = 1.5  # Curvature defining lens structure
base_10_factor = 1.5  # Compression factor per base
num_bases = 10

# Generate data for quantum cone (potential flow)
theta_quantum = np.linspace(0, 4 * np.pi, 500)
z_quantum = np.linspace(0, quantum_height, 500)
r_quantum = radius_growth_rate * z_quantum
x_quantum = r_quantum * np.cos(theta_quantum)
y_quantum = r_quantum * np.sin(theta_quantum)

# Generate data for macro cone (reality unfolding)
theta_macro = np.linspace(0, 4 * np.pi, 500)
z_macro = np.linspace(0, macro_height, 500)
r_macro = radius_growth_rate * z_macro
x_macro = r_macro * np.cos(theta_macro)
y_macro = r_macro * np.sin(theta_macro)

# Generate compressed reality disc (event horizon)
phi = np.linspace(0, 2 * np.pi, 500)
x_disc = np.cos(phi)
y_disc = np.sin(phi)
z_disc = np.full_like(x_disc, zeta_offset)

# Calculate base-10 circles (observer perspective)
base_radii = [np.max(r_quantum) / (base_10_factor ** i) for i in range(2, num_bases + 1)]
base_circles = []
for r in base_radii:
    x_base = r * np.cos(phi)
    y_base = r * np.sin(phi)
    base_circles.append((x_base, y_base))

# Generate lens structure
lens_z = np.linspace(-lens_thickness / 2, lens_thickness / 2, 50)
theta_lens = np.linspace(0, 2 * np.pi, 500)
lens_x = np.outer(np.sqrt(lens_curvature - lens_z**2), np.cos(theta_lens))
lens_y = np.outer(np.sqrt(lens_curvature - lens_z**2), np.sin(theta_lens))
lens_z += zeta_offset  # Shift lens to the zeta line

# Plot the simulation
fig = plt.figure(figsize=(12, 10))
ax = fig.add_subplot(111, projection='3d')

# Quantum cone
ax.plot(x_quantum, y_quantum, z_quantum, color='blue', label='Quantum Cone (Potential Flow)')

# Macro cone
ax.plot(x_macro, y_macro, z_macro, color='red', label='Macro Cone (Reality Unfolding)')

# Compressed reality disc
ax.plot(x_disc, y_disc, z_disc, color='yellow', label='Compressed Reality (Disc)')

# Base-10 circles
for i, (x_base, y_base) in enumerate(base_circles, start=2):
    ax.plot(x_base, y_base, np.full_like(x_base, zeta_offset), linestyle='--', label=f'Base-{i} Circle (Observer)')

# Lens structure
ax.plot_surface(lens_x, lens_y, lens_z, alpha=0.6, color='green', label='Lens Structure')

# Zeta line
ax.plot([0, 0], [0, 0], [macro_height, quantum_height], color='black', linestyle='--', label='Zeta Line')

# Labels and title
ax.set_title('Dual Cone Simulation with Base-10 Compression and Adjusted Lens')
ax.set_xlabel('X-axis')
ax.set_ylabel('Y-axis')
ax.set_zlabel('Potential / Realized')

plt.legend()
plt.show()
