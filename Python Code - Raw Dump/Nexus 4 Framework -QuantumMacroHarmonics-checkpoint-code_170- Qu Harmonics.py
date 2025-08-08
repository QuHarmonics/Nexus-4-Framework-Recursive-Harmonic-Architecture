import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Define constants
quantum_height = 100
macro_height = -100
radius_base = 1.0  # Base radius of the cones and lens
yellow_disc_ratio = 0.9  # 10% inward compression for the yellow disc
lens_thickness = 0.05  # Refined thickness
lens_curvature = 1.5  # Uniform curvature for lens

# Generate data for quantum cone (potential flow)
theta_quantum = np.linspace(0, 4 * np.pi, 500)
z_quantum = np.linspace(0, quantum_height, 500)
r_quantum = (radius_base / quantum_height) * z_quantum
x_quantum = r_quantum * np.cos(theta_quantum)
y_quantum = r_quantum * np.sin(theta_quantum)

# Generate data for macro cone (reality unfolding)
theta_macro = np.linspace(0, 4 * np.pi, 500)
z_macro = np.linspace(0, macro_height, 500)
r_macro = (radius_base / abs(macro_height)) * z_macro
x_macro = r_macro * np.cos(theta_macro)
y_macro = r_macro * np.sin(theta_macro)

# Generate compressed reality disc (event horizon)
phi = np.linspace(0, 2 * np.pi, 500)
x_disc = yellow_disc_ratio * np.cos(phi) * radius_base
y_disc = yellow_disc_ratio * np.sin(phi) * radius_base
z_disc = np.zeros_like(x_disc)

# Generate lens structure
lens_z = np.linspace(-lens_thickness / 2, lens_thickness / 2, 50)
theta_lens = np.linspace(0, 2 * np.pi, 500)
lens_x = np.outer(np.cos(theta_lens), np.full_like(lens_z, radius_base))
lens_y = np.outer(np.sin(theta_lens), np.full_like(lens_z, radius_base))
lens_z = np.outer(np.ones_like(theta_lens), lens_z)

# Plot the simulation
fig = plt.figure(figsize=(12, 10))
ax = fig.add_subplot(111, projection='3d')

# Quantum cone
ax.plot(x_quantum, y_quantum, z_quantum, color='blue', label='Quantum Cone (Potential Flow)')

# Macro cone
ax.plot(x_macro, y_macro, z_macro, color='red', label='Macro Cone (Reality Unfolding)')

# Compressed reality disc
ax.plot(x_disc, y_disc, z_disc, color='yellow', label='Compressed Reality (Disc)')

# Lens structure
ax.plot_surface(
    lens_x, lens_y, lens_z,
    alpha=0.5,
    color='green',
    edgecolor='black',
    linewidth=0.5,
    label='Lens Structure'
)

# Zeta line
ax.plot([0, 0], [0, 0], [macro_height, quantum_height], color='black', linestyle='--', label='Zeta Line')

# Labels and title
ax.set_title('Dual Cone Simulation with Adjusted Lens and Base-10 Compression')
ax.set_xlabel('X-axis')
ax.set_ylabel('Y-axis')
ax.set_zlabel('Potential (Z-axis)')
ax.set_xlim(-1.5, 1.5)
ax.set_ylim(-1.5, 1.5)
ax.set_zlim(macro_height, quantum_height)
ax.legend(loc='upper right')

plt.show()
