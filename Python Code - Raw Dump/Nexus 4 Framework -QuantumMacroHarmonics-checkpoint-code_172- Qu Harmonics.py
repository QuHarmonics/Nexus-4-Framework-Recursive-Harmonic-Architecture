import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Define constants
quantum_height = 60
macro_height = -60
base_radius = 1.5  # Shared base radius for spirals and lens
radius_growth_rate = base_radius / quantum_height  # Uniform growth
angle_increment = np.pi / 15  # Incremental angle for spiral
lens_thickness = 0.1  # Adjustable lens thickness
lens_curvature_macro = base_radius  # Base curvature for macro
lens_curvature_quantum = base_radius  # Base curvature for quantum
zeta_offset = 0  # Zeta line is centered

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
x_disc = base_radius * np.cos(phi)
y_disc = base_radius * np.sin(phi)
z_disc = np.full_like(x_disc, zeta_offset)

# Generate lens structure
lens_z = np.linspace(-lens_thickness / 2, lens_thickness / 2, 50)
theta_lens = np.linspace(0, 2 * np.pi, 500)
lens_x = np.outer(np.sqrt(base_radius**2 - lens_z**2), np.cos(theta_lens))
lens_y = np.outer(np.sqrt(base_radius**2 - lens_z**2), np.sin(theta_lens))
lens_z += zeta_offset  # Shift lens to the zeta line

# Create a polarizing effect
polarization_vectors = np.linspace(-1, 1, 50)
polarization_x = polarization_vectors * base_radius
polarization_y = np.zeros_like(polarization_x)
polarization_z = np.zeros_like(polarization_x) + zeta_offset

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
ax.plot_surface(lens_x, lens_y, lens_z, alpha=0.6, color='green', label='Lens Structure')

# Zeta line
ax.plot([0, 0], [0, 0], [macro_height, quantum_height], color='black', linestyle='--', label='Zeta Line')

# Polarization visualization
ax.quiver(
    polarization_x, polarization_y, polarization_z,
    np.zeros_like(polarization_x), np.zeros_like(polarization_x), np.ones_like(polarization_x),
    color="purple", length=5, arrow_length_ratio=0.2, alpha=0.6, label="Polarization"
)

# Labels and title
ax.set_title('Dual Cone Simulation with Twisting Dynamics and Polarization')
ax.set_xlabel('X-axis')
ax.set_ylabel('Y-axis')
ax.set_zlabel('Potential (Z-axis)')
ax.legend()

plt.show()
