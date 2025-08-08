import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Define constants
quantum_height = 10
macro_height = -10
radius_growth_rate = 0.15
angle_increment = np.pi / 15  # Incremental angle for spiral
zeta_offset = 0  # Zeta line is centered
lens_thickness = 0.02  # 1-bit thickness approximation
lens_curvature_macro = 1.5  # Macro lens curvature
lens_curvature_quantum = 1.618  # Quantum lens curvature (golden ratio)

# Lens radius and deformations
lens_radius = 1  # Base radius of the lens
h_macro = 0.02 * lens_radius  # Macro deformation (2% upward)
h_quantum = 0.02 * lens_radius  # Quantum deformation (2% downward)
lens_height = h_macro + h_quantum

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
x_disc = np.cos(phi) * lens_radius
y_disc = np.sin(phi) * lens_radius
z_disc = np.full_like(x_disc, zeta_offset)

# Generate lens structure (biconvex)
theta_lens = np.linspace(0, 2 * np.pi, 500)
lens_x = np.outer(lens_radius * np.cos(theta_lens), np.ones(50))
lens_y = np.outer(lens_radius * np.sin(theta_lens), np.ones(50))
lens_z = np.linspace(-h_quantum, h_macro, 50)[:, None]  # Concave + Convex

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
ax.plot_surface(lens_x, lens_y, lens_z + zeta_offset, alpha=0.6, color='green', label='Lens Structure')

# Zeta line
ax.plot([0, 0], [0, 0], [macro_height, quantum_height], color='black', linestyle='--', label='Zeta Line')

# Labels and title
ax.set_title('Dual Cone Simulation with Reflected Lens and Calculated Height')
ax.set_xlabel('X-axis')
ax.set_ylabel('Y-axis')
ax.set_zlabel('Z-axis')
ax.legend(loc='upper right')

plt.show()
