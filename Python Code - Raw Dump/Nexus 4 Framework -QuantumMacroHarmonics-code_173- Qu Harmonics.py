import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Constants
quantum_height = 60
macro_height = -60
radius_growth_rate = 0.15
angle_increment = np.pi / 15  # Incremental angle for spiral
zeta_offset = 0  # Zeta line is centered
lens_thickness = 0.02  # Thickness of the lens (approximation)
lens_curvature_macro = 1.5  # Macro lens curvature
lens_curvature_quantum = 1.618  # Quantum lens curvature (golden ratio)
base10_offset = 10  # Base-10 offset for observer

# Generate data for quantum cone
theta_quantum = np.linspace(0, 4 * np.pi, 500)
z_quantum = np.linspace(0, quantum_height, 500)
r_quantum = radius_growth_rate * z_quantum
x_quantum = r_quantum * np.cos(theta_quantum)
y_quantum = r_quantum * np.sin(theta_quantum)

# Generate data for macro cone
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

# Generate base-10 observer circle
x_base10 = base10_offset * np.cos(phi)
y_base10 = base10_offset * np.sin(phi)
z_base10 = np.full_like(x_base10, zeta_offset)

# Generate lens structure
lens_z = np.linspace(-lens_thickness / 2, lens_thickness / 2, 50)
theta_lens = np.linspace(0, 2 * np.pi, 500)
lens_x = np.outer(np.sqrt(lens_curvature_macro - lens_z**2), np.cos(theta_lens))
lens_y = np.outer(np.sqrt(lens_curvature_macro - lens_z**2), np.sin(theta_lens))
lens_z = np.outer(lens_z, np.ones_like(theta_lens)) + zeta_offset

# Plot the simulation
fig = plt.figure(figsize=(12, 10))
ax = fig.add_subplot(111, projection='3d')

# Quantum cone
ax.plot(x_quantum, y_quantum, z_quantum, color='blue', label='Quantum Cone (Potential Flow)')

# Macro cone
ax.plot(x_macro, y_macro, z_macro, color='red', label='Macro Cone (Reality Unfolding)')

# Compressed reality disc
ax.plot(x_disc, y_disc, z_disc, color='yellow', label='Compressed Reality (Disc)')

# Base-10 observer circle
ax.plot(x_base10, y_base10, z_base10, color='orange', linestyle='--', label='Base-10 Circle (Observer)')

# Lens structure
ax.plot_surface(lens_x, lens_y, lens_z, alpha=0.6, color='green', label='Lens Structure')

# Zeta line
ax.plot([0, 0], [0, 0], [macro_height, quantum_height], color='black', linestyle='--', label='Zeta Line')

# Labels and title
ax.set_title('Dual Cone Simulation with Adjusted Lens and Zeta Line')
ax.set_xlabel('X-axis')
ax.set_ylabel('Y-axis')
ax.set_zlabel('Potential/Realized (Z-axis)')
ax.legend()

plt.show()
