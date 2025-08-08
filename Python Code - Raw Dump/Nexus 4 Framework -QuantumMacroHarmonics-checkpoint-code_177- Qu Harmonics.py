# Adjusting the simulation based on refined requirements

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Constants
quantum_height = 60
macro_height = -60
radius_growth_rate = 0.15
zeta_offset = 0  # Zeta line is centered
lens_thickness = 2  # Base-2 thickness
lens_base_10_offset = 10  # Base-10 observer ring

# Generate data for Quantum and Macro cones
theta = np.linspace(0, 4 * np.pi, 500)
z_quantum = np.linspace(0, quantum_height, 500)
r_quantum = radius_growth_rate * z_quantum
x_quantum = r_quantum * np.cos(theta)
y_quantum = r_quantum * np.sin(theta)

z_macro = np.linspace(0, macro_height, 500)
r_macro = radius_growth_rate * z_macro
x_macro = r_macro * np.cos(theta)
y_macro = r_macro * np.sin(theta)

# Lens dimensions
lens_radius = r_quantum[-1]  # Same as wave's widest radius
lens_z = np.linspace(-lens_thickness / 2, lens_thickness / 2, 50)
theta_lens = np.linspace(0, 2 * np.pi, 500)
lens_x = np.outer(np.sqrt(lens_radius**2 - (lens_z**2 / (lens_thickness / 2)**2) * lens_radius**2), np.cos(theta_lens))
lens_y = np.outer(np.sqrt(lens_radius**2 - (lens_z**2 / (lens_thickness / 2)**2) * lens_radius**2), np.sin(theta_lens))
lens_z = np.outer(lens_z, np.ones_like(theta_lens))

# Base-2 Circle (outer edge of lens)
phi = np.linspace(0, 2 * np.pi, 500)
x_base_2 = lens_radius * np.cos(phi)
y_base_2 = lens_radius * np.sin(phi)

# Base-10 Circle (inside the lens)
x_base_10 = (lens_radius - lens_base_10_offset) * np.cos(phi)
y_base_10 = (lens_radius - lens_base_10_offset) * np.sin(phi)

# Plotting
fig = plt.figure(figsize=(12, 10))
ax = fig.add_subplot(111, projection='3d')

# Quantum Cone
ax.plot(x_quantum, y_quantum, z_quantum, color='blue', label='Quantum Cone (Potential Flow)')

# Macro Cone
ax.plot(x_macro, y_macro, z_macro, color='red', label='Macro Cone (Reality Unfolding)')

# Lens structure
ax.plot_surface(lens_x, lens_y, lens_z + zeta_offset, alpha=0.4, color='green', label='Lens Structure')

# Base-2 Circle
ax.plot(x_base_2, y_base_2, np.full_like(x_base_2, zeta_offset), color='yellow', linestyle='-', label='Base-2 Circle (Lens Edge)')

# Base-10 Circle
ax.plot(x_base_10, y_base_10, np.full_like(x_base_10, zeta_offset), color='orange', linestyle='--', label='Base-10 Circle (Our Reality)')

# Zeta Line
ax.plot([0, 0], [0, 0], [macro_height, quantum_height], color='black', linestyle='--', label='Zeta Line')

# Labels and title
ax.set_title('Dual Cone Simulation with Lens Dynamically Linked to Wave')
ax.set_xlabel('X-axis')
ax.set_ylabel('Y-axis')
ax.set_zlabel('Potential/Realized Z-axis')
ax.legend()

plt.show()
