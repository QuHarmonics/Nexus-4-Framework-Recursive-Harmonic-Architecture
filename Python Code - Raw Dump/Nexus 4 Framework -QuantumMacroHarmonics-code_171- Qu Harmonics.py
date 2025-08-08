import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Constants
golden_ratio = 1.618
termination_ratio = 1 / golden_ratio  # Limit for wave expansion
max_height = 100  # Maximum potential height for the cones
radius_base = 1.0  # Base radius for both cones
twist_maintenance = 0.02  # 2% twist for tension

# Generate quantum cone (potential flow) wave
theta_quantum = np.linspace(0, 4 * np.pi, 500)
z_quantum = np.linspace(0, max_height * termination_ratio, 500)
r_quantum = (radius_base / (max_height * termination_ratio)) * z_quantum
x_quantum = r_quantum * np.cos(theta_quantum)
y_quantum = r_quantum * np.sin(theta_quantum)

# Generate macro cone (reality unfolding) wave
theta_macro = np.linspace(0, 4 * np.pi, 500)
z_macro = np.linspace(0, -max_height * termination_ratio, 500)
r_macro = (radius_base / (max_height * termination_ratio)) * z_macro
x_macro = r_macro * np.cos(-theta_macro)  # Opposite twist
y_macro = r_macro * np.sin(-theta_macro)

# Lens parameters
lens_radius = radius_base
lens_height = max_height * termination_ratio
lens_thickness = 0.05  # Thin lens

# Generate lens
lens_theta = np.linspace(0, 2 * np.pi, 500)
lens_z = np.linspace(-lens_thickness / 2, lens_thickness / 2, 50)
lens_x = np.outer(lens_radius * np.cos(lens_theta), np.ones_like(lens_z))
lens_y = np.outer(lens_radius * np.sin(lens_theta), np.ones_like(lens_z))
lens_z = np.outer(np.ones_like(lens_theta), lens_z)

# Plot the structure
fig = plt.figure(figsize=(12, 10))
ax = fig.add_subplot(111, projection='3d')

# Quantum cone
ax.plot(x_quantum, y_quantum, z_quantum, color='blue', label='Quantum Cone (Potential Flow)')

# Macro cone
ax.plot(x_macro, y_macro, z_macro, color='red', label='Macro Cone (Reality Unfolding)')

# Lens
ax.plot_surface(lens_x, lens_y, lens_z, color='green', alpha=0.5, edgecolor='black', linewidth=0.2, label='Lens Structure')

# Compressed reality disc (event horizon)
phi = np.linspace(0, 2 * np.pi, 500)
x_disc = radius_base * 0.9 * np.cos(phi)  # Compressed inward by 10%
y_disc = radius_base * 0.9 * np.sin(phi)
z_disc = np.zeros_like(x_disc)
ax.plot(x_disc, y_disc, z_disc, color='yellow', label='Base 10 back, Our reality and time')

# Zeta Line
ax.plot([0, 0], [0, 0], [-lens_height, lens_height], color='black', linestyle='--', label='Zeta Line')

# Labels and titles
ax.set_title('Dual Cone Simulation with Twisting Dynamics and Lens Adjustment')
ax.set_xlabel('X-axis')
ax.set_ylabel('Y-axis')
ax.set_zlabel('Potential (Z-axis)')
ax.set_xlim(-1.5, 1.5)
ax.set_ylim(-1.5, 1.5)
ax.set_zlim(-lens_height, lens_height)
ax.legend(loc='upper right')

plt.show()
