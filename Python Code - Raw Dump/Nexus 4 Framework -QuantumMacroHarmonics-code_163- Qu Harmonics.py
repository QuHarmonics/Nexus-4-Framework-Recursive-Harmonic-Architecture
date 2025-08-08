import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Define constants
quantum_height = 10
macro_height = -10
radius_growth_rate = 0.15
angle_increment = np.pi / 15  # Incremental angle for spiral
zeta_offset = 0.5  # Zeta line offset
lens_thickness = 2  # Base-2 thickness
lens_curvature = 1.5  # Curvature defining lens arch
golden_ratio = 1.618  # Approximation of the golden ratio

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

# Generate lens structure
x_lens = np.linspace(-1, 1, 500)
y_lens = np.sqrt(lens_curvature - (x_lens ** 2)) - lens_thickness / 2
z_lens = np.full_like(x_lens, zeta_offset)

# Plot the simulation
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# Quantum cone
ax.plot(x_quantum, y_quantum, z_quantum, color='blue', label='Quantum Cone (Potential Flow)')

# Macro cone
ax.plot(x_macro, y_macro, z_macro, color='red', label='Macro Cone (Reality Unfolding)')

# Compressed reality disc
ax.plot(x_disc, y_disc, z_disc, color='yellow', label='Compressed Reality (Disc)')

# Zeta line
ax.plot([0, 0], [0, 0], [macro_height, quantum_height], color='black', linestyle='--', label='Zeta Line')

# Lens structure
ax.plot(x_lens, y_lens, z_lens, color='green', label='Lens Structure')

# Labels and title
ax.set_title('Dual Cone Simulation with Zeta Line and Lens')
ax.set_xlabel('X-axis')
ax.set_ylabel('Y-axis')
ax.set_zlabel('Z-axis')
ax.legend()

plt.show()
