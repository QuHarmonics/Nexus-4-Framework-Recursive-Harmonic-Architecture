import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Define parameters
zeta_line = 0  # Zeta line at the center
theta = np.linspace(0, 2 * np.pi, 100)  # Angular resolution
z = np.linspace(-10, 10, 100)  # Height along the cones
lens_radius = 1.5  # Radius of the lens

# Quantum Cone (blue)
x_quantum = z * np.cos(theta)
y_quantum = z * np.sin(theta)
z_quantum = z

# Macro Cone (red)
x_macro = -z * np.cos(theta)
y_macro = -z * np.sin(theta)
z_macro = -z

# Lens (yellow disc)
lens_theta, lens_z = np.meshgrid(theta, np.linspace(-lens_radius, lens_radius, 100))
lens_x = np.sqrt(lens_radius**2 - lens_z**2) * np.cos(lens_theta)
lens_y = np.sqrt(lens_radius**2 - lens_z**2) * np.sin(lens_theta)
lens_z = np.sqrt(lens_radius**2 - lens_z**2) * np.sign(lens_z)  # Concave bottom

# Plot
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# Plot Quantum Cone
ax.plot(z, x_quantum, y_quantum, color='blue', label='Quantum Cone (Potential Flow)')

# Plot Macro Cone
ax.plot(z, x_macro, y_macro, color='red', label='Macro Cone (Reality Unfolding)')

# Plot Lens
ax.plot_surface(lens_x, lens_y, zeta_line + lens_z, alpha=0.4, color='yellow')
ax.text(0, 0, zeta_line, "Lens (Compressed Reality)", color='black')

# Zeta Line
ax.plot([0, 0], [0, 0], [-10, 10], color='green', linestyle='--', label='Zeta Line')

# Labels
ax.set_xlabel("X-axis")
ax.set_ylabel("Y-axis")
ax.set_zlabel("Z-axis")
ax.set_title("Dual Cone Simulation with Zeta Line and Lens")
ax.legend()

plt.show()
