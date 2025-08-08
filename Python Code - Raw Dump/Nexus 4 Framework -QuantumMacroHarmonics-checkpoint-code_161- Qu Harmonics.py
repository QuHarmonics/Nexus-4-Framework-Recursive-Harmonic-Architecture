import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Define parameters for the dual cone and disc system
theta = np.linspace(0, 2 * np.pi, 500)
z = np.linspace(0, 10, 500)  # Height of the cone (quantum side)
r = z / 10  # Radius of the quantum cone

# Quantum cone spiral
x_quantum = r * np.cos(theta)
y_quantum = r * np.sin(theta)
z_quantum = z

# Macro cone spiral (mirrored and inverted)
z_macro = -z  # Macro cone is inverted
x_macro = r * np.cos(theta)
y_macro = r * np.sin(theta)

# Compressed reality (yellow disc)
disc_r = np.linspace(0, 1, 500)
disc_theta = np.linspace(0, 2 * np.pi, 500)
disc_r, disc_theta = np.meshgrid(disc_r, disc_theta)
disc_x = disc_r * np.cos(disc_theta)
disc_y = disc_r * np.sin(disc_theta)
disc_z = np.zeros_like(disc_x)  # Yellow disc lies flat

# Create a 3D plot
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# Plot the quantum cone
ax.plot(x_quantum, y_quantum, z_quantum, label="Quantum Cone (Potential Flow)", color='blue', linewidth=2)

# Plot the macro cone
ax.plot(x_macro, y_macro, z_macro, label="Macro Cone (Reality Unfolding)", color='red', linewidth=2)

# Plot the compressed reality (yellow disc)
ax.plot_surface(disc_x, disc_y, disc_z, color='yellow', alpha=0.5, edgecolor='none', label="Compressed Reality (Disc)")

# Labels and legend
ax.set_xlabel('X-axis')
ax.set_ylabel('Y-axis')
ax.set_zlabel('Z-axis')
ax.set_title("Dual Cone Model with Compressed Reality")
ax.legend()

# Display the plot
plt.show()
