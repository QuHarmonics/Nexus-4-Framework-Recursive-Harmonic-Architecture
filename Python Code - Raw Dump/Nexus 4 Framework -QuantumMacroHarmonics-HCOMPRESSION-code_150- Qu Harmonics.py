# Lens structure generation with proper dimensions
lens_z = np.linspace(-lens_thickness / 2, lens_thickness / 2, 50)  # Vertical range for lens
theta_lens = np.linspace(0, 2 * np.pi, 500)  # Circular range for lens

# Generate the grid for the lens
lens_z, theta_lens = np.meshgrid(lens_z, theta_lens)

# Define the lens radii based on curvature
lens_x = np.sqrt(lens_curvature_macro - lens_z**2) * np.cos(theta_lens)
lens_y = np.sqrt(lens_curvature_macro - lens_z**2) * np.sin(theta_lens)

# Ensure lens_z is correctly offset to the zeta line
lens_z += zeta_offset

# Plot the updated simulation
fig = plt.figure(figsize=(12, 10))
ax = fig.add_subplot(111, projection='3d')

# Quantum cone
ax.plot(x_quantum, y_quantum, z_quantum, color='blue', label='Quantum Cone (Potential Flow)')

# Macro cone
ax.plot(x_macro, y_macro, z_macro, color='red', label='Macro Cone (Reality Unfolding)')

# Compressed reality disc
ax.plot(x_disc, y_disc, z_disc, color='yellow', label='Compressed Reality (Disc)')

# Base-10 circles
for i, r_base in enumerate(base_radii, start=2):
    x_base = r_base * np.cos(phi)
    y_base = r_base * np.sin(phi)
    ax.plot(x_base, y_base, np.full_like(x_base, zeta_offset), linestyle='--', label=f'Base-{i} Circle (Observer)')

# Lens structure
ax.plot_surface(lens_x, lens_y, lens_z, alpha=0.6, color='green', label='Lens Structure')

# Zeta line
ax.plot([0, 0], [0, 0], [macro_height, quantum_height], color='black', linestyle='--', label='Zeta Line')

# Labels and title
ax.set_title('Dual Cone Simulation with Adjusted Lens and Base-10 Circle')
ax.set_xlabel('X-axis')
ax.set_ylabel('Y-axis')
ax.set_zlabel('Potential (Z-axis)')
ax.legend()

plt.show()
