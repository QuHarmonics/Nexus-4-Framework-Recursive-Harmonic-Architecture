# Correct plotting to reattempt visualization
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')

# Quantum cone (blue)
ax.plot3D(x_quantum, y_quantum, z, color='blue', label='Quantum Cone (Potential Flow)')

# Macro cone (red)
ax.plot3D(x_macro, y_macro, z, color='red', label='Macro Cone (Reality Unfolding)')

# Compressed reality disc (yellow)
ax.plot_surface(X_disc, Y_disc, Z_disc, color='yellow', alpha=0.3, edgecolor='none')

# Zeta zeros (green dots)
ax.scatter(zeta_zeros_x, zeta_zeros_y, zeta_zeros_z, color='green', s=50, label='Zeta Zeros')

# Labels and legend
ax.set_title('Dual Cone Model with Zeta Zeros')
ax.set_xlabel('X-axis')
ax.set_ylabel('Y-axis')
ax.set_zlabel('Z-axis')
ax.legend()

plt.show()
