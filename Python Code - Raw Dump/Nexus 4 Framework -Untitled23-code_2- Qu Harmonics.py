# Investigating the missing 39% kinetic motion

# Define the missing energy fraction
missing_fraction = 1 - mean_ratio  # 39% missing energy

# Compute the missing energy contribution for each data point
missing_energy_values = np.array(theoretical_ke_values) * missing_fraction

# Hypothesis: If the missing energy corresponds to a swirling kinetic motion, it may follow a geometric shape.
# Let's examine circular, spiral, or vortex-like scaling patterns.

# Define a test function for missing energy based on different geometric scaling laws

def circular_motion_energy(missing_energy, ω, r):
    """ Assumes missing energy is stored in a circular motion kinetic term """
    return (1/2) * m * (ω**2 * r**2)

def spiral_motion_energy(missing_energy, ω, r):
    """ Assumes missing energy follows a logarithmic spiral path """
    return missing_energy * np.log1p(ω * r)

def vortex_energy(missing_energy, ω, r):
    """ Assumes missing energy is distributed in a vortex-like pattern """
    return missing_energy * (ω**2 / (1 + r**2))

# Compute missing energy fits
circular_fit = np.array([circular_motion_energy(e, ω, r) for e, ω, r in zip(missing_energy_values, real_omega_values, real_r_values)])
spiral_fit = np.array([spiral_motion_energy(e, ω, r) for e, ω, r in zip(missing_energy_values, real_omega_values, real_r_values)])
vortex_fit = np.array([vortex_energy(e, ω, r) for e, ω, r in zip(missing_energy_values, real_omega_values, real_r_values)])

# Compare fits against missing energy values
plt.figure(figsize=(10, 6))
plt.plot(real_omega_values, missing_energy_values, 'o-', label="Observed Missing Energy", color='black')
plt.plot(real_omega_values, circular_fit, 's--', label="Circular Motion Fit", color='red')
plt.plot(real_omega_values, spiral_fit, 'd--', label="Spiral Motion Fit", color='blue')
plt.plot(real_omega_values, vortex_fit, 'x--', label="Vortex Energy Fit", color='green')
plt.xlabel("ω (Angular Velocity)")
plt.ylabel("Missing Energy (J)")
plt.title("Comparing Missing Energy with Geometric Models")
plt.legend()
plt.grid()
plt.show()

# Compute error metrics (Mean Squared Error) for each model
mse_circular = np.mean((missing_energy_values - circular_fit) ** 2)
mse_spiral = np.mean((missing_energy_values - spiral_fit) ** 2)
mse_vortex = np.mean((missing_energy_values - vortex_fit) ** 2)

# Determine the best-fitting shape by lowest error
best_fit = min(
    [("Circular", mse_circular), ("Spiral", mse_spiral), ("Vortex", mse_vortex)],
    key=lambda x: x[1]
)

# Display results
{
    "MSE Circular Fit": mse_circular,
    "MSE Spiral Fit": mse_spiral,
    "MSE Vortex Fit": mse_vortex,
    "Best Fitting Shape": best_fit[0]
}
