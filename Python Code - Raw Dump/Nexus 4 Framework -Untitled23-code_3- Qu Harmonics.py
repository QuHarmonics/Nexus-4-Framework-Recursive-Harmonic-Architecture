
def modified_classical_ke(m, v, ω, r):
    linear_ke = (1/2) * m * v**2
    circular_ke = (1/2) * m * ω**2 * r**2  # Circular motion energy component
    return linear_ke + circular_ke

# New Relativistic KE Formula incorporating circular motion
def modified_relativistic_ke(m, γ, c, v, ω, r):
    base_ke = (γ - 1) * m * c**2
    circular_ke = (1/2) * m * ω**2 * r**2  # Circular motion correction
    return base_ke + circular_ke

# Compute updated KE values
updated_ke_classical = np.array([modified_classical_ke(m, v, ω, r) for ω, r in zip(real_omega_values, real_r_values)])
updated_ke_relativistic = np.array([modified_relativistic_ke(m, γ, c, v, ω, r) for ω, r in zip(real_omega_values, real_r_values)])

# 2. Investigating potential connections to quantum rotational effects (spin) or frame-dragging
# Using Quantum Angular Momentum Approximation: L = Iω

def quantum_spin_energy(m, ω, r):
    """ Estimates quantum rotational energy using angular momentum principles """
    I = (1/2) * m * r**2  # Moment of inertia for a rotating mass
    return (1/2) * I * ω**2  # Rotational kinetic energy

quantum_rotational_energy = np.array([quantum_spin_energy(m, ω, r) for ω, r in zip(real_omega_values, real_r_values)])

# 3. Comparing old and new KE models
plt.figure(figsize=(8, 6))
plt.plot(real_omega_values, real_ke_values, label="Experimental KE", marker='o', linestyle='-', color='red')
plt.plot(real_omega_values, theoretical_ke_values, label="Theoretical KE (Old)", marker='x', linestyle='--', color='blue')
plt.plot(real_omega_values, updated_ke_classical, label="Updated Classical KE", marker='s', linestyle='-', color='green')
plt.plot(real_omega_values, updated_ke_relativistic, label="Updated Relativistic KE", marker='d', linestyle='-', color='purple')
plt.xlabel("ω (Angular Velocity)")
plt.ylabel("Kinetic Energy (J)")
plt.title("Comparison of Old and New KE Models with Experimental Data")
plt.legend()
plt.grid()
plt.show()

# Compute new MSE to check improvement in fit
mse_updated_ke = np.mean((np.array(real_ke_values) - np.array(updated_ke_classical))**2)

# Display updated KE equations and findings in tabular format
df_updated_ke = pd.DataFrame({
    "ω": real_omega_values,
    "r": real_r_values,
    "Experimental KE": real_ke_values,
    "Old Theoretical KE": theoretical_ke_values,
    "Updated Classical KE": updated_ke_classical,
    "Updated Relativistic KE": updated_ke_relativistic,
    "Quantum Spin Energy": quantum_rotational_energy
})

tools.display_dataframe_to_user(name="Updated Kinetic Energy Models", dataframe=df_updated_ke)

# Display new error metric for accuracy assessment
mse_updated_ke
