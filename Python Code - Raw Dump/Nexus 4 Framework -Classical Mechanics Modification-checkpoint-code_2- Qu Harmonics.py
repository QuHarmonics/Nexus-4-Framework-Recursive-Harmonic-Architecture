import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ─────────────── Physical constants & dummy data ───────────────
m     = 1.0                          # mass (arb units)
v     = 10.0                         # linear speed
c     = 299_792_458.0               # speed of light
γ     = 1.0 / np.sqrt(1 - (v**2)/(c**2))  # Lorentz factor

# Suppose you have 20 measurements at 20 (ω, r) pairs:
real_omega_values = np.linspace(0, 10, 20)
real_r_values     = np.linspace(0.1, 5, 20)

# THEORETICAL KE (old model) for those points:
theoretical_ke_values = 0.5 * m * (v**2 + (real_omega_values**2)*(real_r_values**2))

# EXPERIMENTAL KE (replace this with your real data!)
# Here we simulate ~90% of theory plus a bit of noise:
np.random.seed(0)
real_ke_values = theoretical_ke_values * 0.9 + np.random.normal(0, theoretical_ke_values*0.05)

# ─────────────── Define your updated models ───────────────

def modified_classical_ke(m, v, ω, r):
    linear_ke   = 0.5 * m * v**2
    circular_ke = 0.5 * m * (ω**2) * (r**2)
    return linear_ke + circular_ke

def modified_relativistic_ke(m, γ, c, v, ω, r):
    base_ke     = (γ - 1) * m * c**2
    circular_ke = 0.5 * m * (ω**2) * (r**2)
    return base_ke + circular_ke

def quantum_spin_energy(m, ω, r):
    I = 0.5 * m * (r**2)  # moment of inertia
    return 0.5 * I * (ω**2)

# ─────────────── Compute updated arrays ───────────────

updated_ke_classical     = np.array([
    modified_classical_ke(m, v, ω, r)
    for ω, r in zip(real_omega_values, real_r_values)
])
updated_ke_relativistic  = np.array([
    modified_relativistic_ke(m, γ, c, v, ω, r)
    for ω, r in zip(real_omega_values, real_r_values)
])
quantum_rotational_energy = np.array([
    quantum_spin_energy(m, ω, r)
    for ω, r in zip(real_omega_values, real_r_values)
])

# ─────────────── Plot comparison ───────────────

plt.figure(figsize=(8,6))
plt.plot(real_omega_values, real_ke_values,            'o-', color='red',    label="Experimental KE")
plt.plot(real_omega_values, theoretical_ke_values,     'x--',color='blue',   label="Old Theoretical KE")
plt.plot(real_omega_values, updated_ke_classical,      's-', color='green',  label="Updated Classical KE")
plt.plot(real_omega_values, updated_ke_relativistic,   'd-', color='purple', label="Updated Relativistic KE")
plt.xlabel("ω (Angular Velocity)")
plt.ylabel("Kinetic Energy (arb. units)")
plt.title("Old vs. Updated KE Models vs. Experimental Data")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

# ─────────────── Compute & print MSE ───────────────

mse_updated_ke = np.mean((real_ke_values - updated_ke_classical)**2)
print(f"Mean Squared Error (Experimental vs. Updated Classical KE): {mse_updated_ke:.4e}")

# ─────────────── Build & display DataFrame ───────────────

df_updated_ke = pd.DataFrame({
    "ω":                    real_omega_values,
    "r":                    real_r_values,
    "Experimental KE":      real_ke_values,
    "Old Theoretical KE":   theoretical_ke_values,
    "Updated Classical KE": updated_ke_classical,
    "Updated Relativistic KE": updated_ke_relativistic,
    "Quantum Spin Energy":  quantum_rotational_energy
})

print("\nUpdated KE Models:")
print(df_updated_ke.to_string(index=False))
