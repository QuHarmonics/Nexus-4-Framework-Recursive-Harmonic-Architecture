import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ─────────── Parameter sweep for classical KE ───────────
m = 1.0    # mass (arbitrary units)
v = 10.0   # linear speed (same units)

omega_values = np.linspace(0, 10, 20)
r_values     = np.linspace(0.1, 5, 20)

ke_classical_results = []
for ω in omega_values:
    for r in r_values:
        ke_cl = 0.5 * m * (v**2 + (ω**2)*(r**2))
        ke_classical_results.append([ω, r, ke_cl])

df_ke = pd.DataFrame(ke_classical_results, columns=["ω", "r", "KE_Theoretical"])

# ─────────── Define missing‐energy fraction ───────────
# Suppose you observe only 61% of this KE in your experiments:
mean_ratio = 0.61  
missing_fraction = 1 - mean_ratio  # = 0.39

# Extract arrays for plotting
real_omega_values       = df_ke["ω"].values
real_r_values           = df_ke["r"].values
theoretical_ke_values   = df_ke["KE_Theoretical"].values
missing_energy_values   = theoretical_ke_values * missing_fraction

# ─────────── Geometric‐model functions ───────────
def circular_motion_energy(missing_energy, ω, r):
    return 0.5 * m * (ω**2) * (r**2)

def spiral_motion_energy(missing_energy, ω, r):
    return missing_energy * np.log1p(ω * r)

def vortex_energy(missing_energy, ω, r):
    return missing_energy * (ω**2 / (1 + r**2))

# Compute fits
circular_fit = np.array([
    circular_motion_energy(e, ω, r)
    for e, ω, r in zip(missing_energy_values, real_omega_values, real_r_values)
])
spiral_fit = np.array([
    spiral_motion_energy(e, ω, r)
    for e, ω, r in zip(missing_energy_values, real_omega_values, real_r_values)
])
vortex_fit = np.array([
    vortex_energy(e, ω, r)
    for e, ω, r in zip(missing_energy_values, real_omega_values, real_r_values)
])

# ─────────── Plot observed vs. fits ───────────
plt.figure(figsize=(10,6))
plt.plot(real_omega_values, missing_energy_values, 'o-', color='black', label="Observed Missing Energy")
plt.plot(real_omega_values, circular_fit,      's--', color='red',    label="Circular Motion Fit")
plt.plot(real_omega_values, spiral_fit,        'd--', color='blue',   label="Spiral Motion Fit")
plt.plot(real_omega_values, vortex_fit,        'x--', color='green',  label="Vortex Energy Fit")
plt.xlabel("ω (Angular Velocity)")
plt.ylabel("Missing Energy (arb. units)")
plt.title("Missing Energy (39%) vs. Geometric Models")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

# ─────────── Compute & report MSEs ───────────
mse_circular     = np.mean((missing_energy_values - circular_fit)**2)
mse_spiral       = np.mean((missing_energy_values - spiral_fit)**2)
mse_vortex       = np.mean((missing_energy_values - vortex_fit)**2)

best_fit = min(
    [("Circular", mse_circular), ("Spiral", mse_spiral), ("Vortex", mse_vortex)],
    key=lambda x: x[1]
)

print("MSE Circular Fit:    ", mse_circular)
print("MSE Spiral Fit:      ", mse_spiral)
print("MSE Vortex Fit:      ", mse_vortex)
print("Best Fitting Shape:  ", best_fit[0])
