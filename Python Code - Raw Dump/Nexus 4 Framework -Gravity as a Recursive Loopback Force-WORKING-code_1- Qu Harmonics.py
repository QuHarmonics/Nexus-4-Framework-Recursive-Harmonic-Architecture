import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Physical constants
G    = 6.67430e-11          # gravitational constant, m^3 kg^-1 s^-2
hbar = 1.054571817e-34      # reduced Planck constant, J·s
c    = 299792458            # speed of light, m/s

# Compute the "critical radius" r_c where ΔE/(ħc) = 1
def critical_radius(m1, m2):
    return np.sqrt(G * m1 * m2 / (hbar * c))

# Mass-pairs to analyze
pairs = [
    ("Proton–Proton",      1.6726219e-27, 1.6726219e-27),
    ("Electron–Electron",  9.10938356e-31, 9.10938356e-31),
    ("Earth–Sun",          5.97219e24,     1.9885e30),
]

# Build a DataFrame of critical radii
rcs = [(name, critical_radius(m1, m2)) for name, m1, m2 in pairs]
df_rc = pd.DataFrame(rcs, columns=["System", "r_c (m)"])

# Print the table
print(df_rc.to_string(index=False))

# Suppression ratio: G_loop/G_newton = exp[-(r_c / r)^2]
def suppression(r, m1, m2):
    rc = critical_radius(m1, m2)
    return np.exp(-(rc / r)**2)

# Plot suppression vs. scaled radius
r_rel = np.logspace(-2, 2, 500)  # r/r_c from 0.01 to 100
plt.figure(figsize=(8,5))
for name, m1, m2 in pairs:
    rc = critical_radius(m1, m2)
    plt.semilogx(r_rel, suppression(r_rel * rc, m1, m2), label=name)

plt.axhline(0.5, color='gray', linestyle='--', label='50% suppression')
plt.xlabel("r / r_c")
plt.ylabel("G_loop / G_newton")
plt.title("Loopback‐Gravity Suppression vs Scaled Radius")
plt.legend()
plt.tight_layout()
plt.show()
