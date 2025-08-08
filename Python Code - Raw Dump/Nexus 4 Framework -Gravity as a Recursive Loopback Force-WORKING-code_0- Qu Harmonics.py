import numpy as np
import matplotlib.pyplot as plt

# Constants
c = 3e8  # Speed of light (m/s)
hbar = 1.0545718e-34  # Reduced Planck's constant (J·s)
k = 1e-27  # Scaling constant for entanglement equation

# Mass range (kg) from small quantum scale to large macro scale
masses = np.logspace(-30, 5, 500)  # Spanning from quantum to macro

# Traditional General Relativity gravity (Einstein's prediction)
E_mc2 = masses * c**2

# Entangled energy from Nexus2 recursive model
E_entangled = k * np.sum(masses) * (np.sum(masses))**0.35

# Gravity as a loopback force
G_loop = (E_mc2 - E_entangled) * np.exp(-(E_mc2 - E_entangled) / (hbar * c))

# Plot results
plt.figure(figsize=(10, 6))

# Plot Einstein's Gravity (E=mc²)
plt.plot(masses, E_mc2, 'b--', label="Einstein's Gravity (E=mc²)")

# Plot Recursive Gravity Loopback (Nexus2 Model)
plt.plot(masses, G_loop, 'r-', label="Recursive Gravity Loopback (Nexus2)")

# Log scale for better visualization
plt.xscale("log")
plt.yscale("log")

plt.xlabel("Mass (kg)")
plt.ylabel("Energy (J)")
plt.title("Comparison of Einstein's Gravity and Recursive Gravity (Nexus2)")
plt.legend()
plt.grid(True, which="both", linestyle="--", linewidth=0.5)

# Show the plot
plt.show()

# Check if the recursive gravity deviates significantly from Einstein's at high mass levels
mass_threshold = masses[np.argmax(np.abs(E_mc2 - G_loop))]
mass_threshold
