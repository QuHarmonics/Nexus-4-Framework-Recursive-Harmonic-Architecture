import numpy as np
import matplotlib.pyplot as plt

# Simulation parameters
timesteps = 500
beta = 0.65  # Reflective memory coefficient
H = 0.35     # Harmonic constant

# Initialize recursive states and existence
R = np.zeros(timesteps)
E = np.zeros(timesteps)

# Set initial condition
R[0] = 1.0  # Initial recursive state

# Recursive loop based on REL formula
for t in range(1, timesteps):
    delta_R = R[t-1] - R[t-2] if t > 1 else R[t-1]
    E[t] = delta_R
    R[t] = E[t] + beta * R[t-1]

# Plotting
plt.figure(figsize=(10, 6))
plt.plot(R, label="Recursive State R(t)", linewidth=2)
plt.plot(E, label="Existence E(t) - Folded Delta", linestyle='--', linewidth=2)
plt.title("Recursive Existence Loop (REL) Simulation")
plt.xlabel("Time Step (t)")
plt.ylabel("Value")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
