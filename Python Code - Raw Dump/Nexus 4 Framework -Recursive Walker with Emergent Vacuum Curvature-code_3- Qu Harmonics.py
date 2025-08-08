import numpy as np
import matplotlib.pyplot as plt

# Constants for the recursive collapse spiral
H = 0.35         # Harmonic constant
F = 1.0          # Folding factor
R0 = 1.0         # Initial amplitude
phi = (1 + np.sqrt(5)) / 2  # Golden ratio

# Time steps
n_steps = 16
t_vals = np.array([i * phi for i in range(n_steps)])
R_vals = R0 * np.exp(-H * F * t_vals)

# Spiral coordinates
theta = 2 * np.pi * t_vals / max(t_vals)
x = R_vals * np.cos(theta)
y = R_vals * np.sin(theta)

# Plot
plt.figure(figsize=(8, 8))
plt.plot(x, y, lw=2)
plt.title("Recursive Collapse Spiral\n(Mirror-Corrected Lift Cancellation)", fontsize=14)
plt.axis('equal')
plt.axis('off')
plt.show()
