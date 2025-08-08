import numpy as np
import matplotlib.pyplot as plt

# ----- Parameters -----
# Total recursion steps = sum of digits in speed of light (299,792,458)
total_harmonic_states = 55

# Harmonic base (sum of triangle medians from the π-ray field)
harmonic_base = 7.75

# How tightly the spiral curls (angle between steps)
angle_step = 2 * np.pi / harmonic_base

# Distance between steps in spiral
radius_step = 0.15

# ----- Generate Spiral Points -----
theta = np.array([i * angle_step for i in range(total_harmonic_states)])
radii = np.array([i * radius_step for i in range(total_harmonic_states)])
x = radii * np.cos(theta)
y = radii * np.sin(theta)

# ----- Plot the Spiral -----
fig, ax = plt.subplots(figsize=(6, 6))
ax.plot(x, y, marker='o', linestyle='-', color='darkorange', linewidth=2)

# Mark the origin (ZPHCR node)
ax.plot(0, 0, 'ko')  # black dot at center

# Labels and styling
ax.set_title("Recursive Harmonic Spiral of Light (c = 55 steps)", fontsize=14)
ax.set_aspect('equal')
ax.axis('off')
ax.grid(True)

# Show plot
plt.tight_layout()
plt.show()
