import matplotlib.pyplot as plt
import numpy as np

# Recursive tone parameters
phi = 0.61803398875
H = 0.35
F = 1.0
base_freq = 220.0

# Generate polar (spiral) coordinates
theta_values = []
radius_values = []
labels = []

for i in range(32):
    t = phi * i
    R = np.exp(H * F * t)
    freq = base_freq * R
    theta = 2 * np.pi * phi * i  # angle spreads by golden angle
    radius = np.log(freq)        # use log scale for spiral
    theta_values.append(theta)
    radius_values.append(radius)
    labels.append(f"{round(freq)} Hz")

# Convert to Cartesian coordinates for plotting
x = [r * np.cos(theta) for r, theta in zip(radius_values, theta_values)]
y = [r * np.sin(theta) for r, theta in zip(radius_values, theta_values)]

# Plot
plt.figure(figsize=(8, 8))
plt.plot(x, y, marker='o', linestyle='-', linewidth=2)
for i, label in enumerate(labels):
    plt.text(x[i], y[i], label, fontsize=9, ha='center', va='center')
plt.title("Recursive Harmonic Spiral (Frequency Fractal)", fontsize=14)
plt.axis('equal')
plt.grid(True)
plt.show()
