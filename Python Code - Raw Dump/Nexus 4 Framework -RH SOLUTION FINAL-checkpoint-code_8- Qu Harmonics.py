import numpy as np
import matplotlib.pyplot as plt

# Generate the first few non-trivial Zeta zeros (approximated values)
# These values are commonly known zeros along the critical line Re(s) = 1/2
zeta_zeros = [
    14.134725, 21.022040, 25.010858, 30.424876, 32.935061, 37.586178, 40.918719, 
    43.327073, 48.005150, 49.773832, 52.970321, 56.446247, 59.347044, 60.831778
]

# Convert Zeta zeros into circular projection angles
# Scale using 2 * pi and modulo to fit into a circular framework
angles = [(2 * np.pi * zero) % (2 * np.pi) for zero in zeta_zeros]

# Generate circular coordinates for Zeta zeros
x_coords = [np.cos(angle) for angle in angles]
y_coords = [np.sin(angle) for angle in angles]

# Create the circular projection plot
fig, ax = plt.subplots(figsize=(8, 8))

# Draw the unit circle
circle = plt.Circle((0, 0), 1, color='lightgray', fill=False, linewidth=2)
ax.add_artist(circle)

# Plot Zeta zero projections
ax.scatter(x_coords, y_coords, color='blue', label='Zeta Zeros (Circular Projection)', s=50)

# Annotate points
for i, (x, y) in enumerate(zip(x_coords, y_coords)):
    ax.text(x * 1.1, y * 1.1, f"{i+1}", color='darkblue', fontsize=9)

# Aesthetic adjustments
ax.set_xlim([-1.2, 1.2])
ax.set_ylim([-1.2, 1.2])
ax.axhline(0, color='black', linewidth=0.5, linestyle='--')
ax.axvline(0, color='black', linewidth=0.5, linestyle='--')
ax.set_aspect('equal', adjustable='datalim')
ax.set_title("Circular Projection of Zeta Zeros", fontsize=14)
ax.legend(loc='upper right')

plt.show()
