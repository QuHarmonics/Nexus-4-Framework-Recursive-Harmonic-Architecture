import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Step 1: Initialize π digits and zeta zeros
pi_digits = [ 4, 1, 5, 9, 2, 6, 5, 3, 5, 8, 9, 7, 9, 3, 2, 3, 8, 4]
zeta_real_parts = [14.134725, 21.02204, 25.010857, 30.424876, 32.935062,
                   37.586178, 40.918719, 43.327073, 48.00515, 49.773832]

# Step 2: Calculate differences and cumulative effects
pi_differences = np.diff(pi_digits)
zeta_differences = np.diff(zeta_real_parts)

# Cumulative sums
pi_cumulative = np.cumsum(pi_digits)
zeta_cumulative = np.cumsum(zeta_real_parts)

# Ratios of differences
pi_ratio_of_differences = np.divide(
    pi_differences[1:].astype(float),
    pi_differences[:-1].astype(float),
    out=np.zeros_like(pi_differences[1:], dtype=float),
    where=pi_differences[:-1] != 0
)

zeta_ratio_of_differences = np.divide(
    zeta_differences[1:].astype(float),
    zeta_differences[:-1].astype(float),
    out=np.zeros_like(zeta_differences[1:], dtype=float),
    where=zeta_differences[:-1] != 0
)

# Combined effects: cumulative + ratio of differences (adjusted shapes)
pi_combined_effect = pi_cumulative[1:-1] * pi_ratio_of_differences
zeta_combined_effect = zeta_cumulative[1:-1] * zeta_ratio_of_differences

# Step 3: Visualize in 3D
fig = plt.figure(figsize=(15, 10))

# π dynamics in 3D
ax1 = fig.add_subplot(121, projection='3d')
ax1.scatter(pi_cumulative[1:-1], pi_ratio_of_differences, pi_combined_effect, c='blue', label='π Dynamics')
ax1.set_xlabel('Cumulative Sum (π)')
ax1.set_ylabel('Ratio of Differences (π)')
ax1.set_zlabel('Combined Effect (π)')
ax1.set_title('3D Visualization of π Dynamics')
ax1.legend()

# Zeta dynamics in 3D
ax2 = fig.add_subplot(122, projection='3d')
ax2.scatter(zeta_cumulative[1:-1], zeta_ratio_of_differences, zeta_combined_effect, c='red', label='Zeta Dynamics')
ax2.set_xlabel('Cumulative Sum (Zeta)')
ax2.set_ylabel('Ratio of Differences (Zeta)')
ax2.set_zlabel('Combined Effect (Zeta)')
ax2.set_title('3D Visualization of Zeta Dynamics')
ax2.legend()

plt.tight_layout()
plt.show()
