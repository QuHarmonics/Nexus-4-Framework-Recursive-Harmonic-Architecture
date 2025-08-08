import numpy as np
import matplotlib.pyplot as plt

# Sample linear data (e.g. atomic numbers of the first 20 elements)
linear_data =  [1, 4, 1, 5, 9, 2, 6, 5, 3, 5,
                8, 9, 7, 9, 3, 2, 3, 8, 4, 6,
                2, 6, 4, 3, 3, 8, 3, 2, 7, 9,
                5, 0, 2, 8, 8, 4, 1, 9, 7, 1,
                6, 9, 3, 9, 9, 3, 7, 5, 1, 0]

# Calculate differences between consecutive numbers
diffs = np.diff(linear_data)

# Calculate ratios of differences
ratios = diffs / np.max(diffs)

# Function to convert linear data to spiral coordinates
def linear_to_spiral(data, angle_step=np.pi / 6, radius_step=0.2):
    theta = np.arange(len(data)) * angle_step
    r = np.arange(len(data)) * radius_step
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    return x, y

# Generate spiral coordinates for the ratios
x_spiral, y_spiral = linear_to_spiral(ratios)

# Plot the spiral
plt.figure(figsize=(10, 10))
plt.plot(x_spiral, y_spiral, linestyle='--', color='gray')
for i, val in enumerate(ratios):
    plt.plot(x_spiral[i], y_spiral[i], 'o', markersize=10)
    plt.text(x_spiral[i], y_spiral[i], f"{val:.2f}", ha='center', va='center', color='white', fontsize=8)

plt.title("Ratios of Differences in Spiral")
plt.axis("off")
plt.grid(False)
plt.show()