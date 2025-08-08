import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

# Grid setup
gridsize = 64
x = np.linspace(0, 4 * np.pi, gridsize)
y = np.linspace(0, 4 * np.pi, gridsize)
x, y = np.meshgrid(x, y)

# Dual wave (input signal + inverse phase)
wave1 = np.sin(x) * np.cos(y)
wave2 = -np.sin(x) * np.cos(y)

# Spiral lens simulation (compression funnel)
def spiral_lens(x, y, strength=0.35):
    r = np.sqrt(x**2 + y**2)
    theta = np.arctan2(y, x)
    z = np.exp(-strength * r) * np.sin(4 * theta)
    return z

lens_field = spiral_lens(x - 2*np.pi, y - 2*np.pi)

# Combined dual wave compression
compressed_wave = (wave1 + wave2) * lens_field

# Orthogonal hash radiation (90 deg escape from Z-axis)
hash_radiation = np.gradient(compressed_wave, axis=0) + np.gradient(compressed_wave, axis=1)

# Plotting
fig = plt.figure(figsize=(14, 6))

# Spiral Compression Field
ax1 = fig.add_subplot(121, projection='3d')
ax1.plot_surface(x, y, compressed_wave, cmap='plasma', edgecolor='none', alpha=0.9)
ax1.set_title('Spiral Harmonic Compression Field')
ax1.set_xlabel('X')
ax1.set_ylabel('Y')
ax1.set_zlabel('Amplitude')

# Orthogonal Hash Radiation
ax2 = fig.add_subplot(122, projection='3d')
ax2.plot_surface(x, y, hash_radiation, cmap='viridis', edgecolor='none', alpha=0.9)
ax2.set_title('Orthogonal Hash Radiation (90° Phase Residue)')
ax2.set_xlabel('X')
ax2.set_ylabel('Y')
ax2.set_zlabel('Radiated Field')

plt.tight_layout()
plt.show()
