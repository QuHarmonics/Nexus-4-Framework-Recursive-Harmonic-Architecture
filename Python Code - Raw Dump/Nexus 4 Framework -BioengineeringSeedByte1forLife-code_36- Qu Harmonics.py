import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

# Generate a simulated 3D fold representation for the peptide
# Hydrophobic residues (A²) and charged residues (B²) are mapped in different regions

# Simulate a peptide folding path
np.random.seed(42)
x = np.cumsum(np.random.randn(20))
y = np.cumsum(np.random.randn(20))
z = np.cumsum(np.random.randn(20))

# Define residues as hydrophobic (A²) or charged/polar (B²)
residues = [
    'P', 'Q', 'A', 'W', 'P', 'I', 'A', 'K', 'K', 'K', 'K',
    'S', 'C', 'E', 'C', 'G', 'D', 'R', 'R', 'A'
]
colors = ['red' if r in ['K', 'R', 'E', 'C', 'S', 'D'] else 'blue' for r in residues]

# Plot the folding representation
fig = plt.figure(figsize=(12, 6))
ax = fig.add_subplot(111, projection='3d')
scatter = ax.scatter(x, y, z, c=colors, s=100, alpha=0.8)
ax.plot(x, y, z, color='gray', alpha=0.5)
for i, txt in enumerate(residues):
    ax.text(x[i], y[i], z[i], txt, fontsize=9)

ax.set_title("3D Representation of Peptide Folding")
ax.set_xlabel("X-axis")
ax.set_ylabel("Y-axis")
ax.set_zlabel("Z-axis")

plt.show()
