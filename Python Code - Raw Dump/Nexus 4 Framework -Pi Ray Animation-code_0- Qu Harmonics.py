import matplotlib.pyplot as plt
import numpy as np

# Triangle data: (label, A, B, C, inradius r, circumradius R, centroid (Gx, Gy))
triangles = [
    ("π-ray", np.array([0, 0]), np.array([4, 0]), np.array([1, 0]), 0.0, 1.0, (1.6667, 0.0)),
    ("2-3-4", np.array([0, 0]), np.array([4, 0]), np.array([2.625, 1.45237]), 0.6455, 2.06559, (2.20833, 0.48412)),
    ("2-2-3", np.array([0, 0]), np.array([3, 0]), np.array([1.5, 1.32288]), 0.66144, 1.52059, (1.5, 0.44096)),
    ("3-4-5", np.array([0, 0]), np.array([4, 0]), np.array([4, 3]), 1.0, 2.5, (2.6667, 1.0))
]

# Extract data for plots
labels = []
inradii = []
circumradii = []
centroids_x = []
centroids_y = []

for label, A, B, C, r, R, G in triangles:
    labels.append(label)
    inradii.append(r)
    circumradii.append(R)
    centroids_x.append(G[0])
    centroids_y.append(G[1])

# Plot 1: Inradius and Circumradius over Triangle Sequence
fig1, ax1 = plt.subplots()
x_vals = np.arange(len(triangles))
ax1.plot(x_vals, inradii, 'o-', label='Inradius', linewidth=2)
ax1.plot(x_vals, circumradii, 's--', label='Circumradius', linewidth=2)
ax1.set_title("Recursive Triangle Pulse: Inradius vs Circumradius")
ax1.set_xticks(x_vals)
ax1.set_xticklabels(labels, rotation=45)
ax1.set_ylabel("Radius")
ax1.legend()
ax1.grid(True)

# Plot 2: Centroid Spiral Path
fig2, ax2 = plt.subplots()
ax2.plot(centroids_x, centroids_y, 'o-', linewidth=2)
for i, label in enumerate(labels):
    ax2.text(centroids_x[i], centroids_y[i], label, fontsize=9, ha='right')
ax2.set_title("Centroid Spiral Evolution")
ax2.set_xlabel("X")
ax2.set_ylabel("Y")
ax2.grid(True)
ax2.set_aspect('equal')

plt.tight_layout()
plt.show()
