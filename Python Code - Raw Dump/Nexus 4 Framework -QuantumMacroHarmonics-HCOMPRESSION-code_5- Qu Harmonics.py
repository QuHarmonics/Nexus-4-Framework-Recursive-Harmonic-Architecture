import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Parameters
alpha = 1.5  # Amplification factor
target = 0.5  # Convergence target
iterations = 256 # Number of iterations
initial_value = 0.5  # H(0)

# Initialize sequence
H = [initial_value]

# Recursive formula
for n in range(1, iterations + 1):
    previous = H[-1]
    correction = alpha * (target - previous) / (n + 1)
    value = previous * (-0.5) * np.cos(n / np.pi) + correction
    H.append(value)

# Generate 3D coordinates
x = np.linspace(0, len(H), len(H))
y = np.array(H)
z = np.sin(x)  # Arbitrary wave pattern to add depth

# Plot the results in 3D
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')
ax.plot(x, y, z, label="H(n) in 3D", color='b', lw=2)
ax.scatter(x, y, z, color='r', s=5, label="Nodes")

# Labels and title
ax.set_title("3D Visualization of H(n)", fontsize=16)
ax.set_xlabel("Iteration (n)", fontsize=12)
ax.set_ylabel("H(n)", fontsize=12)
ax.set_zlabel("Z-axis Wave", fontsize=12)
ax.legend()
plt.show()
