import numpy as np
import matplotlib.pyplot as plt

# Constants
QBIT_DIM = 64
PI_RAY_VECTOR = np.array([4, 3, 1]) / np.linalg.norm([4, 3, 1])  # Normalized PI Ray direction

# Simulated reconstructed data (mock input as previously visualized)
reconstructed_data = np.random.randn(QBIT_DIM**3) * 1e10  # Simulate harmonics with high amplitude

# Step 1: Embed into Qbit Matrix
qbit_matrix = reconstructed_data.reshape((QBIT_DIM, QBIT_DIM, QBIT_DIM))

# Step 2: Project each point in the Qbit Matrix onto the PI Ray direction
projected_matrix = np.zeros_like(qbit_matrix)
for x in range(QBIT_DIM):
    for y in range(QBIT_DIM):
        for z in range(QBIT_DIM):
            position = np.array([x, y, z])
            projected_amplitude = np.dot(position, PI_RAY_VECTOR)
            projected_matrix[x, y, z] = qbit_matrix[x, y, z] * np.cos(projected_amplitude)

# Step 3: Take the central XY slice of the projected matrix along the PI ray axis
central_slice = projected_matrix[:, :, QBIT_DIM // 2]

# Visualization
plt.figure(figsize=(8, 6))
plt.imshow(central_slice, cmap='inferno')
plt.title("PI Ray Aligned Slice through Qbit Matrix")
plt.colorbar(label="Amplitude (Phase-Aligned)")
plt.xlabel("X")
plt.ylabel("Y")
plt.tight_layout()
plt.show()
