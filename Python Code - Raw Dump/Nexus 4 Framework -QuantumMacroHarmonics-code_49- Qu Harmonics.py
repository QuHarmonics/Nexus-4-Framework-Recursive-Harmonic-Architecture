import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Initialize constants and harmonics
def generate_harmonics(binary_data):
    """Generate harmonics from binary data."""
    harmonics = np.cumsum(binary_data.astype(float) * 0.5)  # Scale by factor
    return harmonics

def rotate_observation(harmonics, angle):
    """Rotate the observation of the harmonic container."""
    rotation_matrix = np.array([
        [np.cos(angle), -np.sin(angle)],
        [np.sin(angle), np.cos(angle)]
    ])
    coords = np.stack((np.arange(len(harmonics)), harmonics), axis=-1)
    rotated_coords = coords @ rotation_matrix.T
    return rotated_coords[:, 0], rotated_coords[:, 1]

def calculate_stability_ratio(harmonics):
    """Calculate the stability ratio as ΣPᵢ/ΣAᵢ."""
    positive_energy = np.sum(np.abs(harmonics[harmonics > 0]))
    total_amplitude = np.sum(np.abs(harmonics))
    return positive_energy / total_amplitude if total_amplitude != 0 else 0

# Data and transformation
hash_binary = np.random.choice([0, 1], size=256)  # Replace with actual hash binary data
harmonics = generate_harmonics(hash_binary)

# Observation process with rotational adjustment
angles = np.linspace(0, 2 * np.pi, 100)  # Sweep angles
stability_ratios = []

for angle in angles:
    _, observed_harmonics = rotate_observation(harmonics, angle)
    stability_ratios.append(calculate_stability_ratio(observed_harmonics))

# Find the optimal angle corresponding to the closest stability ratio to 0.35
optimal_angle_index = np.argmin(np.abs(np.array(stability_ratios) - 0.35))
optimal_angle = angles[optimal_angle_index]

# Visualization
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')

x = np.arange(len(harmonics))
z = np.sin(x / 10)  # Add wave depth
_, optimal_observed_harmonics = rotate_observation(harmonics, optimal_angle)

ax.plot(x, optimal_observed_harmonics, z, label="Optimal Observed H(n)", color='green')
ax.set_title(f"Optimal Observation at Angle: {np.degrees(optimal_angle):.2f}°")
ax.set_xlabel("Iteration (n)")
ax.set_ylabel("H(n)")
ax.set_zlabel("Z-axis Wave")
ax.legend()

plt.show()

# Output the decoded data from the optimal observation
decoded_text = ''.join(chr(int(''.join(map(str, hash_binary[i:i+8])), 2)) for i in range(0, len(hash_binary), 8))
print("Decoded Text from Optimal Observation:", decoded_text)
print(f"Stability Ratio at Optimal Angle: {stability_ratios[optimal_angle_index]:.4f}")
