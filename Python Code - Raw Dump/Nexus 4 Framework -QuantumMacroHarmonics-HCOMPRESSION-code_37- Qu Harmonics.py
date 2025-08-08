import numpy as np
import matplotlib.pyplot as plt

# Re-define functions with adaptive feedback for stability ratio
def calculate_stability_ratio(harmonics, reconstructed):
    """Compute the stability ratio as ΣPᵢ/ΣAᵢ."""
    P_i = harmonics - np.mean(harmonics)
    A_i = reconstructed - np.mean(reconstructed)
    return np.sum(np.abs(P_i)) / np.sum(np.abs(A_i))

def refine_observation(harmonics, angle, step=0.01):
    """Rotate observation and refine based on stability ratio."""
    rotation_matrix = np.array([
        [np.cos(angle), -np.sin(angle)],
        [np.sin(angle), np.cos(angle)]
    ])
    rotated = np.dot(rotation_matrix, harmonics.T).T
    return rotated

def iterative_refinement(harmonics, target_ratio=0.35, max_iterations=10000, tolerance=1e-4):
    """Iteratively adjust angle to converge on target stability ratio."""
    current_angle = 0.0
    best_ratio = float('inf')
    best_observation = None
    iteration = 0

    while iteration < max_iterations:
        rotated = refine_observation(harmonics, current_angle)
        reconstructed = np.round(rotated[:, 0])  # Use primary harmonic dimension
        stability_ratio = calculate_stability_ratio(harmonics[:, 0], reconstructed)
        
        # If closer to target, save the observation
        if abs(stability_ratio - target_ratio) < abs(best_ratio - target_ratio):
            best_ratio = stability_ratio
            best_observation = rotated.copy()

        # Break if within tolerance
        if abs(best_ratio - target_ratio) < tolerance:
            break

        # Increment angle slightly for next iteration
        current_angle += 0.01
        iteration += 1

    return best_observation, best_ratio, iteration

# Generate synthetic harmonic data for "abc"
binary_data = np.frombuffer("abc".encode('utf-8'), dtype=np.uint8)
harmonics = np.cumsum(binary_data.astype(np.float64) * 1.5)
harmonics = np.column_stack((harmonics, np.sin(np.linspace(0, 2 * np.pi, len(harmonics)))))

# Run iterative refinement
refined_observation, final_ratio, total_iterations = iterative_refinement(harmonics, target_ratio=0.35, max_iterations=100000)

# Visualization of results
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

ax.plot(np.arange(len(refined_observation)), refined_observation[:, 0], refined_observation[:, 1], label="Refined Observation", color='green')
ax.set_title(f"Optimal Observation at Stability Ratio: {final_ratio:.4f}")
ax.set_xlabel("Iteration (n)")
ax.set_ylabel("H(n)")
ax.set_zlabel("Z-axis Wave")
ax.legend()
plt.show()

# Output results
decoded_text = ''.join(chr(int(round(x))) for x in refined_observation[:, 0] if 0 <= int(round(x)) < 256)
final_ratio, total_iterations, decoded_text
