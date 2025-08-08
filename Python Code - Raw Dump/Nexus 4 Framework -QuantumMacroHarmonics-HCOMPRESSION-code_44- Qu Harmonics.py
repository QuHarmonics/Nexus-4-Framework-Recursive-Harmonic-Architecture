import numpy as np
import matplotlib.pyplot as plt

# Convert hash into quantum state (waveform)
def hash_to_waveform(hash_hex):
    binary_data = np.array([int(bit) for bit in ''.join(format(int(c, 16), '04b') for c in hash_hex)])
    x = np.linspace(0, 2 * np.pi, len(binary_data))
    return np.sin(x) * (binary_data - 0.5)

# Recursive feedback expansion
def quantum_expand(waveform, iterations=1000, alpha=1.5, target_ratio=0.35, tolerance=1e-6):
    expanded_waveform = waveform.copy()
    for i in range(iterations):
        # Calculate current stability ratio
        positive_sum = np.sum(expanded_waveform[expanded_waveform > 0])
        total_sum = np.sum(np.abs(expanded_waveform))
        stability_ratio = positive_sum / total_sum if total_sum != 0 else 0

        # Check convergence
        if abs(stability_ratio - target_ratio) < tolerance:
            print(f"Converged at iteration {i} with stability ratio: {stability_ratio:.6f}")
            break

        # Expand using recursive feedback
        correction = alpha * (target_ratio - stability_ratio)
        expanded_waveform += correction * np.sin(np.linspace(0, 2 * np.pi, len(expanded_waveform)))

    return expanded_waveform

# Visualization
def plot_waveforms(original, expanded, iterations):
    plt.figure(figsize=(12, 6))
    plt.plot(original, label="Original Hash Waveform", color='orange')
    plt.plot(expanded, label=f"Expanded Waveform (After {iterations} Iterations)", color='green')
    plt.title("Quantum Expansion of Hash Waveform")
    plt.legend()
    plt.grid()
    plt.show()

# Input hash (e.g., SHA-256 of 'abc')
hash_hex = "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"

# Convert hash to waveform
hash_waveform = hash_to_waveform(hash_hex)

# Perform quantum expansion
expanded_waveform = quantum_expand(hash_waveform)

# Visualize results
plot_waveforms(hash_waveform, expanded_waveform, iterations=1000)
