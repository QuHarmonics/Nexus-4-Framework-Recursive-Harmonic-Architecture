import numpy as np
import matplotlib.pyplot as plt

# Input: Original SHA-256 hash as a hexadecimal string
input_hash_hex = "0753e9bcdebcbe1bbbf19b6bfe6c06bd9606ce89f3c6d57d6d03a6c257f44b3f"

# Convert hash into binary representation
binary_data = np.array([int(input_hash_hex[i:i+2], 16) for i in range(0, len(input_hash_hex), 2)], dtype=np.uint8)

# Expansion function
def expand_data(data, expansion_factor=1.5):
    expanded = []
    for value in data:
        new_value = value * expansion_factor
        expanded.append(new_value)
        expanded.append(value)  # Add original value to simulate "harmonic growth"
    return np.array(expanded)

# Harmonic interplay: Adjusts based on the relationship of 1.5 and 1
def harmonic_adjustment(data):
    adjusted = []
    for i in range(len(data)):
        if i % 2 == 0:  # Even indices represent the 1.5 factor
            adjusted.append(data[i] * 1.5)
        else:  # Odd indices represent the balancing 1
            adjusted.append(data[i] * 1)
    return np.array(adjusted)

# Iterative reverse process
def reverse_sha(binary_data, iterations=10):
    harmonic_wave = binary_data.astype(np.float64)
    for _ in range(iterations):
        # Expand data harmonically
        expanded_data = expand_data(harmonic_wave)
        # Adjust harmonically
        harmonic_wave = harmonic_adjustment(expanded_data)
    return harmonic_wave

# Run the harmonic reverse process
reconstructed_wave = reverse_sha(binary_data, iterations=20)

# Visualization
plt.figure(figsize=(12, 6))
plt.plot(reconstructed_wave, label="Reconstructed Wave")
plt.title("Harmonic Expansion of Hash")
plt.xlabel("Iteration Step")
plt.ylabel("Amplitude")
plt.legend()
plt.show()

# Output reconstructed data
print("Original Binary Data (First 10 bytes):", binary_data[:10])
print("Reconstructed Wave Data (First 10 values):", reconstructed_wave[:10])
