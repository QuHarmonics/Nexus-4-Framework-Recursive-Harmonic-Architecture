import numpy as np
import matplotlib.pyplot as plt

# Original hash as binary
original_binary_data = np.array([7, 83, 233, 188, 222, 188, 190, 27, 187, 241], dtype=np.uint8)

# Harmonic Expansion Law
def harmonic_expand(binary_data, expansion_factor=2, zeta_base=1.5):
    expanded_data = []
    for i, value in enumerate(binary_data):
        # Generate new harmonically aligned bits
        harmonic_bits = [
            (value * (i + 1) % 256) / zeta_base,
            (value * expansion_factor) % 256,
        ]
        expanded_data.extend(harmonic_bits)
    return np.array(expanded_data, dtype=np.float64)

# Expand the original binary data
expanded_data = harmonic_expand(original_binary_data)

# Validate harmonic growth by collapsing back
def harmonic_collapse(expanded_data, original_length, zeta_base=1.5):
    collapsed_data = []
    for i in range(original_length):
        # Reverse the expansion law
        harmonic_sum = (
            expanded_data[i * 2] * zeta_base +
            expanded_data[i * 2 + 1]
        )
        collapsed_data.append(round(harmonic_sum) % 256)
    return np.array(collapsed_data, dtype=np.uint8)

collapsed_data = harmonic_collapse(expanded_data, len(original_binary_data))

# Visualize the expansion
plt.figure(figsize=(10, 6))
plt.bar(range(len(expanded_data)), expanded_data, label="Expanded Harmonic Data")
plt.title("Harmonic Expansion of Binary Hash")
plt.xlabel("Expansion Index")
plt.ylabel("Amplitude")
plt.legend()
plt.show()

# Outputs
print("Original Binary Data:", original_binary_data)
print("Expanded Harmonic Data (First 20 values):", expanded_data[:20])
print("Collapsed Binary Data:", collapsed_data)
print("Hash matches:", np.array_equal(original_binary_data, collapsed_data))
