import numpy as np
import matplotlib.pyplot as plt

# Define the golden ratio
GOLDEN_RATIO = (1 + 5**0.5) / 2

# Helper function to generate a golden ratio spiral
def golden_ratio_spiral(index):
    return (index * GOLDEN_RATIO) % 1  # Fractional part for normalized values

# Expand the hash harmonically
def harmonic_expand(hash_bytes, iterations=256):
    # Convert hash bytes to integers
    initial_values = np.frombuffer(hash_bytes, dtype=np.uint8)
    expanded_data = np.zeros(iterations, dtype=np.uint64)

    # Initialize the expansion
    previous_value = initial_values.sum()  # Use the sum of the bytes as a seed
    for i in range(iterations):
        spiral_value = int(golden_ratio_spiral(i) * 0xFFFFFFFFFFFFFFFF)  # Scale to 64-bit range
        if i % 2 == 0:
            # Value: XOR the spiral value with the previous value
            expanded_data[i] = previous_value ^ spiral_value
        else:
            # Space: Mirror the previous value
            expanded_data[i] = ~previous_value & 0xFFFFFFFFFFFFFFFF  # Bitwise NOT within 64 bits
        # Update the previous value
        previous_value = expanded_data[i]
    return expanded_data

# Visualize the expanded hash growth
def visualize_growth(expanded_hash):
    plt.figure(figsize=(10, 6))
    plt.plot(expanded_hash, label="Expanded Hash")
    plt.title("Harmonic Expansion of Hash Data")
    plt.xlabel("Iterations")
    plt.ylabel("Values")
    plt.legend()
    plt.show()

# Example hash input (seed hash)
seed_hash = bytes.fromhex("0753e9bcdebcbe1bbbf19b6bfe6c06bd9606ce89f3c6d57d6d03a6c257f44b3f")
seed_bytes = np.frombuffer(seed_hash, dtype=np.uint8)

# Expand the hash
expanded_hash = harmonic_expand(seed_bytes, iterations=256)

# Visualize the expansion
visualize_growth(expanded_hash)
