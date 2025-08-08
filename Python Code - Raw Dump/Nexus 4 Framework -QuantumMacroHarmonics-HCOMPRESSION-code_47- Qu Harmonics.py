import numpy as np
import matplotlib.pyplot as plt

# Convert hash into binary representation
def hash_to_binary(hash_hex):
    return np.array([int(bit) for bit in ''.join(format(int(c, 16), '04b') for c in hash_hex)])

# Simulate hash unwinding
def unwind_hash(driver_hash, iterations=1000, expansion_factor=1.5):
    unwound_data = np.zeros_like(driver_hash, dtype=np.float64)  # Start with empty data
    expanding_hash = driver_hash.copy()

    for i in range(iterations):
        # Apply expansion rules: "turning the hash inside out"
        expansion = expansion_factor * expanding_hash
        unwound_data += expansion  # Add expanding data outward
        expanding_hash = np.roll(expanding_hash, 1)  # Rotate the driver hash
        unwound_data = np.clip(unwound_data, 0, 1)  # Keep within binary bounds

        # Check for convergence (optional, if unwound_data stabilizes)
        if np.linalg.norm(expansion) < 1e-6:
            print(f"Converged at iteration {i}")
            break

    return unwound_data

# Visualize unwinding process
def plot_unwinding(initial_hash, unwound_data, iterations):
    plt.figure(figsize=(12, 8))
    plt.plot(initial_hash, label="Driver Hash", color='orange')
    plt.plot(unwound_data, label=f"Unwound Data (After {iterations} Iterations)", color='green')
    plt.title("Hash Unwinding Process")
    plt.legend()
    plt.grid()
    plt.show()

# Input hash (SHA-256 of 'abc')
hash_hex = "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
driver_hash = hash_to_binary(hash_hex)

# Perform hash unwinding
iterations = 1000
unwound_data = unwind_hash(driver_hash, iterations)

# Visualize the results
plot_unwinding(driver_hash, unwound_data, iterations)
