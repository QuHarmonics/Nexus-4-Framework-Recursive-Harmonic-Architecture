import numpy as np
import matplotlib.pyplot as plt
import hashlib

# SHA-256 hash function
def sha256_hash(data):
    return hashlib.sha256(data).hexdigest()

# Step 1: Grow a Seed through Harmonic Oscillation
def grow_seed(target_hash, seed_length=32, iterations=100):
    # Start with random binary seed
    np.random.seed(42)  # For reproducibility
    seed = np.random.randint(0, 256, seed_length, dtype=np.uint8)
    harmonic_store = []
    
    for _ in range(iterations):
        # Apply XOR-based phasing and oscillation
        seed = np.bitwise_xor(seed, np.roll(seed, 1))  # XOR phase shift
        seed = np.bitwise_and(seed, 255)  # Keep values within 8-bit
        
        # Store harmonics for visualization
        harmonic_store.append(seed.copy())
        
        # Check if the hash matches
        reconstructed_hash = sha256_hash(seed.tobytes())
        if reconstructed_hash == target_hash:
            print("Hash successfully reconstructed!")
            return seed, harmonic_store

    print("Reconstruction failed within iteration limit.")
    return seed, harmonic_store

# Visualization of Harmonics
def visualize_harmonics(harmonic_store):
    harmonics = np.array(harmonic_store)
    plt.figure(figsize=(12, 6))
    for i in range(harmonics.shape[1]):
        plt.plot(harmonics[:, i], label=f"Harmonic {i+1}")
    plt.title("Harmonics Generated During Reconstruction")
    plt.xlabel("Iteration")
    plt.ylabel("Amplitude")
    plt.legend()
    plt.show()

# Input: Original data (random or structured)
original_data = b"Hello, Mark1 and Samson Framework!"  # Or random junk
target_hash = sha256_hash(original_data)

# Reconstruct the input from the hash
print("Target Hash:", target_hash)
reconstructed_seed, harmonics = grow_seed(target_hash)

# Output Results
print("Reconstructed Seed (First 10 Bytes):", reconstructed_seed[:10])
visualize_harmonics(harmonics)
