import numpy as np
import matplotlib.pyplot as plt

# Mark1 Quantum Harmonic Formula for Spiral to Hash
def harmonic_projection(data, threshold, tau=1.0):
    hash_projection = np.zeros(len(data))
    for n in range(len(data)):
        if n < threshold:
            hash_projection[n] = data[n] * (-0.5) * np.cos(n / np.pi)
        else:
            hash_projection[n] = data[n] * (-0.5) * np.cos(n / np.pi) * np.exp(-(n - threshold) / tau)
    return hash_projection

# Reverse Projection: Hash to Spiral
def reverse_projection(hash_projection, threshold, tau=1.0):
    reconstructed = np.zeros(len(hash_projection))
    for n in range(len(hash_projection)-1, -1, -1):  # Reverse loop
        if n >= threshold:
            reconstructed[n] = hash_projection[n] * np.exp((n - threshold) / tau)
        else:
            reconstructed[n] = hash_projection[n]
    return reconstructed

# Generate a Spiral
def generate_spiral(data_length):
    return np.random.randint(0, 2, data_length)  # Random binary data

# Visualization
def visualize_process(original, hash_projection, reconstructed, title):
    plt.figure(figsize=(12, 8))

    # Original Spiral
    plt.subplot(3, 1, 1)
    plt.plot(original, label="Original Spiral", color="blue")
    plt.title(f"{title}: Original Spiral")
    plt.grid()

    # Hash Projection
    plt.subplot(3, 1, 2)
    plt.plot(hash_projection, label="Hash Projection", color="orange")
    plt.title(f"{title}: Hash Projection (Side View)")
    plt.grid()

    # Reconstructed Spiral
    plt.subplot(3, 1, 3)
    plt.plot(reconstructed, label="Reconstructed Spiral", color="green")
    plt.title(f"{title}: Reconstructed Spiral")
    plt.grid()

    plt.tight_layout()
    plt.show()

# Parameters
data_length = 512
threshold = 256
tau = 50.0

# Generate data and project
original_spiral = generate_spiral(data_length)
hash_projection = harmonic_projection(original_spiral, threshold, tau)
reconstructed_spiral = reverse_projection(hash_projection, threshold, tau)

# Visualize the results
visualize_process(original_spiral, hash_projection, reconstructed_spiral, "Mark1 Quantum Harmonic Spiral")
