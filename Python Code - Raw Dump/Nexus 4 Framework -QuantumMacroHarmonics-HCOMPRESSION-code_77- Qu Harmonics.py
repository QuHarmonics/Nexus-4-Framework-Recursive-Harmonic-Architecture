import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Constants
EXPANSION_FACTOR = 1.5
QUANTUM_ADJUSTMENT = 0.35

# Generate a quantum wave using Samson principles (from hash)
def quantum_wave_samson(hash_value):
    binary_data = ''.join(format(int(char, 16), '04b') for char in hash_value)
    n = len(binary_data)
    theta = np.linspace(0, 2 * np.pi, n)
    radius = np.linspace(0, 1, n)
    x = radius * np.cos(theta)
    y = radius * np.sin(theta)
    z = np.sin(2 * theta) + np.cos(theta / 2) + QUANTUM_ADJUSTMENT
    return x, y, z

# Generate a wave from seed + padding binary data
def macro_wave_from_binary(binary_data):
    n = len(binary_data)
    theta = np.linspace(0, 2 * np.pi, n)
    radius = np.linspace(0, 1, n)
    x = radius * np.cos(theta)
    y = radius * np.sin(theta)
    z = np.sin(2 * theta) + np.cos(theta / 2) + QUANTUM_ADJUSTMENT
    return x, y, z



# Visualization of two waveforms
def visualize_waves(x1, y1, z1, x2, y2, z2):
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111, projection='3d')

    # Plot Quantum Wave
    ax.plot(x1, y1, z1, label="Quantum Wave (Samson)", color='blue')
    ax.scatter(x1, y1, z1, color='red', s=5)

    # Plot Macro Wave
    ax.plot(x2, y2, z2, label="Macro Wave (Binary)", color='green')
    ax.scatter(x2, y2, z2, color='yellow', s=5)

    ax.set_title("Comparison of Quantum and Macro Waves")
    ax.set_xlabel("X-axis")
    ax.set_ylabel("Y-axis")
    ax.set_zlabel("Z-axis")
    ax.legend()
    plt.show()

# Input hash and seed + padding binary data
hash_value = "3d54e9e342156569b583ac407ffab9ed400568ab718df38fefe19ee3970b02b7"
def sha_padding(seed):
    # Convert seed to binary
    seed_binary = ''.join(format(ord(c), '08b') for c in seed)

    # Append a '1' bit
    seed_binary += '1'

    # Pad with '0' bits until length is 448 (mod 512)
    while len(seed_binary) % 512 != 448:
        seed_binary += '0'

    # Append original message length as a 64-bit big-endian integer
    seed_length = len(seed) * 8
    seed_binary += format(seed_length, '064b')

    return seed_binary

# Input seed (original input to SHA process)
seed = "hello"  # Replace with actual seed

# Generate seed + padding binary
seed_padding_binary = sha_padding(seed)
print("Seed + Padding Binary:", seed_padding_binary)

# Generate Quantum Wave
x_q, y_q, z_q = quantum_wave_samson(hash_value)

# Generate Macro Wave
x_m, y_m, z_m = macro_wave_from_binary(seed_padding_binary)

# Visualize both
visualize_waves(x_q, y_q, z_q, x_m, y_m, z_m)
