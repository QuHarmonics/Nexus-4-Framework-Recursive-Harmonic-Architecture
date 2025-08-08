import numpy as np
import matplotlib.pyplot as plt

# Set random seed for reproducibility
np.random.seed(0)

# Simulate bit flipping over time
num_bits = 1000  # Number of bits
num_steps = 100  # Number of time steps
flip_prob = 0.01  # Probability of flipping a bit at each step

# Initialize bits as a binary array (0s and 1s)
bits = np.random.randint(0, 2, size=num_bits)

# Array to store the state of bits over time
bit_history = np.zeros((num_steps, num_bits))

# Simulate bit flipping
for t in range(num_steps):
    # Copy the current state of bits
    bit_history[t, :] = bits
    # Flip bits with probability flip_prob
    flips = np.random.random(num_bits) < flip_prob
    bits[flips] = 1 - bits[flips]  # Flip 0 to 1 or 1 to 0

# Visualize the bit flipping dynamics
plt.figure(figsize=(10, 5))
plt.imshow(bit_history, cmap='binary', interpolation='nearest', aspect='auto')
plt.colorbar(label='Bit Value')
plt.xlabel('Bit Index')
plt.ylabel('Time Step')
plt.title('Bit Flip Dynamics Over Time')
plt.show()