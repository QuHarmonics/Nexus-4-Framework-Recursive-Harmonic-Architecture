import matplotlib.pyplot as plt
import numpy as np

# Define parameters
k_bits = np.arange(1, 21)  # First 1 to 20 bits
collision_space_full = 2**256  # Full hash space
collision_space_reduced = 2**k_bits  # Reduced space for k bits

# Collision probabilities
prob_match_k_bits = 1 / collision_space_reduced  # Probability of matching k bits
prob_collision_full = 1 / collision_space_full  # Full collision probability

# Plot
plt.figure(figsize=(12, 6))
plt.plot(k_bits, np.log10(prob_match_k_bits), label="Log10(Probability of Matching First k Bits)", marker='o')
plt.axhline(np.log10(prob_collision_full), color='r', linestyle='--', label="Log10(Full Collision Probability)")

# Formatting the plot
plt.title("Relationship Between Starting Bits and Collision Probability in SHA-256")
plt.xlabel("Number of Starting Bits (k)")
plt.ylabel("Log10(Probability)")
plt.legend()
plt.grid(True)
plt.show()
