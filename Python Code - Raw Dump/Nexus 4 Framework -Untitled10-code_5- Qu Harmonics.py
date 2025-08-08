
import numpy as np
import matplotlib.pyplot as plt

# Simulated SHA digit stream (as example from prior)
digit_stream = [3, 3, 5, 6, 1, 6, 1, 6, 2, 3, 3, 6, 3, 8, 3, 6, 4, 3, 5, 3]

# Normalize to [0, 1] scale
norm = np.array(digit_stream) / 9.0

# Convert to RGB encoding using phase -> chroma mapping
# R: high phase (>0.6), G: mid (≈0.4–0.6), B: low (<0.3)
rgb = np.zeros((len(norm), 3))
rgb[:, 0] = (norm > 0.6) * norm         # Red
rgb[:, 1] = ((norm > 0.35) & (norm <= 0.6)) * norm  # Green
rgb[:, 2] = (norm <= 0.35) * norm       # Blue

# Visualize as color band
plt.figure(figsize=(12, 1))
plt.imshow([rgb], aspect='auto')
plt.axis('off')
plt.title("Harmonic Lattice Phase Mapping to RGB")
plt.show()

# Optional: print RGB bands numerically
print("RGB Encoded Vector:")
for i, vec in enumerate(rgb):
    print(f"Index {i:2}: R={vec[0]:.2f}, G={vec[1]:.2f}, B={vec[2]:.2f}")

