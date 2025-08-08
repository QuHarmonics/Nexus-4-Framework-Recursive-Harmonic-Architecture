import numpy as np
import matplotlib.pyplot as plt

# Simulated data layers: each is a 1D array with random values and an embedded pattern
np.random.seed(42)
layers = []
num_layers = 64
length = 64

# Generate random layers and inject a phase-aligned pattern (a sine wave)
base_pattern = np.sin(np.linspace(0, 3 * np.pi, length))
for i in range(num_layers):
    noise = np.random.normal(0, 0.35, length)
    rotation = np.roll(base_pattern, i * 3)  # rotate pattern differently per layer
    layer = noise + rotation
    layers.append(layer)

# Try aligning all layers by rotating them and finding max overlap
def align_layers(layers):
    best_score = -np.inf
    best_shift = 0
    composite = np.zeros(length)

    for shift in range(length):
        aligned_sum = np.zeros(length)
        for layer in layers:
            aligned_sum += np.roll(layer, -shift)
        score = np.sum(aligned_sum ** 2)  # coherence metric
        if score > best_score:
            best_score = score
            best_shift = shift
            composite = aligned_sum

    return best_shift, composite

# Execute alignment
best_shift, composite = align_layers(layers)

# Plot the results
plt.figure(figsize=(12, 6))
plt.title(f"Phase-Aligned Composite at Shift = {best_shift}")
plt.plot(composite, label="Aligned Composite Signal", color='cyan')
plt.axhline(0, color='gray', linestyle='--')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
