import numpy as np
import matplotlib.pyplot as plt

# Define constants
GOLDEN_RATIO = 1.618
EXPANSION_RATE = 1.5

# Initial hash input
input_hash = "0753e9bcdebcbe1bbbf19b6bfe6c06bd9606ce89f3c6d57d6d03a6c257f44b3f"
binary_data = np.array([int(input_hash[i:i+2], 16) for i in range(0, len(input_hash), 2)], dtype=np.uint8)

# Step 1: Generate Golden Spiral Path
def generate_golden_spiral(n_points, scale=1.0):
    theta = np.linspace(0, 4 * np.pi, n_points)  # Angular positions
    r = scale * np.exp(GOLDEN_RATIO * theta)  # Radial positions (exponential growth)
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    return x, y

# Step 2: Expand Hash Using Golden Spiral
def expand_hash(binary_data, expansion_rate=EXPANSION_RATE):
    expanded_hash = np.zeros(len(binary_data) * 2, dtype=np.uint8)  # Initialize larger array
    for i in range(len(binary_data)):
        index = int((i * expansion_rate) % len(binary_data))  # Determine complementary index
        expanded_hash[i] = binary_data[i]
        expanded_hash[index] = int((binary_data[i] * expansion_rate) % 256)  # Modulo for wrapping
    return expanded_hash

# Step 3: Visualize Expansion
def visualize_expansion(original_data, expanded_data):
    x_original, y_original = generate_golden_spiral(len(original_data), scale=1.0)
    x_expanded, y_expanded = generate_golden_spiral(len(expanded_data), scale=1.5)
    
    plt.figure(figsize=(12, 6))
    plt.subplot(1, 2, 1)
    plt.scatter(x_original, y_original, c=original_data, cmap='viridis', label="Original Hash")
    plt.title("Original Hash (Golden Spiral)")
    plt.colorbar()
    
    plt.subplot(1, 2, 2)
    plt.scatter(x_expanded, y_expanded, c=expanded_data, cmap='plasma', label="Expanded Hash")
    plt.title("Expanded Hash (Golden Spiral)")
    plt.colorbar()
    
    plt.tight_layout()
    plt.show()

# Run expansion
expanded_hash = expand_hash(binary_data)
visualize_expansion(binary_data, expanded_hash)

# Output results
print("Original Hash (Binary Data):", binary_data)
print("Expanded Hash:", expanded_hash)
