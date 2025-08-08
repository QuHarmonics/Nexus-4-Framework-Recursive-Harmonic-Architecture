import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.ndimage import gaussian_filter

# Constants
HARMONIC_CONSTANT = 0.35  # Scaling factor for harmonics
MAX_ITERATIONS = 100  # Number of feedback iterations

# Step 1: Encode Data into a 3D Lattice
def store_in_lattice(binary_data, harmonic_constant=HARMONIC_CONSTANT):
    """
    Store binary data into a 3D lattice using harmonic properties.
    """
    normalized_data = binary_data / 255.0  # Normalize binary data to [0, 1]
    
    # Calculate lattice size
    lattice_size = int(np.cbrt(len(normalized_data))) + 1  # Ensure enough space
    total_lattice_elements = lattice_size ** 3
    
    # Pad or truncate the data to match the lattice size
    if len(normalized_data) < total_lattice_elements:
        normalized_data = np.pad(normalized_data, (0, total_lattice_elements - len(normalized_data)))
    elif len(normalized_data) > total_lattice_elements:
        normalized_data = normalized_data[:total_lattice_elements]
    
    # Create and populate the lattice
    lattice = normalized_data.reshape((lattice_size, lattice_size, lattice_size)) * harmonic_constant
    return lattice

# Step 2: Retrieve Data from the 3D Lattice
def retrieve_from_lattice(lattice, harmonic_constant=HARMONIC_CONSTANT):
    """
    Retrieve binary data from a 3D lattice.
    """
    flattened_data = lattice.flatten() / harmonic_constant  # Scale back by harmonic constant
    return np.round(flattened_data * 255).astype(np.uint8)  # Ensure integer byte output

# Step 3: Feedback Correction
def feedback_correction(lattice, binary_data, harmonic_constant=HARMONIC_CONSTANT):
    """
    Apply feedback correction to the lattice to minimize data differences.
    """
    retrieved_data = retrieve_from_lattice(lattice, harmonic_constant)
    # Ensure binary_data and retrieved_data have the same size
    min_size = min(len(binary_data), len(retrieved_data))
    error = (binary_data[:min_size] - retrieved_data[:min_size]) / 255.0  # Normalize the error
    
    for idx, value in enumerate(error):
        x, y, z = idx % lattice.shape[0], (idx // lattice.shape[0]) % lattice.shape[1], idx // (lattice.shape[0] ** 2)
        lattice[x, y, z] += value * harmonic_constant  # Apply feedback correction
    return lattice

# Step 4: Reflective Gain
def apply_reflective_gain(lattice, gain_factor=0.1):
    """
    Apply reflective gain to balance the center of the lattice.
    """
    center = np.array(lattice.shape) // 2
    for x in range(lattice.shape[0]):
        for y in range(lattice.shape[1]):
            for z in range(lattice.shape[2]):
                distance = np.linalg.norm(np.array([x, y, z]) - center)
                scaling = 1 + gain_factor / (1 + distance)
                lattice[x, y, z] *= scaling
    return lattice

# Step 5: Smooth Harmonic Transitions
def smooth_lattice(lattice, sigma=1.0):
    """
    Smooth the lattice using a Gaussian filter.
    """
    return gaussian_filter(lattice, sigma=sigma)

# Visualization of the 3D Lattice
def visualize_lattice(lattice, title="3D Lattice Visualization of Harmonics"):
    """
    Visualize the 3D harmonic lattice.
    """
    x, y, z = np.nonzero(lattice)  # Get non-zero positions
    values = lattice[x, y, z]  # Corresponding harmonic values
    
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    scatter = ax.scatter(x, y, z, c=values, cmap='viridis', s=10)

    # Add labels and title
    ax.set_title(title, fontsize=16)
    ax.set_xlabel("X-axis", fontsize=12)
    ax.set_ylabel("Y-axis", fontsize=12)
    ax.set_zlabel("Z-axis", fontsize=12)
    plt.colorbar(scatter, ax=ax, label="Harmonic Values")
    plt.show()

# Main Execution
if __name__ == "__main__":
    # Load binary bits from your file
    with open(r'd:\\colecovision.rom', 'rb') as file:
        binary_data = np.frombuffer(file.read(), dtype=np.uint8)  # Read binary as bytes

    # Store binary data in harmonic lattice
    lattice = store_in_lattice(binary_data)

    # Iterative correction with reflective gain and smoothing
    for iteration in range(MAX_ITERATIONS):
        print(f"Iteration {iteration + 1}")
        
        # Apply feedback correction
        lattice = feedback_correction(lattice, binary_data)
        
        # Apply reflective gain to center values
        lattice = apply_reflective_gain(lattice, gain_factor=0.1)
        
        # Smooth the lattice
        lattice = smooth_lattice(lattice, sigma=1.0)
        
        # Retrieve and compare data
        retrieved_data = retrieve_from_lattice(lattice)
        visualize_lattice(lattice, title=f"3D Lattice Visualization of Harmonics - Iteration {iteration + 1}")
        
        # Check for convergence
        if np.array_equal(binary_data, retrieved_data):
            print("Data successfully recovered!")
            break

    # Outputs: Compare Original and Retrieved Data
    print("Lattice Shape:", lattice.shape)
    print("Original Data (First 10 Bytes):", binary_data[:10])
    print("Retrieved Data (First 10 Bytes):", retrieved_data[:10])
    
    # Validate if original and retrieved data match
    data_matches = np.array_equal(binary_data, retrieved_data)
    print("Data matches:", data_matches)
    if not data_matches:
        print("Differences detected in data.")
    else:
        print("Data successfully recovered!")
