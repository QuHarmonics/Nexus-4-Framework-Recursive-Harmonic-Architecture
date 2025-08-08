import numpy as np
import plotly.graph_objects as go

# Constants
HARMONIC_CONSTANT = 0.35  # Scaling factor for harmonics
MAX_ITERATIONS = 100
GAIN_FACTOR = 0.05
GRID_SIZE = 20

# Initialize Lattice with Fixed Size
def initialize_lattice_fixed_size(binary_data, harmonic_constant=HARMONIC_CONSTANT, grid_size=GRID_SIZE):
    normalized_data = binary_data / 255.0  # Normalize binary data to [0, 1]
    lattice = np.zeros((grid_size, grid_size, grid_size), dtype=np.float64)
    
    center = grid_size // 2
    offset = 0
    for idx, value in enumerate(normalized_data):
        x, y, z = (center + offset) % grid_size, (center - offset) % grid_size, (center + offset) % grid_size
        lattice[x, y, z] += value * harmonic_constant
        offset += 1
    return lattice, len(normalized_data)

# Retrieve Data from Lattice
def retrieve_from_lattice(lattice, harmonic_constant=HARMONIC_CONSTANT, data_length=None):
    flattened_data = lattice.flatten() / harmonic_constant  # Scale back by harmonic constant
    binary_data = np.round(flattened_data[:data_length] * 255).astype(np.uint8)  # Crop to original size
    return binary_data

# Apply Feedback Correction
def feedback_correction(lattice, binary_data, harmonic_constant=HARMONIC_CONSTANT):
    retrieved_data = retrieve_from_lattice(lattice, harmonic_constant, data_length=len(binary_data))
    error = (binary_data - retrieved_data) / 255.0  # Normalize the error
    
    for idx, value in enumerate(error):
        x, y, z = idx % lattice.shape[0], (idx // lattice.shape[0]) % lattice.shape[1], idx // (lattice.shape[0] ** 2)
        lattice[x, y, z] += value * harmonic_constant  # Correct the lattice harmonically
    return lattice

# Apply Reflective Gain
def apply_reflective_gain(lattice, gain_factor=GAIN_FACTOR):
    center = lattice.shape[0] // 2
    for x in range(lattice.shape[0]):
        for y in range(lattice.shape[1]):
            for z in range(lattice.shape[2]):
                distance = np.sqrt((x - center)**2 + (y - center)**2 + (z - center)**2)
                lattice[x, y, z] += gain_factor / (1 + distance)
    return lattice

# Generate Interactive Visualization
def generate_visualization(lattice_history):
    frames = []
    for i, lattice in enumerate(lattice_history):
        x, y, z = np.nonzero(lattice)
        values = lattice[x, y, z]
        
        frames.append(
            go.Scatter3d(
                x=x, y=y, z=z, mode='markers',
                marker=dict(size=4, color=values, colorscale='Viridis', opacity=0.8),
                name=f"Iteration {i+1}"
            )
        )
    
    fig = go.Figure(
        data=frames[0],
        layout=go.Layout(
            title="3D Lattice Visualization Over Iterations",
            scene=dict(
                xaxis=dict(title="X-axis"),
                yaxis=dict(title="Y-axis"),
                zaxis=dict(title="Z-axis")
            ),
            updatemenus=[
                dict(
                    buttons=[
                        dict(label=f"Iteration {i+1}",
                             method="restyle",
                             args=["visible", [j == i for j in range(len(frames))]]
                        )
                        for i in range(len(frames))
                    ],
                    direction="down",
                    showactive=True
                )
            ]
        ),
        frames=[go.Frame(data=[frame]) for frame in frames]
    )
    
    fig.show()

# Main Execution
if __name__ == "__main__":
    # Simulated binary data for testing
    binary_data = np.random.randint(0, 256, size=8000, dtype=np.uint8)  # Test dataset

    # Initialize the lattice with fixed size
    lattice, data_length = initialize_lattice_fixed_size(binary_data, grid_size=GRID_SIZE)

    # Store lattice history for visualization
    lattice_history = [lattice.copy()]
    
    for iteration in range(1, MAX_ITERATIONS + 1):
        # Apply feedback correction
        lattice = feedback_correction(lattice, binary_data)
        
        # Apply reflective gain
        lattice = apply_reflective_gain(lattice, gain_factor=GAIN_FACTOR)
        
        # Store the lattice state for this iteration
        lattice_history.append(lattice.copy())
        
        # Retrieve data
        retrieved_data = retrieve_from_lattice(lattice, data_length=data_length)
        
        # Check for data match
        data_matches = np.array_equal(binary_data, retrieved_data)
        if data_matches:
            break

    # Generate the interactive visualization with slider
    generate_visualization(lattice_history)
