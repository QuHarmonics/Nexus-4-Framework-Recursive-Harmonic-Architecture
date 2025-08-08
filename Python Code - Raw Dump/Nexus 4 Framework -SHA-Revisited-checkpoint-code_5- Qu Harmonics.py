import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import plotly.graph_objects as go

# Recursive harmonic folding functions
def fold_data_with_metadata(binary_data, block_size):
    padding_length = (block_size - (len(binary_data) % block_size)) % block_size
    padded_binary = binary_data + '0' * padding_length
    blocks = [padded_binary[i:i + block_size] for i in range(0, len(padded_binary), block_size)]
    folded_hex = ''.join(hex(int(block, 2))[2:].zfill(block_size // 4).upper() for block in blocks)
    metadata = f"{block_size:04X}{padding_length:04X}"
    return folded_hex + metadata

def unfold_data_with_metadata(folded_hex, block_size=None):
    metadata_start = -8
    metadata = folded_hex[metadata_start:]
    folded_hex = folded_hex[:metadata_start]
    block_size = block_size or int(metadata[:4], 16)
    padding_length = int(metadata[4:], 16)
    binary_blocks = [
        bin(int(folded_hex[i:i + block_size // 4], 16))[2:].zfill(block_size)
        for i in range(0, len(folded_hex), block_size // 4)
    ]
    unfolded_binary = ''.join(binary_blocks)
    return unfolded_binary[:len(unfolded_binary) - padding_length]

# Visualization: 2D Harmonic spiral glyph
def plot_folded_spiral(binary_data, title="Harmonic Spiral Fold"):
    steps = len(binary_data)
    angles = np.linspace(0, 4 * np.pi, steps)
    radii = np.linspace(0.5, 1.5, steps)
    x = radii * np.cos(angles)
    y = radii * np.sin(angles)
    colors = ['black' if bit == '1' else 'lightgray' for bit in binary_data]

    plt.figure(figsize=(8, 8))
    plt.scatter(x, y, c=colors, s=20)
    plt.title(title)
    plt.axis('off')
    plt.gca().set_aspect('equal', 'box')
    legend_elements = [
        mpatches.Patch(color='black', label='1'),
        mpatches.Patch(color='lightgray', label='0')
    ]
    plt.legend(handles=legend_elements, loc='upper right')
    plt.show()

# 3D Harmonic spiral using Plotly
def plot_3d_spiral(binary_data, title="3D Harmonic Spiral"):
    steps = len(binary_data)
    theta = np.linspace(0, 6 * np.pi, steps)
    z = np.linspace(0, 1, steps)
    r = np.linspace(0.5, 1.5, steps)
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    color_map = ['black' if bit == '1' else 'lightgray' for bit in binary_data]

    fig = go.Figure(data=[
        go.Scatter3d(
            x=x, y=y, z=z,
            mode='markers',
            marker=dict(
                size=3,
                color=color_map,
            )
        )
    ])
    fig.update_layout(title=title, showlegend=False, scene=dict(xaxis=dict(visible=False), yaxis=dict(visible=False), zaxis=dict(visible=False)))
    fig.show()

# Test the upgraded system
def test_harmonic_fold():
    original_binary = '110100101011001010101001110100101011001010101001' * 3
    block_size = 32

    print("Original Binary:", original_binary)
    folded_hex = fold_data_with_metadata(original_binary, block_size)
    print("\nFolded Hex:", folded_hex)

    unfolded_binary = unfold_data_with_metadata(folded_hex)
    print("\nUnfolded Binary:", unfolded_binary)
    print("\nMatch with Original:", unfolded_binary == original_binary)

    plot_folded_spiral(original_binary, title="Original Binary Spiral")
    plot_folded_spiral(unfolded_binary, title="Unfolded Binary Spiral")
    plot_3d_spiral(original_binary, title="Original Binary 3D Spiral")
    plot_3d_spiral(unfolded_binary, title="Unfolded Binary 3D Spiral")

# Run test
test_harmonic_fold()
