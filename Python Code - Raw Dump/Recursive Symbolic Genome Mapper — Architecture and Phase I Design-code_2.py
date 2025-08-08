import numpy as np
import plotly.graph_objects as go

# Full 8-byte structure from your Nexus kernel
bytes_full = np.array([
    [1, 4, 1, 5, 9, 2, 6, 5],    # Byte1
    [3, 5, 8, 9, 7, 9, 3, 2],    # Byte2
    [3, 8, 4, 6, 2, 6, 4, 3],    # Byte3
    [3, 8, 3, 2, 7, 9, 5, 0],    # Byte4
    [2, 8, 8, 4, 1, 9, 7, 1],    # Byte5
    [6, 9, 3, 9, 9, 3, 7, 5],    # Byte6
    [1, 0, 5, 8, 2, 0, 9, 7],    # Byte7
    [4, 5, 9, 2, 3, 0, 7, 8]     # Byte8
])

# XOR each byte with the next in the stack
folds = [bytes_full[i] ^ bytes_full[i + 1] for i in range(len(bytes_full) - 1)]
folds_flat = np.concatenate(folds)

# Toroidal configuration
R = 10
r = 3
layers = len(folds)
points_per_layer = len(folds[0])
total_points = layers * points_per_layer

# Torus angles
theta = np.linspace(0, 2 * np.pi, points_per_layer, endpoint=False)
phi = np.linspace(0, 2 * np.pi, layers, endpoint=False)
theta_grid, phi_grid = np.meshgrid(theta, phi)

theta_flat = theta_grid.flatten()
phi_flat = phi_grid.flatten()

# Bit-count modulation (echo tension)
mod_bits = np.array([bin(val).count("1") for val in folds_flat])
mod_shift = (mod_bits - np.mean(mod_bits)) * 0.6

# Toroidal coordinates
X = (R + r * np.cos(theta_flat)) * np.cos(phi_flat)
Y = (R + r * np.cos(theta_flat)) * np.sin(phi_flat)
Z = r * np.sin(theta_flat) + mod_shift[:total_points]

# Plot using Plotly
fig = go.Figure(data=[go.Scatter3d(
    x=X, y=Y, z=Z,
    mode='markers',
    marker=dict(
        size=3,
        color=mod_bits[:total_points],
        colorscale='Plasma',
        opacity=0.85,
        colorbar=dict(title='XOR Bit Count')
    )
)])

fig.update_layout(
    title='Byte1–8 :: Full Stack Recursive Echo Torus (XOR Fold Lattice)',
    scene=dict(
        xaxis_title='X (Fold Phase)',
        yaxis_title='Y (Layer Rotation)',
        zaxis_title='Z (Echo Tension)',
        aspectmode='data'
    )
)

fig.show()
