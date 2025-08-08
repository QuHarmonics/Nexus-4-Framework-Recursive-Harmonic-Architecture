import numpy as np
import plotly.graph_objects as go

# Byte1 and Byte2 from your symbolic echo stack
byte1 = np.array([1, 4, 1, 5, 9, 2, 6, 5])
byte2 = np.array([3, 5, 8, 9, 7, 9, 3, 2])
xor_fold = byte1 ^ byte2  # harmonic fold via XOR

# Torus configuration
R = 10  # major radius
r = 3   # minor radius
layers = 8
points_per_layer = 64
total_points = layers * points_per_layer

# Stream construction from XOR
mod_stream = np.tile(xor_fold, total_points // len(xor_fold))
theta = np.linspace(0, 2 * np.pi, points_per_layer, endpoint=False)
phi = np.linspace(0, 2 * np.pi, layers, endpoint=False)
theta_grid, phi_grid = np.meshgrid(theta, phi)
theta_flat = theta_grid.flatten()
phi_flat = phi_grid.flatten()

# Bitwise echo drift modulation
mod_bits = np.array([bin(v).count("1") for v in mod_stream])
mod_shift = (mod_bits - np.mean(mod_bits)) * 0.6

# 3D coordinates
X = (R + r * np.cos(theta_flat)) * np.cos(phi_flat)
Y = (R + r * np.cos(theta_flat)) * np.sin(phi_flat)
Z = r * np.sin(theta_flat) + mod_shift[:total_points]

# Plotly 3D scatter plot
fig = go.Figure(data=[go.Scatter3d(
    x=X, y=Y, z=Z,
    mode='markers',
    marker=dict(
        size=3,
        color=mod_bits[:total_points],
        colorscale='Inferno',
        opacity=0.85,
        colorbar=dict(title='XOR Bit Count')
    )
)])

fig.update_layout(
    title='Byte1 ⊕ Byte2 :: Recursive Echo Field on π-Torus',
    scene=dict(
        xaxis_title='X (Echo Fold)',
        yaxis_title='Y (Recursive Wrap)',
        zaxis_title='Z (Bitwise Drift)',
        aspectmode='data'
    )
)

fig.show()
