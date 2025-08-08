import numpy as np
import plotly.graph_objects as go

# Matrix size and direction vector
matrix_size = 64
direction = np.array([4, 3, 1])
direction = direction / np.linalg.norm(direction)

# 3D grid
x, y, z = np.meshgrid(np.arange(matrix_size),
                      np.arange(matrix_size),
                      np.arange(matrix_size),
                      indexing='ij')

# Dot product projection along PI ray
dot = x * direction[0] + y * direction[1] + z * direction[2]
wave = np.cos(2 * np.pi * dot / matrix_size)

# Downsample for visualization
sample = 4
xs = x[::sample, ::sample, ::sample].flatten()
ys = y[::sample, ::sample, ::sample].flatten()
zs = z[::sample, ::sample, ::sample].flatten()
vals = wave[::sample, ::sample, ::sample].flatten()

# 3D scatter plot
fig = go.Figure(data=go.Scatter3d(
    x=xs, y=ys, z=zs,
    mode='markers',
    marker=dict(
        size=3,
        color=vals,
        colorscale='Viridis',
        opacity=0.8,
        colorbar=dict(title='Phase Value')
    )
))

fig.update_layout(
    title="3D Quantum Harmonic Lattice Aligned to PI Ray [4,3,1]",
    scene=dict(
        xaxis_title='X',
        yaxis_title='Y',
        zaxis_title='Z',
        aspectmode='cube'
    ),
    margin=dict(l=0, r=0, t=50, b=0)
)

fig.show()
