import plotly.graph_objects as go
from mpl_toolkits.mplot3d import Axes3D

# Generate 3D toroidal coordinates from the spiral data
def spiral_to_torus(radii, angles, minor_radius=0.3):
    # Torus major and minor radii
    R = 1.5  # major
    r = minor_radius  # minor

    x_3d = (R + radii * np.cos(angles)) * np.cos(angles)
    y_3d = (R + radii * np.cos(angles)) * np.sin(angles)
    z_3d = radii * np.sin(angles)
    return x_3d, y_3d, z_3d

# Generate toroidal coordinates using the same theta_256 and bit_driven_radii_256
x3d, y3d, z3d = spiral_to_torus(bit_driven_radii_256, theta_256)

# Plot 3D torus-based spiral
fig3d = go.Figure(data=go.Scatter3d(
    x=x3d,
    y=y3d,
    z=z3d,
    mode='markers+lines',
    marker=dict(
        size=4,
        color=diffs_256,
        colorscale='Plasma',
        colorbar=dict(title='Δ'),
        opacity=0.8
    ),
    line=dict(color='rgba(120,120,120,0.3)', width=2)
))

fig3d.update_layout(
    title="3D Toroidal Harmonic Spiral (Pi 256 Digits)",
    scene=dict(
        xaxis=dict(visible=False),
        yaxis=dict(visible=False),
        zaxis=dict(visible=False)
    ),
    height=900,
    width=900,
    showlegend=False
)

fig3d.show()
