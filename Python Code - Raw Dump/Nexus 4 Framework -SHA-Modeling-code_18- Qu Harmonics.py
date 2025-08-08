# All-in-one working code with example usage
import numpy as np
from mpmath import mp
import plotly.graph_objects as go

# Function to generate spiral and plot it based on number of Pi digits
def generate_pi_spiral(num_digits=256):
    # Set precision based on requested digits
    mp.dps = num_digits + 2
    pi_str = str(mp.pi)[2:num_digits+2]

    # Convert to integer list
    pi_digits = [int(d) for d in pi_str]
    diffs = np.diff(pi_digits)

    # Bit length as harmonic driver
    bit_lengths = [len(bin(d)[2:]) for d in pi_digits[:-1]]
    bit_driven_radii = np.array(bit_lengths) * 0.2

    # Golden ratio angle shift
    phi = (1 + np.sqrt(5)) / 2
    theta = np.arange(len(diffs)) * (2 * np.pi / phi)

    # Toroidal transformation
    def spiral_to_torus(radii, angles, minor_radius=0.3):
        R = 1.5
        x = (R + radii * np.cos(angles)) * np.cos(angles)
        y = (R + radii * np.cos(angles)) * np.sin(angles)
        z = radii * np.sin(angles)
        return x, y, z

    x3d, y3d, z3d = spiral_to_torus(bit_driven_radii, theta)

    # Plotly 3D scatter
    fig = go.Figure(data=go.Scatter3d(
        x=x3d,
        y=y3d,
        z=z3d,
        mode='markers+lines',
        marker=dict(
            size=4,
            color=diffs,
            colorscale='Plasma',
            colorbar=dict(title='Δ'),
            opacity=0.85
        ),
        line=dict(color='rgba(120,120,120,0.3)', width=2)
    ))

    fig.update_layout(
        title=f"3D Toroidal Harmonic Spiral (First {num_digits} Digits of Pi)",
        scene=dict(xaxis=dict(visible=False), yaxis=dict(visible=False), zaxis=dict(visible=False)),
        height=900,
        width=900,
        showlegend=False
    )

    fig.show()

# 🚀 Example usage
generate_pi_spiral(911512)  # Just change the number to any digit count you want!

