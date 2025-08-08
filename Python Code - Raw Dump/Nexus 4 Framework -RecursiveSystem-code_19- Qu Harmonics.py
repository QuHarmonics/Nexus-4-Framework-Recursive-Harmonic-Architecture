# Updated version: streamed construction using Python generators to handle large π datasets

import numpy as np
from mpmath import mp
import plotly.graph_objects as go

def stream_pi_digits(n):
    """Generate n digits of pi as a stream of integers."""
    mp.dps = n + 2
    pi_str = str(mp.pi)[2:n+2]
    for ch in pi_str:
        yield int(ch)

def generate_pi_spiral_streamed(num_digits=100000, chunk_size=10000):
    """Stream-compatible spiral generator for very large digit counts."""
    pi_stream = stream_pi_digits(num_digits)

    diffs = []
    radii = []
    angles = []

    prev_digit = next(pi_stream)
    index = 0

    phi = (1 + np.sqrt(5)) / 2

    for current_digit in pi_stream:
        diff = current_digit - prev_digit
        diffs.append(diff)
        bit_length = len(bin(current_digit)[2:])
        radii.append(bit_length * 0.2)
        angles.append(index * (2 * np.pi / phi))
        prev_digit = current_digit
        index += 1

    radii = np.array(radii)
    angles = np.array(angles)
    diffs = np.array(diffs)

    def spiral_to_torus(radii, angles, minor_radius=0.3):
        R = 1.5
        x = (R + radii * np.cos(angles)) * np.cos(angles)
        y = (R + radii * np.cos(angles)) * np.sin(angles)
        z = radii * np.sin(angles)
        return x, y, z

    x3d, y3d, z3d = spiral_to_torus(radii, angles)

    fig = go.Figure(data=go.Scatter3d(
        x=x3d,
        y=y3d,
        z=z3d,
        mode='markers',
        marker=dict(
            size=2,
            color=diffs,
            colorscale='Plasma',
            colorbar=dict(title='Δ'),
            opacity=0.8
        )
    ))

    fig.update_layout(
        title=f"3D Toroidal Harmonic Spiral (First {num_digits} Digits of Pi)",
        scene=dict(xaxis=dict(visible=False), yaxis=dict(visible=False), zaxis=dict(visible=False)),
        height=900,
        width=900,
        showlegend=False
    )

    fig.show()

# Example usage for streamed pi visualization
generate_pi_spiral_streamed(100000)

