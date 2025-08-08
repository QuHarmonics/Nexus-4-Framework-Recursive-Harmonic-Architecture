import numpy as np
from mpmath import mp
import plotly.graph_objects as go

# 🔧 Set this to control how many digits of Pi you want to use
pi_length = 256  # Change this to 512, 1024, etc.

# Set precision with a bit of buffer
mp.dps = pi_length + 2
pi_str = str(mp.pi)[2:2 + pi_length]  # skip "3."

# Convert to list of integers
pi_digits = [int(d) for d in pi_str]

# Calculate differences between consecutive digits
diffs = np.diff(pi_digits)

# Golden ratio for angular stepping
phi = (1 + np.sqrt(5)) / 2

# Bit lengths (binary) as radial weights
bit_lengths = [len(bin(d)[2:]) for d in pi_digits[:-1]]  # match diff length
bit_driven_radii = np.array(bit_lengths) * 0.2

# Angular step via golden ratio spiral
theta = np.arange(len(diffs)) * (2 * np.pi / phi)

# Polar to Cartesian
x = bit_driven_radii * np.cos(theta)
y = bit_driven_radii * np.sin(theta)

# 🌀 Plot with Plotly
fig = go.Figure(data=go.Scatter(
    x=x,
    y=y,
    mode='markers+lines+text',
    text=[f"Δ={d}" for d in diffs],
    marker=dict(
        size=8,
        color=diffs,
        colorscale='Plasma',
        colorbar=dict(title="Δ"),
        showscale=True
    )
))

fig.update_layout(
    title=f"Golden Ratio Spiral Driven by Bit-Length (First {pi_length} Digits of Pi)",
    xaxis=dict(visible=False),
    yaxis=dict(visible=False),
    showlegend=False,
    height=900,
    width=900
)

fig.show()
