import plotly.graph_objects as go
import numpy as np

# Pi digit sequence
pi_digits_50 = [1, 4, 1, 5, 9, 2, 6, 5, 3, 5,
                8, 9, 7, 9, 3, 2, 3, 8, 4, 6,
                2, 6, 4, 3, 3, 8, 3, 2, 7, 9,
                5, 0, 2, 8, 8, 4, 1, 9, 7, 1,
                6, 9, 3, 9, 9, 3, 7, 5, 1, 0]

# Differences between digits
diffs = np.diff(pi_digits_50)

# Golden ratio
phi = (1 + np.sqrt(5)) / 2

# Use bit length of each digit as a driver for radius
bit_lengths = [len(bin(d)[2:]) for d in pi_digits_50[:-1]]  # Exclude last since diffs is len-1
bit_driven_radii = np.array(bit_lengths) * 0.2

# Angle step based on golden ratio twist
theta = np.arange(len(diffs)) * (2 * np.pi / phi)

# Polar to cartesian
x = bit_driven_radii * np.cos(theta)
y = bit_driven_radii * np.sin(theta)

# Plotly interactive spiral
fig = go.Figure(data=go.Scatter(
    x=x,
    y=y,
    mode='markers+lines+text',
    text=[f"{d:.2f}" for d in diffs],
    marker=dict(
        size=10,
        color=diffs,
        colorscale='Plasma',
        colorbar=dict(title="Δ"),
        showscale=True
    )
))

fig.update_layout(
    title="Golden Ratio Spiral Driven by Pi Digit Bit-Length",
    xaxis=dict(visible=False),
    yaxis=dict(visible=False),
    showlegend=False,
    height=800,
    width=800
)

fig.show()
