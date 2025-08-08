# Extended pi digits (256) for input
from mpmath import mp

# Set precision and get digits
mp.dps = 258  # extra precision for slicing
pi_str = str(mp.pi)[2:2050]  # skip "3."

# Convert to list of integers
pi_digits_256 = [int(d) for d in pi_str]

# Calculate differences between consecutive digits
diffs_256 = np.diff(pi_digits_256)

# Golden ratio
phi = (1 + np.sqrt(5)) / 2

# Use bit length of each digit as a driver for radius
bit_lengths_256 = [len(bin(d)[2:]) for d in pi_digits_256[:-1]]  # same length as diffs
bit_driven_radii_256 = np.array(bit_lengths_256) * 0.2

# Angle step via golden ratio
theta_256 = np.arange(len(diffs_256)) * (2 * np.pi / phi)

# Polar to Cartesian coordinates
x_256 = bit_driven_radii_256 * np.cos(theta_256)
y_256 = bit_driven_radii_256 * np.sin(theta_256)

# Interactive Plotly spiral
fig = go.Figure(data=go.Scatter(
    x=x_256,
    y=y_256,
    mode='markers+lines+text',
    text=[f"{d:.2f}" for d in diffs_256],
    marker=dict(
        size=8,
        color=diffs_256,
        colorscale='Plasma',
        colorbar=dict(title="Δ"),
        showscale=True
    )
))

fig.update_layout(
    title="Golden Ratio Spiral Driven by Bit-Length (First 256 Digits of Pi)",
    xaxis=dict(visible=False),
    yaxis=dict(visible=False),
    showlegend=False,
    height=900,
    width=900
)

fig.show()
