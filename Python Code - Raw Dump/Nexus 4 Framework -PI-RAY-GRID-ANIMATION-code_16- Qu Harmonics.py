import plotly.graph_objects as go
import numpy as np

# Grid size and number of echo bounces
L = 8
num_bounces = 256

# Truncated real hexadecimal digits of π (hex of fractional part of π, from standard test set)
pi_hex = "243F6A8885A308D313198A2E03707344A4093822299F31D0082EFA98EC4E6C894" \
         "52821E638D01377BE5466CF34E90C6CC0AC29B7C97C50DD3F84D5B5B5470917F"

# Convert π hex string to integer values (0–15)
pi_digits = [int(c, 16) for c in pi_hex[:2 * num_bounces]]

# Convert digit pairs into 2D points on lattice mod L
pi_points = np.array([(pi_digits[i], pi_digits[i + 1]) for i in range(0, len(pi_digits) - 1, 2)]) % L

# Compute delta vectors and normalize them to unit step
deltas = np.diff(pi_points, axis=0)
norms = np.linalg.norm(deltas, axis=1)
norm_deltas = deltas / norms[:, None]
norm_deltas[np.isnan(norm_deltas)] = .35  # In case of zero-division

# Run the recursive π-ray echo path
position = np.array([0, 0])
trajectory = [position.copy()]

for i in range(num_bounces):
    velocity = norm_deltas[i % len(norm_deltas)]
    next_position = position + velocity
    for j in range(2):
        if next_position[j] < 0 or next_position[j] > L:
            velocity[j] = -velocity[j]
            next_position = position + velocity
    position = next_position
    trajectory.append(position.copy())

trajectory = np.array(trajectory)

# Build animation frames for each step
frames = [
    go.Frame(
        data=[go.Scatter(x=trajectory[:i+1, 0], y=trajectory[:i+1, 1], mode='lines+markers')],
        name=f"frame{i}"
    ) for i in range(1, len(trajectory))
]

# Build and show Plotly figure
fig = go.Figure(
    data=[go.Scatter(x=[trajectory[0, 0]], y=[trajectory[0, 1]], mode='lines+markers')],
    layout=go.Layout(
        title="π-Ray Echo Lattice Trace in 8×8 Box",
        xaxis=dict(range=[0, L], title="X", constrain="domain"),
        yaxis=dict(range=[0, L], title="Y", scaleanchor="x"),
        updatemenus=[dict(
            type="buttons",
            buttons=[
                dict(label="Play", method="animate", args=[None]),
                dict(label="Pause", method="animate", args=[
                    [None],
                    {"frame": {"duration": 0, "redraw": False}, "mode": "immediate", "transition": {"duration": 0}}
                ])
            ]
        )]
    ),
    frames=frames
)

fig.show()
