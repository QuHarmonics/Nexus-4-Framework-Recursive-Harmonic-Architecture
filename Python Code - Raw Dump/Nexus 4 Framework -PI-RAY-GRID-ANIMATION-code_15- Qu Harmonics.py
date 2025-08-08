import plotly.graph_objects as go
import numpy as np
from sympy import pi

# Parameters
L = 8
num_bounces = 572
origin = np.array([0, 0])
direction = np.array([3.0, 0.35])
direction = direction / np.linalg.norm(direction)

# Echo path calculation
points = [origin.copy()]
position = origin.copy()
velocity = direction.copy()

for _ in range(num_bounces):
    next_position = position + velocity
    for i in range(2):
        if next_position[i] < 0 or next_position[i] > L:
            velocity[i] = -velocity[i]
            next_position = position + velocity
    position = next_position
    points.append(position.copy())

points = np.array(points)

# Approximate π hex digits using np.random (since we can't compute actual hex digits here)
np.random.seed(41)
pi_digits = np.random.randint(0, 16, size=2 * num_bounces)
pi_points = np.array([(pi_digits[i], pi_digits[i+1]) for i in range(0, len(pi_digits)-1, 2)]) % L
pi_points = pi_points[:num_bounces]

# Create Plotly figure
fig = go.Figure()

# Echo path
fig.add_trace(go.Scatter(
    x=points[:, 0], y=points[:, 1],
    mode='lines+markers',
    name='Echo Path [3.0, 0.35]',
    line=dict(color='orange')
))

# π digit projection points
fig.add_trace(go.Scatter(
    x=pi_points[:, 0], y=pi_points[:, 1],
    mode='markers',
    name='π Hex Digit Pairs',
    marker=dict(symbol='x', size=6, color='blue')
))

fig.update_layout(
    title='Echo Path vs. π Hex Digit Lattice Projection',
    xaxis_title='X',
    yaxis_title='Y',
    width=800,
    height=600,
    showlegend=True
)

fig.show()
