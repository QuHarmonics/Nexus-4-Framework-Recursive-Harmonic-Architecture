import plotly.graph_objects as go
import numpy as np

# Grid size
L = 8
num_bounces = 572

# Initial position and normalized direction
origin = np.array([0, 0])
input_direction = np.array([2, .35])
direction = input_direction / np.linalg.norm(input_direction)
direction = direction / np.linalg.norm(direction)

# Calculate reflected path
points = [origin]
position = origin.copy()
velocity = direction.copy()

for _ in range(num_bounces):
    next_position = position + velocity
    for i in range(2):
        if next_position[i] < 0 or next_position[i] > L:
            velocity[i] = -velocity[i]
            next_position = position + velocity
    points.append(next_position)
    position = next_position

points = np.array(points)

# Initial figure
fig = go.Figure(
    data=[go.Scatter(x=[points[0, 0]], y=[points[0, 1]], mode='lines+markers')],
    layout=go.Layout(
        title=(f"Recursive{L}×{L} Box — Steps: {num_bounces}, "
           f"Input Dir: {input_direction[0]:.4f}, {input_direction[1]:.4f} → "
           f"Norm Dir: {direction[0]:.4f}, {direction[1]:.4f}"),
        xaxis=dict(range=[0, L], title="X", constrain="domain"),
        yaxis=dict(range=[0, L], title="Y", scaleanchor="x"),
        updatemenus=[dict(
            type="buttons",
            buttons=[
                dict(label="Play", method="animate", args=[None]),
                dict(label="Pause", method="animate", args=[
                    [None],
                    {"frame": {"duration": 0, "redraw": False},
                     "mode": "immediate",
                     "transition": {"duration": 0}}
                ])
            ]
        )]
    ),
    frames=frames
)

fig.show()
