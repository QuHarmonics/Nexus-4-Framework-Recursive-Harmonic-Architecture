import plotly.graph_objects as go
import numpy as np

# Grid size
L = 8
num_bounces = 572
origin = np.array([0, 0])
input_direction = np.array([3, 0.35])
direction = input_direction / np.linalg.norm(input_direction)

# Settings
enable_mobius = False
echo_tolerance = 1e-2

# Calculate reflected path
points = [origin]
echo_indices = []
position = origin.copy()
velocity = direction.copy()

for bounce in range(num_bounces):
    next_position = position + velocity
    for i in range(2):
        if next_position[i] < 0 or next_position[i] > L:
            velocity[i] = -velocity[i]
            if enable_mobius:
                velocity[(i + 1) % 2] *= -1
            next_position = position + velocity
    # Echo detection
    for j, past in enumerate(points):
        if np.linalg.norm(past - next_position) < echo_tolerance:
            echo_indices.append(len(points))
            break
    points.append(next_position)
    position = next_position

points = np.array(points)

# Build animation frames
frames = []
for i in range(1, len(points)):
    trail_opacity = np.linspace(0.2, 1.0, i+1)
    trace_line = go.Scatter(
        x=points[:i+1, 0],
        y=points[:i+1, 1],
        mode='lines+markers',
        marker=dict(
            color=['red' if j in echo_indices else 'blue' for j in range(i+1)],
            size=6,
            opacity=trail_opacity
        ),
        line=dict(color='blue', width=2)
    )
    frame = go.Frame(data=[trace_line], name=f"frame{i}")
    frames.append(frame)

# Initial figure
fig = go.Figure(
    data=[go.Scatter(x=[points[0, 0]], y=[points[0, 1]], mode='lines+markers')],
    layout=go.Layout(
        title=f"Recursive Ray With Echo {L}×{L} Box — Steps: {num_bounces}, "
              f"Input Dir: {input_direction[0]:.4f}, {input_direction[1]:.4f} → "
              f"Norm Dir: {direction[0]:.4f}, {direction[1]:.4f}",
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
