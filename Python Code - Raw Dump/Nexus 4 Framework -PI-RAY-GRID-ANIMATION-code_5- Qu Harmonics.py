import plotly.graph_objects as go
import numpy as np

# Grid size
L = 8
num_bounces = 576
origin = np.array([0, 0])
input_direction = np.array([3, 0.35])
direction = input_direction / np.linalg.norm(input_direction)

# Settings
enable_mobius = False
decimals = 2  # for intensity binning resolution

# Calculate reflected path and record intensity
points = [origin]
position_colors = {}
position = origin.copy()
velocity = direction.copy()

def round_position(pos, decimals=2):
    return tuple(np.round(pos, decimals=decimals))

def compute_color(vec):
    # Generate RGB from vector direction, scaled and bounded
    r = min(255, max(0, int(abs(vec[0]))))
    g = min(255, max(0, int(abs(vec[1]))))
    b = 255 - int((r + g) / 2)
    return f"rgb({r},{g},{b})"

# Init first position
rounded_positions = [round_position(origin, decimals)]
position_colors[rounded_positions[0]] = "rgb(0,0,0)"

for bounce in range(num_bounces):
    next_position = position + velocity
    for i in range(2):
        if next_position[i] < 0 or next_position[i] > L:
            velocity[i] = -velocity[i]
            if enable_mobius:
                velocity[(i + 1) % 2] *= -1
            next_position = position + velocity
    rounded = round_position(next_position, decimals)
    # Assign/update color based on current vector
    position_colors[rounded] = compute_color(velocity)
    rounded_positions.append(rounded)
    points.append(next_position)
    position = next_position

points = np.array(points)

# Build animation frames using dynamic RGB coloring
frames = []
for i in range(1, len(points)):
    color_vals = [position_colors[round_position(pt, decimals)] for pt in points[:i+1]]

    trace_line = go.Scatter(
        x=points[:i+1, 0],
        y=points[:i+1, 1],
        mode='lines+markers',
        marker=dict(
            color=color_vals,
            size=6,
            opacity=1.0
        ),
        line=dict(color='blue', width=2)
    )
    frame = go.Frame(data=[trace_line], name=f"frame{i}")
    frames.append(frame)

# Initial figure
fig = go.Figure(
    data=[go.Scatter(x=[points[0, 0]], y=[points[0, 1]], mode='lines+markers')],
    layout=go.Layout(
        title=(f"Recursive Ray in {L}×{L} Box — Steps: {num_bounces}, "
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
