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

# Prepare the node positions and their respective colors
node_x = []
node_y = []
node_colors = []

for pt in rounded_positions:
    node_x.append(pt[0])
    node_y.append(pt[1])
    node_colors.append(position_colors[pt])

# Create the heatmap
fig = go.Figure(
    data=[
        go.Scatter(
            x=node_x,
            y=node_y,
            mode='markers',
            marker=dict(
                size=10,
                color=node_colors,  # Assign the color based on position colors
                opacity=1.0,
                colorscale='Viridis',  # You can also change to another scale like 'Cividis'
                colorbar=dict(title='Intensity')
            )
        )
    ],
    layout=go.Layout(
        title=f"Recursive Ray in {L}×{L} Box — Steps: {num_bounces}, "
              f"Input Dir: {input_direction[0]:.4f}, {input_direction[1]:.4f} → "
              f"Norm Dir: {direction[0]:.4f}, {direction[1]:.4f}",
        xaxis=dict(range=[0, L], title="X", constrain="domain"),
        yaxis=dict(range=[0, L], title="Y", scaleanchor="x"),
        showlegend=False
    )
)

fig.show()
