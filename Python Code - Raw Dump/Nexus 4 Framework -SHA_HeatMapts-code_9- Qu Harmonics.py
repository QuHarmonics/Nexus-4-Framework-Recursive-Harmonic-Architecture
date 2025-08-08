import plotly.graph_objects as go
import numpy as np

# Input 8x8 array (not used in the plot, but defined for completeness)
a = np.array([
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0]
])

# 8x8 grid size
N = 8

# Initial position and normalized direction
origin = np.array([0, 0])
direction = np.array([1, 1]) / np.sqrt(2)

# Initialize figure
fig = go.Figure()

# Plot the 8x8 grid
for i in range(N + 1):
    fig.add_trace(go.Scatter(x=[i, i], y=[0, N], mode='lines', line=dict(color='black', width=1)))
    fig.add_trace(go.Scatter(x=[0, N], y=[i, i], mode='lines', line=dict(color='black', width=1)))

# Plot the ray
t_max = 2 * N  # Maximum distance to trace the ray
t = 0
current_pos = origin.copy()

while t < t_max:
    fig.add_trace(go.Scatter(x=[current_pos[0]], y=[current_pos[1]], mode='markers', marker=dict(color='red', size=5)))
    t += 0.1
    current_pos = origin + t * direction

    # Check if the ray hits a grid boundary
    if current_pos[0] < 0 or current_pos[0] > N or current_pos[1] < 0 or current_pos[1] > N:
        break

# Update layout
fig.update_layout(
    title=f"Recursive Ray 8x8 Boxes: {N*N} Dir: ({direction[0]:.4f}, {direction[1]:.4f}) Norm Dir: ({np.linalg.norm(direction):.4f}, {np.linalg.norm(direction):.4f})",
    xaxis=dict(title="X", range=[-1, N+1], tickmode='linear', tick0=0, dtick=1),
    yaxis=dict(title="Y", range=[-1, N+1], tickmode='linear', tick0=0, dtick=1),
    showlegend=False
)

# Show the plot
fig.show()