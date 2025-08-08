# Re-import necessary modules after code execution state reset
import plotly.graph_objects as go
import numpy as np
import plotly.io as pio
import os

# Set the output path on the D: drive
output_path = "D:\\recursive_pi_triangle.html"

# Triangle vertices
triangle = np.array([
    [0, 0],
    [1, 0],
    [0.5, np.sqrt(3)/2],
    [0, 0]  # closing the triangle
])

# Bouncing inside the triangle
def reflect(p, v, a, b):
    """Reflect vector v at point p off the edge from a to b."""
    edge = b - a
    edge_norm = np.array([-edge[1], edge[0]])
    edge_norm = edge_norm / np.linalg.norm(edge_norm)
    v_reflected = v - 2 * np.dot(v, edge_norm) * edge_norm
    return v_reflected

# Initial point and velocity
point = np.array([0.1, 0.1])
velocity = np.array([0.03, 0.05])
points = [point.copy()]
num_bounces = 500

for _ in range(num_bounces):
    next_point = point + velocity

    for i in range(3):
        a = triangle[i]
        b = triangle[i+1]
        edge = b - a
        edge_norm = np.array([-edge[1], edge[0]])
        edge_norm = edge_norm / np.linalg.norm(edge_norm)

        to_line = np.dot((next_point - a), edge_norm)
        from_line = np.dot((point - a), edge_norm)

        if np.sign(to_line) != np.sign(from_line):
            velocity = reflect(point, velocity, a, b)
            next_point = point + velocity
            break

    point = next_point
    points.append(point.copy())

points = np.array(points)

# Create animation with plotly
fig = go.Figure()

fig.add_trace(go.Scatter(x=triangle[:, 0], y=triangle[:, 1],
                         mode="lines", name="Triangle"))

fig.add_trace(go.Scatter(x=points[:, 0], y=points[:, 1],
                         mode="lines", name="Path",
                         line=dict(color="cyan")))

fig.update_layout(
    title="Recursive π Triangle Reflection",
    xaxis=dict(visible=False),
    yaxis=dict(visible=False),
    showlegend=False,
    plot_bgcolor='black',
    paper_bgcolor='black',
    width=600,
    height=600
)

# Save as HTML
pio.write_html(fig, file=output_path, auto_open=False)
output_path
