
import plotly.graph_objects as go
import numpy as np

def plot_triangle_animated(triangle_sets):
    """
    Animates triangle sets with tape-style controls: Play, Pause, Forward, Backward, Stop.
    """
    def triangle_coords(a, b, c):
        A = np.array([0, 0])
        B = np.array([c, 0])
        try:
            cos_angle = (a**2 + b**2 - c**2) / (2 * a * b)
            cos_angle = np.clip(cos_angle, -1.0, 1.0)
            angle_C = np.arccos(cos_angle)
            C_x = b * np.cos(angle_C)
            C_y = b * np.sin(angle_C)
            C = np.array([C_x, C_y])
        except:
            C = np.array([0, 0])
        return A, B, C

    frames = []
    for i, (a, b, c) in enumerate(triangle_sets):
        A, B, C = triangle_coords(a, b, c)
        x = [A[0], B[0], C[0], A[0]]
        y = [A[1], B[1], C[1], A[1]]
        frames.append(go.Frame(data=[go.Scatter(x=x, y=y, mode='lines+markers', name=f'Triangle {i+1}')],
                               name=str(i)))

    A, B, C = triangle_coords(*triangle_sets[0])
    x = [A[0], B[0], C[0], A[0]]
    y = [A[1], B[1], C[1], A[1]]

    fig = go.Figure(
        data=[go.Scatter(x=x, y=y, mode='lines+markers')],
        layout=go.Layout(
            title="Triangle Animator - Tape Control Style",
            xaxis=dict(range=[-10, 10], autorange=False),
            yaxis=dict(range=[-10, 10], autorange=False),
            updatemenus=[dict(
                type="buttons",
                buttons=[
                    dict(label="Play ▶️", method="animate", args=[None, {"frame": {"duration": 1000}, "fromcurrent": True}]),
                    dict(label="Pause ⏸️", method="animate", args=[[None], {"frame": {"duration": 0}, "mode": "immediate"}]),
                    dict(label="<< Back ⏮️", method="animate", args=[[str(i) for i in range(len(frames)-1, -1, -1)], {"frame": {"duration": 300}, "mode": "immediate"}]),
                    dict(label="Forward ⏭️", method="animate", args=[[str(i) for i in range(len(frames))], {"frame": {"duration": 300}, "mode": "immediate"}]),
                    dict(label="Stop ⏹️", method="animate", args=[[None], {"frame": {"duration": 0, "redraw": True}, "mode": "immediate"}])
                ],
                direction="right",
                pad={"r": 10, "t": 87},
                showactive=True,
                x=0.1,
                xanchor="left",
                y=1.2,
                yanchor="top"
            )]
        ),
        frames=frames
    )

    fig.show()

# Example usage:
triangle_sets = [(3, 1, 4), (1,4,2), (4, 2, 2), (2, 2, 5), (2, 1, 5), (4,1,5),(1,5,9),(5,9,2),(9,2,6), (2,6,5)]
plot_triangle_animated(triangle_sets)
