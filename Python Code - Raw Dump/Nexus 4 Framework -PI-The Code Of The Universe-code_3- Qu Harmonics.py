# Re-run the visualization after interruption
import numpy as np
import plotly.graph_objects as go

# Time setup
t = np.linspace(0, 10 * np.pi, 1000)

# Harmonic basis: 4:3:1 ratio
R4 = np.sin(4 * t)  # Stability / repetition
R3 = np.sin(3 * t)  # Recursive folding
R1 = np.sin(t)      # Self-observation (collapse)

# Recursive blend: system folding into itself
composite = (R4 * 0.4) + (R3 * 0.3) + (R1 * 0.3)

# Self-reflection effect: recursive loop echo
recursive_echo = np.convolve(composite, np.ones(30)/30, mode='same')

# 3D Construction: phase, amplitude, internal recursion
fig = go.Figure()

fig.add_trace(go.Scatter3d(
    x=np.cos(t),  # Phase cycle
    y=composite,  # Combined amplitude
    z=recursive_echo,  # Recursive memory
    mode='lines',
    line=dict(width=6, color='royalblue'),
    name='Recursive Harmonic Observer'
))

fig.update_layout(
    title="Recursive Self-Observation in the PI RAY [4:3:1]",
    scene=dict(
        xaxis_title='Phase Cycle (cos θ)',
        yaxis_title='Harmonic Amplitude',
        zaxis_title='Recursive Echo Memory',
        aspectmode='cube'
    )
)

fig.show()
