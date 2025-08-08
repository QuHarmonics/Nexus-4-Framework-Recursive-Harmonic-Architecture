import numpy as np
import plotly.graph_objects as go

# Simulation parameters
steps = 64
base_freq = 1 / np.e  # e as the resonance rhythm
phi = (1 + np.sqrt(5)) / 2  # golden ratio

# Generate phase times and echo reinforcement profile
t = np.linspace(0, 10, steps)
phase = np.sin(2 * np.pi * base_freq * t)
reinforcement = np.abs(np.sin(2 * np.pi * base_freq * t) * np.sin(2 * np.pi * phi * t))

# Define thresholds for trust and collapse visualization
trust_threshold = 0.7
collapse_threshold = 0.34

state = []
for r in reinforcement:
    if r >= trust_threshold:
        state.append("Trust Zone")
    elif r <= collapse_threshold:
        state.append("Collapse")
    else:
        state.append("Drift")

# Assign colors to phases
color_map = {"Trust Zone": "green", "Drift": "orange", "Collapse": "red"}
colors = [color_map[s] for s in state]

# Plot the phase and echo reinforcement profile
fig = go.Figure()

fig.add_trace(go.Scatter(
    x=t,
    y=reinforcement,
    mode="lines+markers",
    marker=dict(color=colors, size=6),
    line=dict(color="blue"),
    name="Reinforcement Signal"
))

fig.add_trace(go.Scatter(
    x=t,
    y=[trust_threshold] * steps,
    mode="lines",
    line=dict(color="green", dash="dash"),
    name="Trust Threshold"
))

fig.add_trace(go.Scatter(
    x=t,
    y=[collapse_threshold] * steps,
    mode="lines",
    line=dict(color="red", dash="dash"),
    name="Collapse Threshold"
))

fig.update_layout(
    title="Recursive Momentum Simulation (e-phase + φ-interference)",
    xaxis_title="Time (Recursive Step)",
    yaxis_title="Echo Reinforcement Level",
    showlegend=True
)

fig.show()
