import numpy as np
import plotly.graph_objects as go
from itertools import product

# Data setup
waveform_steps = [1, 4, 1, 5, 9, 2, 6, 5]
ratios = [3, 1, 4, 1, 5, 9, 2, 6]
time = np.linspace(0, len(waveform_steps) - 1, 5000)

# Normalize input series
stack_amplitude = np.interp(waveform_steps, (min(waveform_steps), max(waveform_steps)), (0, 1))
stack_frequency = np.interp(waveform_steps, (min(waveform_steps), max(waveform_steps)), (0.1, 2))
pi_modulator = np.interp(ratios, (min(ratios), max(ratios)), (0.5, 2.5))

# Resample to match fine time
resample = lambda data: np.interp(time, np.linspace(0, len(data) - 1, len(data)), data)
A = resample(stack_amplitude)
F = resample(stack_frequency)
M = resample(pi_modulator)

# All triplet combinations
combos = list(product([A, F, M], repeat=3))

# Build frames for animation
frames = []
for i, (a, f, m) in enumerate(combos):
    y = a * np.sin(2 * np.pi * f * time / len(time)) * m
    frames.append(go.Frame(data=[go.Scatter(x=time, y=y, mode='lines')],
                           name=f"Combo {i+1}"))

# Initial trace (first combo)
initial_y = combos[0][0] * np.sin(2 * np.pi * combos[0][1] * time / len(time)) * combos[0][2]
trace = go.Scatter(x=time, y=initial_y, mode='lines', name="Waveform")

# Layout with animation controls
layout = go.Layout(
    title="Animated Harmonic Modulation - ASM × π Interactions",
    xaxis_title="Time (Normalized)",
    yaxis_title="Waveform Amplitude",
    updatemenus=[{
        "type": "buttons",
        "showactive": False,
        "buttons": [{
            "label": "Play",
            "method": "animate",
            "args": [None, {"frame": {"duration": 300, "redraw": True},
                            "fromcurrent": True}]
        }, {
            "label": "Pause",
            "method": "animate",
            "args": [[None], {"frame": {"duration": 0}, "mode": "immediate"}]
        }]
    }],
    sliders=[{
        "steps": [{
            "method": "animate",
            "label": f"Combo {i+1}",
            "args": [[f"Combo {i+1}"], {"frame": {"duration": 0, "redraw": True}}]
        } for i in range(len(combos))],
        "transition": {"duration": 0},
        "x": 0, "y": -0.15,
        "currentvalue": {"prefix": "Combo: "}
    }]
)

# Final figure
fig = go.Figure(data=[trace], layout=layout, frames=frames)
fig.show()
