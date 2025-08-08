import math
import plotly.graph_objects as go

# First 8 digits of π as provided: 14159625
pi_byte = [1, 4, 1, 5, 9, 6, 2, 5]

# Convert these digits into a normalized "waveform" for π
pi_wave = [value / 9.0 for value in pi_byte]  # Normalize values (max is 9 for single-digit integers)

# Generate a "math waveform" for comparison: simple operation (e.g., addition of pairs)
math_wave = [(pi_byte[i] + pi_byte[i + 1]) / 18.0 for i in range(len(pi_byte) - 1)]  # Normalize again

# Prepare the 3D wave visualization
fig = go.Figure()

# π waveform
fig.add_trace(go.Scatter(
    x=list(range(len(pi_wave))),
    y=pi_wave,
    mode='lines+markers',
    name="π Waveform (Normalized)",
    line=dict(color='blue', width=2)
))

# Math waveform
fig.add_trace(go.Scatter(
    x=list(range(len(math_wave))),
    y=math_wave,
    mode='lines+markers',
    name="Math Waveform",
    line=dict(color='red', width=2)
))

# Update layout for visualization
fig.update_layout(
    title="Comparison of π Waveform and Math Waveform",
    xaxis_title="Index",
    yaxis_title="Normalized Value",
    height=600,
    width=800
)

# Show the plot
fig.show()
