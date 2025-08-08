import numpy as np
import plotly.graph_objects as go
from itertools import product

# Input data
waveform_steps = [1, 4, 1, 5, 9, 2, 6, 5]  # Sequence of ASM-derived values
ratios = [3, 1, 4, 1, 5, 9, 2, 6]  # π-digit-derived ratios

# Time resolution for plotting
time = np.linspace(0, len(waveform_steps) - 1, 500)

# Normalize the inputs to create amplitude, frequency, and modulation
stack_amplitude = np.interp(waveform_steps, (min(waveform_steps), max(waveform_steps)), (0, 1))
stack_frequency = np.interp(waveform_steps, (min(waveform_steps), max(waveform_steps)), (0.1, 2))
pi_modulator = np.interp(ratios, (min(ratios), max(ratios)), (0.5, 1.5))

# Resample the inputs to match the fine time resolution
stack_amplitude_resampled = np.interp(time, np.linspace(0, len(waveform_steps) - 1, len(waveform_steps)), stack_amplitude)
stack_frequency_resampled = np.interp(time, np.linspace(0, len(waveform_steps) - 1, len(waveform_steps)), stack_frequency)
pi_modulator_resampled = np.interp(time, np.linspace(0, len(ratios) - 1, len(ratios)), pi_modulator)

# Prepare the resampled inputs
inputs_resampled = [stack_amplitude_resampled, stack_frequency_resampled, pi_modulator_resampled]
combinations = list(product(inputs_resampled, repeat=3))  # Generate combinations of amplitude, frequency, modulation

# Create the Plotly figure
fig = go.Figure()

# Iterate through all combinations and add them to the figure
for i, (amplitude, frequency, modulator) in enumerate(combinations[:27]):  # Limit to the first 27 combinations for clarity
    # Generate the waveform
    combined_waveform = amplitude * np.sin(2 * np.pi * frequency * time / len(time)) * modulator

    # Add the waveform as a line trace
    fig.add_trace(go.Scatter(
        x=time,
        y=combined_waveform,
        mode='lines',
        line=dict(width=1.5),
        name=f"Combo {i+1}"
    ))

# Highlight reference lines
fig.add_trace(go.Scatter(
    x=time,
    y=[0.5] * len(time),
    mode='lines',
    line=dict(dash='dash', color='green'),
    name="0.5 Reference Line"
))
fig.add_trace(go.Scatter(
    x=time,
    y=[0] * len(time),
    mode='lines',
    line=dict(dash='dash', color='black'),
    name="Zero Line"
))

# Update layout for better aesthetics
fig.update_layout(
    title="Interactive Visualization of ASM and π Modulation",
    xaxis_title="Time (Normalized)",
    yaxis_title="Waveform Amplitude",
    legend_title="Waveform Combinations",
    template="plotly_dark",  # Dark mode for better contrast
    font=dict(size=12),
    showlegend=True,
    width=1000,
    height=600
)

# Show the interactive plot
fig.show()
