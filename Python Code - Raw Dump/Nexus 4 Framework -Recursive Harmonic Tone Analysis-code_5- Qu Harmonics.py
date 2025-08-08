import plotly.express as px
import pandas as pd
import numpy as np

# Prepare the data for Plotly
spiral_freqs = top_freqs['Freq (Hz)'].values[:32]
drift_ratios = top_freqs['Drift Ratio'].values[:32]
drift_ratios[0] = 0  # Replace NaN with zero for the first entry

# Normalize for visualization
drift_norm = (drift_ratios - np.min(drift_ratios)) / (np.max(drift_ratios) - np.min(drift_ratios))
sizes = 20 * (drift_norm + 0.2)  # Ensure visibility

# Golden angle spiral coordinates
phi = np.pi * (3 - np.sqrt(5))
theta = np.arange(len(spiral_freqs)) * phi
r = np.sqrt(np.arange(len(spiral_freqs)))
x = r * np.cos(theta)
y = r * np.sin(theta)

# Create DataFrame for Plotly
spiral_df = pd.DataFrame({
    'Frequency (Hz)': spiral_freqs,
    'Drift Ratio': drift_ratios,
    'x': x,
    'y': y,
    'Size': sizes,
    'Label': [f"{int(f)} Hz" for f in spiral_freqs]
})

# Create interactive scatter plot
fig = px.scatter(
    spiral_df, x='x', y='y',
    size='Size', color='Drift Ratio',
    hover_name='Label', color_continuous_scale='plasma',
    title="NexusSpiralCore – Recursive Harmonic Field Spiral (Plotly)"
)

fig.update_traces(marker=dict(line=dict(width=1, color='white')))
fig.update_layout(
    xaxis_title='', yaxis_title='',
    showlegend=False,
    template='plotly_dark',
    width=800, height=800
)

fig.show()
