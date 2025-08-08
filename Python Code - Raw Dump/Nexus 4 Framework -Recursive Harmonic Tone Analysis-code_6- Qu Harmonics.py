import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import plotly.express as px
from scipy.signal import hilbert
from scipy.fft import fft, fftfreq

# SHA hash polarity extraction
sha_hex = '110b54243ea7d2c78b0a17a4adec85108379356e506259b7c14db36bc7eb1d32'
polarity_wave = [-1 if bit == '0' else 1 for bit in bin(int(sha_hex, 16))[2:].zfill(256)]

# Create base signal from polarity
t = np.linspace(0, 1, len(polarity_wave), endpoint=False)
signal = np.array(polarity_wave) * np.sin(2 * np.pi * 5 * t)

# Hilbert transform for envelope and phase
analytic_signal = hilbert(signal)
amplitude_envelope = np.abs(analytic_signal)

# Frequency analysis
yf = fft(signal)
xf = fftfreq(len(t), 1 / len(t))
xf = xf[xf > 0]
yf = np.abs(yf[:len(xf)])

# Create and trim DataFrame
freq_df = pd.DataFrame({'Freq (Hz)': xf, 'Amplitude': yf})
top_freqs = freq_df.sort_values('Amplitude', ascending=False).head(32).copy()

# Energy and drift calculations
top_freqs['Energy Vector'] = np.sqrt(top_freqs['Freq (Hz)']**2 + top_freqs['Amplitude']**2)
top_freqs['Vector-Amplitude Difference'] = top_freqs['Energy Vector'] - top_freqs['Amplitude']
top_freqs['Drift Ratio'] = top_freqs['Vector-Amplitude Difference'].diff().fillna(0)

# Spiral coordinates
spiral_freqs = top_freqs['Freq (Hz)'].values
drift_ratios = top_freqs['Drift Ratio'].values
drift_norm = (drift_ratios - np.min(drift_ratios)) / (np.max(drift_ratios) - np.min(drift_ratios))
sizes = 20 * (drift_norm + 0.2)

phi = np.pi * (3 - np.sqrt(5))
theta = np.arange(len(spiral_freqs)) * phi
r = np.sqrt(np.arange(len(spiral_freqs)))
x = r * np.cos(theta)
y = r * np.sin(theta)

# Build DataFrame
spiral_df = pd.DataFrame({
    'Frequency (Hz)': spiral_freqs,
    'Drift Ratio': drift_ratios,
    'x': x,
    'y': y,
    'Size': sizes,
    'Label': [f"{int(f)} Hz" for f in spiral_freqs]
})

# Plotly spiral plot
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
