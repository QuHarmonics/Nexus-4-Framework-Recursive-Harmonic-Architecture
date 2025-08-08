import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import hilbert
from scipy.fft import fft, fftfreq
import pandas as pd

# Sample SHA-256 polarity wave from earlier conversion
# Replace this with your full 256-bit polarity wave if needed
polarity_wave = [-1 if bit == '0' else 1 for bit in bin(int(
    '2cf24dba5fb0a30e26e83b2ac5b9e29e1b161e5c1fa7425e73043362938b9824', 16))[2:].zfill(256)]

# Generate time-domain synthetic harmonic signal
t = np.linspace(0, 1, len(polarity_wave), endpoint=False)
signal = np.array(polarity_wave) * np.sin(2 * np.pi * 5 * t)

# Hilbert transform for amplitude and phase
analytic_signal = hilbert(signal)
amplitude_envelope = np.abs(analytic_signal)
instantaneous_phase = np.unwrap(np.angle(analytic_signal))

# FFT to extract frequencies
yf = fft(signal)
xf = fftfreq(len(t), 1 / len(t))
pos_mask = xf > 0
xf = xf[pos_mask]
yf = np.abs(yf[pos_mask])

# Build DataFrame of frequency + amplitude
freq_df = pd.DataFrame({'Freq (Hz)': xf, 'Amplitude': yf})
top_freqs = freq_df.sort_values('Amplitude', ascending=False).head(20)

# Calculate Pythagorean energy vector (hypotenuse)
top_freqs['Energy Vector'] = np.sqrt(top_freqs['Freq (Hz)']**2 + top_freqs['Amplitude']**2)

# Inferred frequency (reverse from energy vector + amplitude)
top_freqs['Inferred Frequency'] = np.sqrt(top_freqs['Energy Vector']**2 - top_freqs['Amplitude']**2)

# Difference vs original frequency
top_freqs['Frequency Error'] = top_freqs['Freq (Hz)'] - top_freqs['Inferred Frequency']

# Vector-Amplitude difference and Drift Ratio
top_freqs['Vector-Amplitude Difference'] = top_freqs['Energy Vector'] - top_freqs['Amplitude']
top_freqs['Drift Ratio'] = top_freqs['Vector-Amplitude Difference'].diff()

# Plot the results
plt.figure(figsize=(12, 6))
plt.plot(top_freqs['Freq (Hz)'], top_freqs['Vector-Amplitude Difference'], label='Vector-Amplitude Difference', marker='o')
plt.plot(top_freqs['Freq (Hz)'], top_freqs['Drift Ratio'], label='Drift Ratio', marker='x', linestyle='--')
plt.title('Harmonic Vector-Amplitude Dynamics from SHA Polarity Wave')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Value')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
