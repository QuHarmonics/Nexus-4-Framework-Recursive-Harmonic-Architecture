import numpy as np
from scipy.io.wavfile import write
import os

# Audio parameters
sample_rate = 44100  # Standard CD-quality sample rate
duration = 4  # seconds
n_samples = sample_rate * duration

# Time array
t = np.linspace(0, duration, n_samples, endpoint=False)

# Collapse spiral mapped to descending frequency (220Hz down)
H = 0.35
F = 1.0
phi = (1 + np.sqrt(5)) / 2
t_collapse = np.arange(n_samples) * phi / sample_rate
R_vals = np.exp(-H * F * t_collapse)

# Map R_vals to frequency in a symbolic range (e.g., 220Hz to ~20Hz)
freqs = 220 * R_vals  # Collapse from 220Hz downward

# Generate the tone
tone = np.sin(2 * np.pi * freqs * t)

# Normalize to 16-bit PCM
tone_normalized = np.int16(tone / np.max(np.abs(tone)) * 32767)

# Save the file
output_path = "/mnt/data/recursive_collapse_spiral.wav"
write(output_path, sample_rate, tone_normalized)

output_path
