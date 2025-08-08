# Re-importing necessary libraries and re-defining the constants due to execution state reset.
import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, fftfreq

# Constants
K = [
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2
]

# Normalize constants to [0, 1]
K_normalized = [k / (2**32) for k in K]

# Generate x values (angle)
x = np.linspace(0, 2 * np.pi, 1000)

# Individual waveforms
individual_waveforms = [np.sin(x + kn) for kn in K_normalized]

# Composite waveform
composite_waveform = np.sum(individual_waveforms, axis=0)

# Fourier Transform of the composite waveform
fft_result = fft(composite_waveform)
frequencies = fftfreq(len(x), (2 * np.pi) / len(x))

# Visualization
# 1. Individual Waveforms
plt.figure(figsize=(12, 8))
for i, waveform in enumerate(individual_waveforms[:8]):  # Show first 8 waveforms for clarity
    plt.plot(x, waveform, label=f"Waveform {i+1}")
plt.title("Individual Waveforms Derived from SHA-256 Constants")
plt.xlabel("Angle (radians)")
plt.ylabel("Amplitude")
plt.legend()
plt.grid()
plt.show()

# 2. Composite Waveform
plt.figure(figsize=(12, 6))
plt.plot(x, composite_waveform, label="Composite Waveform", color="blue")
plt.title("Composite Waveform of SHA-256 Constants")
plt.xlabel("Angle (radians)")
plt.ylabel("Amplitude")
plt.legend()
plt.grid()
plt.show()

# 3. Frequency Domain (Fourier Transform)
plt.figure(figsize=(12, 6))
plt.plot(frequencies[:len(frequencies)//2], np.abs(fft_result)[:len(frequencies)//2], color="orange")
plt.title("Frequency Spectrum of Composite Waveform")
plt.xlabel("Frequency")
plt.ylabel("Amplitude (Magnitude)")
plt.grid()
plt.show()
