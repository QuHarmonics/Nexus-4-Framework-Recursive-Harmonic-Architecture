import numpy as np
import matplotlib.pyplot as plt

# Define parameters
N = 1000  # Number of self-similar patterns
S = 2     # Scaling factor
omega = 2048 # Frequency
t = np.linspace(0, 10, 1000)  # Time array

# Calculate fractal dimension
D = np.log(N) / np.log(S)

# Generate harmonic oscillation
x = np.sin(omega * t)

# Generate wavelet transform
wavelet_coeffs = np.zeros((N, len(t)))
for i in range(N):
    wavelet_coeffs[i] = np.sin((i+1) * omega * t)

# Plot results
plt.figure(figsize=(12, 6))

plt.subplot(1, 3, 1)
plt.plot(t, x, label='Harmonic Oscillation')
plt.title('Harmonic Oscillation')
plt.xlabel('Time')
plt.ylabel('Amplitude')
plt.legend()

plt.subplot(1, 3, 2)
plt.imshow(wavelet_coeffs, cmap='inferno', aspect='auto')
plt.title('Wavelet Transform')
plt.xlabel('Time')
plt.ylabel('Scale')

plt.subplot(1, 3, 3)
plt.plot(t, np.abs(np.sum(wavelet_coeffs, axis=0)), label='Wavelet Coefficients')
plt.title('Wavelet Coefficients')
plt.xlabel('Time')
plt.ylabel('Amplitude')
plt.legend()

plt.tight_layout()
plt.show()