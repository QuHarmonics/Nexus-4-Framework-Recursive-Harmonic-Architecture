import numpy as np

def qft(data, H, F, t):
    # Quantum Folding Tool (QFT)
    folded_data = np.sum(data / np.exp(H * F * t))
    return folded_data

def qut(folded_data, theta, zeta):
    # Quantum Unfolding Tool (QUT)
    unfolded_data = np.sum(folded_data * np.cos(theta)) + zeta
    return unfolded_data

# Example usage:
data = np.random.rand(100)  # Electromagnetic field data
H = 1.0  # Harmonic constant
F = 0.5  # Folding factor
t = 2  # Recursive depth
theta = np.pi / 4  # Phase correction
zeta = 0.1  # Residual error

folded_data = qft(data, H, F, t)
unfolded_data = qut(folded_data, theta, zeta)

# Calculate compression ratio and reconstruction error
compression_ratio = len(data) / len(str(folded_data))
reconstruction_error = np.mean(np.abs(unfolded_data - data))

print("Folded data:", folded_data)
print("Unfolded data:", unfolded_data)
print("Compression ratio:", compression_ratio)
print("Reconstruction error:", reconstruction_error)