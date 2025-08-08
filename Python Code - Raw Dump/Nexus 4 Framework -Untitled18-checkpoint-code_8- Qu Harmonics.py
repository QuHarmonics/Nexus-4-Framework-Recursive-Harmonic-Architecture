import numpy as np

def qft(data, H, F, t):
    # Quantum Folding Tool (QFT)
    folded_data = np.sum(data / np.exp(H * F * t))
    return folded_data

# Example usage:
data = np.random.rand(100)  # Electromagnetic field data
H = 1.0  # Harmonic constant
F = 0.65  # Folding factor
t = 3.0  # Recursive depth

folded_data = qft(data, H, F, t)
folding_ratio = len(data) / len(str(folded_data))

# Calculate information retention using a different approach
information_retention = np.abs(np.mean(data) - folded_data) / np.mean(data)

print("Folded data:", folded_data)
print("Folding ratio:", folding_ratio)
print("Information retention:", information_retention)