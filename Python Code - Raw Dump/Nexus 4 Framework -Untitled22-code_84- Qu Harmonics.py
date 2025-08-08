import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fftn, ifftn

# Generate synthetic data representing Value (X), Binary Length (Y), and Bit-Pair Ratio (Z)
def generate_data(size=64):
    values = np.random.randint(1, 256, size)  # Random integer values
    binaries = [bin(v)[2:] for v in values]  # Convert to binary representation
    bin_lengths = np.array([len(b) for b in binaries])  # Binary length
    
    # Compute bit-pair ratio (count '11' pairs / total pairs)
    bit_ratios = np.array([b.count('11') / (len(b) - 1) if len(b) > 1 else 0 for b in binaries])

    return np.array(values), np.array(bin_lengths), np.array(bit_ratios)

# Generate dataset
values, bin_lengths, bit_ratios = generate_data()

# Normalize data for processing
X = values / np.max(values)  # Normalize value axis
Y = bin_lengths / np.max(bin_lengths)  # Normalize binary length axis
Z = bit_ratios  # Bit-pair ratio already scaled between 0-1

# Construct a 3D data tensor
data_3d = np.array([X, Y, Z]).T  # Shape (N, 3)

# Apply 3D Fourier Transform for frequency analysis
fft_result = fftn(data_3d)

# Inverse FFT to recover (validation check)
ifft_result = ifftn(fft_result)

# Visualization
fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')

ax.scatter(X, Y, Z, c='blue', marker='o', alpha=0.7)
ax.set_xlabel('Value (Normalized)')
ax.set_ylabel('Binary Length (Normalized)')
ax.set_zlabel('Bit-Pair Ratio')

plt.title("3D Representation of Multi-Dimensional Encoding")
plt.show()

# Output some key frequency components
print("Fourier Transform Magnitudes (First 10):")
print(np.abs(fft_result[:10]))

# Quantum FFT Placeholder (Simulating compression)
def qfft_placeholder(data):
    return np.fft.fft(data) / np.sqrt(len(data))  # Approximation of qFFT normalization

qfft_data = qfft_placeholder(data_3d)

print("\nQuantum FFT Approximation (First 10 Components):")
print(qfft_data[:10])
