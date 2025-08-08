import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from scipy.stats import entropy

# Function to generate binary representation of Pi digits
def generate_pi_binary(n_digits=1000):
    from mpmath import mp
    mp.dps = n_digits + 2  # Extra precision to avoid rounding issues
    pi_str = str(mp.pi)[2:]  # Get pi digits after decimal

    # Convert each digit into a 4-bit binary representation
    pi_binary = [format(int(digit), '04b') for digit in pi_str[:n_digits]]
    return pi_binary

# Function to compute bit-pair ratio
def compute_bit_pair_ratio(binary_seq):
    bit_pairs = [b[i:i+2] for b in binary_seq for i in range(len(b)-1)]
    total_pairs = len(bit_pairs)
    count_11 = sum(1 for pair in bit_pairs if pair == '11')
    return count_11 / total_pairs if total_pairs > 0 else 0

# Generate binary representation of Pi digits
n_digits = 1000  # Number of digits to analyze
pi_binary = generate_pi_binary(n_digits)

# Compute bit-pair ratios for different lengths (14, 15, ...)
bit_pair_ratios = [compute_bit_pair_ratio(pi_binary[:length]) for length in range(14, 30)]

# Apply Fourier Transform to the bit-pair ratio sequence
fft_values = np.abs(fft(bit_pair_ratios))

# Perform PCA and t-SNE for clustering analysis
binary_vectors = np.array([[int(bit) for bit in ''.join(b)] for b in pi_binary[:500]])  # Flattened bit vectors
pca = PCA(n_components=2).fit_transform(binary_vectors)
tsne = TSNE(n_components=2, perplexity=30).fit_transform(binary_vectors)

# Compute entropy of the binary sequence
binary_sequence = ''.join(pi_binary[:n_digits])  # Flatten binary sequence
bit_probs = [binary_sequence.count(bit) / len(binary_sequence) for bit in ['0', '1']]
entropy_value = entropy(bit_probs, base=2)

# Plot Results
fig, axs = plt.subplots(3, 1, figsize=(10, 15))

# Bit-Pair Ratios
axs[0].plot(range(14, 30), bit_pair_ratios, marker='o', linestyle='-')
axs[0].set_title("Bit-Pair Ratio of Pi Binary at Different Lengths")
axs[0].set_xlabel("Length (Digits)")
axs[0].set_ylabel("Bit-Pair Ratio")

# Fourier Transform Spectrum
axs[1].plot(fft_values[:20], marker='x')
axs[1].set_title("Fourier Transform of Bit-Pair Ratio Sequence")
axs[1].set_xlabel("Frequency Component")
axs[1].set_ylabel("Magnitude")

# PCA and t-SNE Scatter Plots
fig_pca = plt.figure(figsize=(6, 6))
plt.scatter(pca[:, 0], pca[:, 1], c='blue', alpha=0.5)
plt.title("PCA Clustering of Pi Binary Representation")

fig_tsne = plt.figure(figsize=(6, 6))
plt.scatter(tsne[:, 0], tsne[:, 1], c='red', alpha=0.5)
plt.title("t-SNE Clustering of Pi Binary Representation")

plt.show()

# Display entropy result
entropy_value
