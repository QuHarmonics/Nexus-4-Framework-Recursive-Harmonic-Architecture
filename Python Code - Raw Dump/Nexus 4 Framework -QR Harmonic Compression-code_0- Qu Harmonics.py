import numpy as np
import matplotlib.pyplot as plt

# --- Constants ---
QBIT_SIZE = 64
HARMONIC_CONSTANT = 0.35
K = 0.15  # Stabilization constant
ITERATIONS = 256  # Recursive refinement iterations

# --- Step 1: Quantum Fourier Transform (QFT) ---
def qft(data):
    n = len(data)
    result = np.zeros(n, dtype=complex)
    for k in range(n):
        for j in range(n):
            result[k] += data[j] * np.exp(2j * np.pi * k * j / n)
    return result / np.sqrt(n)

# --- Step 2: Samsons Law Recursive Stabilization ---
def samsons_law(data, k):
    mixed_data = np.zeros(len(data), dtype=complex)
    for i in range(1, len(data)):
        mixed_data[i] = data[i] + k * (data[i] - data[i-1])
    mixed_data[0] = data[0]  # Boundary condition
    return mixed_data

def recursive_refinement(data, iterations, k):
    for _ in range(iterations):
        data = samsons_law(data, k)
    return data

# --- Step 3: Embed into 64x64x64 Qbit Matrix ---
def embed_into_qbit_matrix(data):
    padded = np.zeros(QBIT_SIZE**3)
    length = min(len(data), QBIT_SIZE**3)
    padded[:length] = data[:length]
    return padded.reshape((QBIT_SIZE, QBIT_SIZE, QBIT_SIZE))

# --- Step 4: Collapse from Qbit Matrix back to 1D ---
def collapse_from_qbit_matrix(qmatrix):
    return qmatrix.flatten()

# --- Generate Random Test Data ---
raw_data = np.random.rand(10000)

# --- Apply QFT and Recursive Refinement ---
qft_data = qft(raw_data)
refined_data = recursive_refinement(qft_data, ITERATIONS, K)
reconstructed_data = np.real(refined_data)

# --- Qbit Matrix Embedding and Collapse ---
qbit_matrix = embed_into_qbit_matrix(reconstructed_data)
collapsed_data = collapse_from_qbit_matrix(qbit_matrix)

# --- Visualization of a slice ---
slice_projection = qbit_matrix[:, :, QBIT_SIZE // 2]
plt.figure(figsize=(6, 5))
plt.imshow(slice_projection, cmap='viridis')
plt.title("Middle Z-Slice of Qbit Matrix")
plt.colorbar(label="Amplitude")
plt.tight_layout()
plt.show()

# Also plot reconstructed waveform
plt.figure(figsize=(5, 4))
plt.plot(reconstructed_data[:1000])
plt.title("Reconstructed Data (First 100)")
plt.xlabel("Index")
plt.ylabel("Value")
plt.tight_layout()
plt.show()
