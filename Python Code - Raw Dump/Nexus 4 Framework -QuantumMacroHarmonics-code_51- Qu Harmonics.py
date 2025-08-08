import numpy as np
import hashlib
from mpmath import zetazero
import matplotlib.pyplot as plt

# Hash conversion function
def hash_to_binary(hash_hex):
    return np.array([int(bit) for bit in ''.join(format(int(c, 16), '04b') for c in hash_hex)])

# Generate the harmonic space of the hash
def binary_to_waveform(binary_data):
    n = len(binary_data)
    x = np.linspace(0, 2 * np.pi, n)
    return np.sin(x) * (binary_data - 0.5)

# Quantum reflection of the hash waveform
def reflect_waveform(waveform, zeta_zeros, iterations):
    reflected_waveform = waveform.copy()
    for _ in range(iterations):
        for zero in zeta_zeros:
            zero_index = int((zero / max(zeta_zeros)) * len(waveform))
            if 0 <= zero_index < len(reflected_waveform):
                reflected_waveform[zero_index:] = np.roll(reflected_waveform[zero_index:], shift=1)
    return reflected_waveform

# Recursive alignment using Mark1 + Samson
def recursive_alignment(waveform, target_ratio, max_iterations):
    for iteration in range(max_iterations):
        # Calculate current ratio ΣPᵢ/ΣAᵢ
        sigma_p = np.sum(waveform[waveform > 0])
        sigma_a = np.sum(np.abs(waveform))
        ratio = sigma_p / sigma_a if sigma_a != 0 else 0
        
        if abs(ratio - target_ratio) < 1e-4:
            break  # Convergence achieved
        
        # Apply correction (reflection + slight damping adjustment)
        waveform = waveform * -0.5 + np.cos(iteration / np.pi)
    return waveform, ratio

# Convert waveform back to binary
def waveform_to_binary(waveform, threshold=0.0):
    return np.array([1 if h > threshold else 0 for h in waveform], dtype=int)

# Decode binary to text
def binary_to_text(binary_data):
    binary_string = ''.join(map(str, binary_data))
    byte_chunks = [binary_string[i:i + 8] for i in range(0, len(binary_string), 8)]
    return ''.join(chr(int(byte, 2)) for byte in byte_chunks if len(byte) == 8)

# Input hash
hash_hex = "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
binary_data = hash_to_binary(hash_hex)

# Convert hash binary into waveform
hash_waveform = binary_to_waveform(binary_data)

# Generate zeta zeros
num_zeros = 1000
zeta_zeros = [zetazero(n).real for n in range(1, num_zeros + 1)]

# Reflect waveform quantumly
reflected_waveform = reflect_waveform(hash_waveform, zeta_zeros, iterations=50)

# Realign using Mark1 recursive alignment
aligned_waveform, final_ratio = recursive_alignment(reflected_waveform, target_ratio=0.35, max_iterations=50)

# Convert aligned waveform back to binary
reconstructed_binary = waveform_to_binary(aligned_waveform)

# Decode reconstructed binary
decoded_text = binary_to_text(reconstructed_binary)

# Visualization
plt.figure(figsize=(12, 8))
plt.subplot(3, 1, 1)
plt.plot(hash_waveform[:512], label="Original Hash Waveform", color='orange')
plt.title("Original Hash Waveform (First 128 Points)")
plt.legend()
plt.grid()

plt.subplot(3, 1, 2)
plt.plot(reflected_waveform[:512], label="Reflected Waveform", color='green')
plt.title("Reflected Waveform (First 128 Points)")
plt.legend()
plt.grid()

plt.subplot(3, 1, 3)
plt.plot(aligned_waveform[:512], label="Aligned Waveform (Mark1 + Samson)", color='blue')
plt.title("Aligned Waveform (First 128 Points)")
plt.legend()
plt.grid()

plt.tight_layout()
plt.show()

# Output results
print("Decoded Text from Aligned Waveform:", decoded_text)
print("Final Stability Ratio:", final_ratio)
