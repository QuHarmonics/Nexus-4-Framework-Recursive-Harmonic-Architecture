import numpy as np
import matplotlib.pyplot as plt
from mpmath import zetazero
import hashlib

# Generate zeta zeros
def generate_zeta_zeros(num_zeros):
    return [zetazero(n).real for n in range(1, num_zeros + 1)]

# Transform binary data into harmonic space
def binary_to_waveform(binary_data):
    n = len(binary_data)
    x = np.linspace(0, 2 * np.pi, n)
    return np.sin(x) * (binary_data - 0.5)

# Overlay zeta zeros on harmonic space
def filter_with_zeta_zeros(waveform, zeta_zeros, harmonic_length):
    filtered_waveform = waveform.copy()
    for zero in zeta_zeros:
        zero_index = int((zero / max(zeta_zeros)) * harmonic_length)
        if 0 <= zero_index < len(filtered_waveform):
            filtered_waveform[zero_index] = 0  # Cancel predicted noise at zeta zero
    return filtered_waveform

# Convert harmonic space back to binary
def waveform_to_binary(waveform, threshold=0.0):
    return np.array([1 if h > threshold else 0 for h in waveform], dtype=int)

# Convert binary to text by grouping into bytes
def binary_to_text(binary_data):
    binary_string = ''.join(map(str, binary_data))
    byte_chunks = [binary_string[i:i + 8] for i in range(0, len(binary_string), 8)]
    text = ''.join(chr(int(byte, 2)) for byte in byte_chunks if len(byte) == 8)
    return text

# Generate hashed binary
def hash_to_binary(string):
    hashed = hashlib.sha256(string.encode()).hexdigest()
    return ''.join(format(int(c, 16), '04b') for c in hashed)

# Input string and hash
input_string = "abc"
hashed_binary = np.array([int(b) for b in hash_to_binary(input_string)])

# Transform the hash into harmonic space
hashed_waveform = binary_to_waveform(hashed_binary)

# Generate zeta zeros
num_zeros = 100  # Number of zeta zeros to use
zeta_zeros = generate_zeta_zeros(num_zeros)

# Filter the hashed waveform using zeta zeros
filtered_harmonics = filter_with_zeta_zeros(hashed_waveform, zeta_zeros, len(hashed_waveform))

# Convert the filtered harmonics back to binary
reconstructed_binary = waveform_to_binary(filtered_harmonics)

# Remove the first bit of the reconstructed binary (assumed padding bit)
reconstructed_binary_trimmed = reconstructed_binary[1:]

# Attempt to decode the reconstructed binary into text
decoded_text = binary_to_text(reconstructed_binary_trimmed)

# Display results
print(f"Input String: {input_string}")
print(f"\nSHA-256 Hashed Binary (First 128 Bits):\n{hashed_binary[:128]}")
print(f"\nReconstructed Binary (First 128 Bits, After Zeta Filtering):\n{reconstructed_binary[:128]}")
print(f"\nDecoded Text from Reconstructed Binary:\n{decoded_text}")

# Plot the comparison
plt.figure(figsize=(12, 10))

# Original hashed binary waveform
plt.subplot(3, 1, 1)
plt.plot(hashed_waveform[:128], label="Hashed Binary Waveform", color='orange')
plt.title("Hashed Binary Waveform (First 128 Points)")
plt.grid()
plt.legend()

# Filtered harmonic waveform
plt.subplot(3, 1, 2)
plt.plot(filtered_harmonics[:128], label="Filtered Harmonic Waveform", color='green')
plt.title("Filtered Harmonic Waveform After Zeta Zeros (First 128 Points)")
plt.grid()
plt.legend()

# Reconstructed binary waveform
plt.subplot(3, 1, 3)
plt.plot(reconstructed_binary[:128], label="Reconstructed Binary Waveform", color='blue')
plt.title("Reconstructed Binary Waveform (First 128 Points)")
plt.grid()
plt.legend()

plt.tight_layout()
plt.show()
