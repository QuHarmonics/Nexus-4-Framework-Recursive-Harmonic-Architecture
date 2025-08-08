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

# Generate padded binary
def generate_padded_binary(input_string):
    binary_representation = ''.join(format(ord(c), '08b') for c in input_string)
    padding_length = 448 - len(binary_representation) - 1
    padded_binary = binary_representation + '1' + '0' * padding_length + format(len(input_string) * 8, '064b')
    return np.array([int(b) for b in padded_binary])

# Generate hashed binary
def hash_to_binary(string):
    hashed = hashlib.sha256(string.encode()).hexdigest()
    return ''.join(format(int(c, 16), '04b') for c in hashed)

# Strip padding to decode back to the original text
def remove_padding_and_decode(binary_data, original_length):
    # Convert binary array to a string
    binary_string = ''.join(map(str, binary_data))
    
    # Extract the original meaningful bits (exclude padding)
    meaningful_bits = binary_string[:original_length * 8]
    
    # Split into 8-bit chunks and convert to characters
    byte_chunks = [meaningful_bits[i:i + 8] for i in range(0, len(meaningful_bits), 8)]
    decoded_text = ''.join(chr(int(byte, 2)) for byte in byte_chunks)
    
    return decoded_text





# Input string
input_string = "abc"

# Generate padded binary and transform into harmonic space  ??????????wtf
padded_binary = generate_padded_binary(input_string)
padded_harmonics = binary_to_waveform(padded_binary)

# Generate hashed binary for visualization purposes
hashed_binary = np.array([int(b) for b in hash_to_binary(input_string)])
hashed_harmonics = binary_to_waveform(hashed_binary)

# Generate zeta zeros
num_zeros = 100  # Number of zeta zeros to use
zeta_zeros = generate_zeta_zeros(num_zeros)

# Apply zeta zero filtering to padded harmonics
filtered_harmonics = filter_with_zeta_zeros(padded_harmonics, zeta_zeros, len(padded_harmonics))

# Convert filtered harmonics back to binary
reconstructed_binary = waveform_to_binary(filtered_harmonics)

# Remove the first bit of the reconstructed binary
reconstructed_binary_trimmed = reconstructed_binary[1:]

# Convert reconstructed binary back to text
decoded_text = binary_to_text(reconstructed_binary_trimmed)

# Assuming the original text's length was 3 (for "abc")
original_length = len(input_string)

# Decode the reconstructed binary
decoded_text_from_reconstructed = remove_padding_and_decode(reconstructed_binary, original_length)



# Display results
print(f"Input String: {input_string}")
print(f"\nOrignal Padded Hash (First 128 Bits):\n{padded_binary[:128]}")
print(f"\nSHA-256 Hashed Binary (First 128 Bits):\n{hashed_binary[:128]}")
print(f"\nReconstructed Binary (First 128 Bits):\n{reconstructed_binary[:128]}")
print(f"\nDecoded Text from Reconstructed Binary:\n{decoded_text_from_reconstructed}")

# Plot the comparison
plt.figure(figsize=(12, 12))

# Initial padded harmonic waveform
plt.subplot(3, 1, 1)
plt.plot(padded_harmonics[:128], label="Padded Harmonic Waveform", color='blue')
plt.title("Padded Harmonic Waveform (First 128 Points)")
plt.ylabel("Amplitude")
plt.grid()
plt.legend()

# Hashed harmonic waveform
plt.subplot(3, 1, 2)
plt.plot(hashed_harmonics[:128], label="Hashed Harmonic Waveform (Visualization Only)", color='orange')
plt.title("Hashed Harmonic Waveform (Visualization Only, First 128 Points)")
plt.ylabel("Amplitude")
plt.grid()
plt.legend()

# Filtered harmonic waveform after zeta zeros
plt.subplot(3, 1, 3)
plt.plot(filtered_harmonics[:128], label="Filtered Harmonic Waveform", color='green', linestyle='--')
plt.title("Filtered Harmonic Waveform (After Zeta Zeros Applied, First 128 Points)")
plt.xlabel("Index")
plt.ylabel("Amplitude")
plt.grid()
plt.legend()

plt.tight_layout()
plt.show()
