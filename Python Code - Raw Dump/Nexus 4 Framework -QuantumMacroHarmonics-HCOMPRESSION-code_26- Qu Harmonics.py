import numpy as np
import matplotlib.pyplot as plt
import hashlib

# Convert a binary string to harmonic space
def harmonic_transform(binary_data):
    n = len(binary_data)
    x = np.linspace(0, 2 * np.pi, n)
    return np.sin(x) * (binary_data - 0.5)

# Recursive harmonic echo cancellation
def trinity_harmonic_cancellation(padded_harmonics, hashed_harmonics, iterations, alpha=1.5, decay_factor=0.9):
    aligned_harmonics = hashed_harmonics.copy()
    trinity_history = []
    for _ in range(iterations):
        phase_shift = np.roll(aligned_harmonics, shift=-1) - aligned_harmonics
        amplitude_correction = alpha * (padded_harmonics - aligned_harmonics)
        echo_correction = decay_factor * phase_shift
        aligned_harmonics += amplitude_correction - echo_correction
        trinity_history.append(aligned_harmonics.copy())
    return aligned_harmonics, trinity_history

# Detect bit flips
def detect_bit_flips(padded_binary, reconstructed_binary):
    flip_mask = padded_binary != reconstructed_binary
    bit_flips = np.sum(flip_mask)
    return flip_mask, bit_flips

# Apply bit flip corrections
def apply_flip_correction(reconstructed_binary, flip_mask):
    corrected_binary = reconstructed_binary.copy()
    corrected_binary[flip_mask] = 1 - corrected_binary[flip_mask]
    return corrected_binary

# Convert harmonics back to binary
def harmonics_to_binary(harmonics, threshold=0.0):
    return np.array([1 if h > threshold else 0 for h in harmonics], dtype=int)

# Match harmonic lengths by interpolation
def scale_harmonics(harmonics, target_length):
    current_length = len(harmonics)
    return np.interp(np.linspace(0, current_length - 1, target_length), np.arange(current_length), harmonics)

# Convert binary to text
def binary_to_text(binary_data):
    """
    Convert binary data back to text.
    Assumes binary data is grouped into 8-bit ASCII characters.
    """
    binary_string = ''.join(str(bit) for bit in binary_data)
    byte_chunks = [binary_string[i:i + 8] for i in range(0, len(binary_string), 8)]
    text = ''.join(chr(int(byte, 2)) for byte in byte_chunks if len(byte) == 8)
    return text

# Generate binary and hash
def string_to_binary(string):
    return ''.join(format(ord(c), '08b') for c in string)

def hash_to_binary(string):
    hashed = hashlib.sha256(string.encode()).hexdigest()
    return ''.join(format(int(c, 16), '04b') for c in hashed)

# Input string for hashing
input_string = "abc"
print(f"Input String: {input_string}")

# Padded and hashed binary data
original_binary = string_to_binary(input_string) + '1' + '0' * (448 - len(string_to_binary(input_string)) - 1) + format(len(input_string) * 8, '064b')
padded_binary = np.array([int(b) for b in original_binary])
hashed_binary = np.array([int(b) for b in hash_to_binary(input_string)])

# Transform into harmonic space
padded_harmonics = scale_harmonics(harmonic_transform(padded_binary), 1024)
hashed_harmonics = scale_harmonics(harmonic_transform(hashed_binary), 1024)

# Perform harmonic echo cancellation
iterations = 1024
final_harmonics, trinity_history = trinity_harmonic_cancellation(padded_harmonics, hashed_harmonics, iterations)

# Convert harmonics back to binary
reconstructed_binary = harmonics_to_binary(final_harmonics)

# Rescale padded binary to match the harmonic sample rate
padded_binary_scaled = scale_harmonics(padded_binary, 1024).round().astype(int)

# Detect and correct bit flips
flip_mask, bit_flips = detect_bit_flips(padded_binary_scaled, reconstructed_binary)
corrected_binary = apply_flip_correction(reconstructed_binary, flip_mask)

# Calculate differences
total_differences = np.sum(padded_binary_scaled != corrected_binary)

# Decode the padded binary, reconstructed binary, and SHA-256 hash
padded_text = binary_to_text(padded_binary)
reconstructed_text = binary_to_text(corrected_binary)
sha256_text = binary_to_text(hashed_binary)

# Display results
print("\nSHA-256 Hexadecimal Hash:")
print(hashlib.sha256(input_string.encode()).hexdigest())
print("\nSHA-256 Binary Hash:")
print(hashed_binary)
print("\Original Text with padding in binary:")
print(original_binary)
print("\nDecoded Binary:")
print(corrected_binary)


# Visualization of binary comparison
plt.figure(figsize=(12, 18))
plt.subplot(3, 1, 1)
plt.plot(padded_binary_scaled[:512], label="Padded Binary (Scaled)", color='blue', linestyle='-')
plt.title("Padded Binary")
plt.grid()
plt.legend()

plt.subplot(3, 1, 2)
plt.plot(hashed_binary[:512], label="SHA-256 Binary Output", color='orange', linestyle='-')
plt.title("SHA-256 Binary Output")
plt.grid()

plt.subplot(3, 1, 3)
plt.plot(corrected_binary[:512], label="Corrected Binary", color='green', linestyle='-')
plt.title("Corrected Binary")
plt.grid()
plt.legend()


plt.legend()

plt.tight_layout()
plt.show()
