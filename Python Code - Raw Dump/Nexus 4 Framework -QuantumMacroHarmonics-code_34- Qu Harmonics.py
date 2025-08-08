import numpy as np
import matplotlib.pyplot as plt

# Harmonic transformation for binary data
def harmonic_transform(binary_data):
    """Transform binary data into harmonic space."""
    n = len(binary_data)
    x = np.linspace(0, 2 * np.pi, n)
    harmonics = np.sin(x) * (binary_data - 0.5)
    return harmonics

# Recursive harmonic echo cancellation
def cancel_harmonic_echoes(padded_harmonics, hashed_harmonics, iterations=256, alpha=1.5, decay_factor=0.9):
    """
    Iteratively cancel harmonic echoes (residual deviations) to refine alignment.
    """
    aligned_harmonics = hashed_harmonics.copy()
    echo_correction_history = []

    for _ in range(iterations):
        # Calculate phase and amplitude differences (echo-like residuals)
        phase_shift = np.roll(aligned_harmonics, shift=-1) - aligned_harmonics
        amplitude_correction = alpha * (padded_harmonics - aligned_harmonics)

        # Apply echo cancellation (simulate decay of recursive reflections)
        echo_correction = decay_factor * phase_shift

        # Combine corrections to refine alignment
        aligned_harmonics += amplitude_correction - echo_correction
        echo_correction_history.append(aligned_harmonics.copy())

    return aligned_harmonics, echo_correction_history

# Truncate or scale harmonics to match lengths
def truncate_or_scale_harmonics(padded_harmonics, target_length):
    """
    Truncate or scale padded harmonics to match the length of hashed harmonics.
    """
    current_length = len(padded_harmonics)
    if current_length > target_length:
        return padded_harmonics[:target_length]  # Truncate
    elif current_length < target_length:
        return np.interp(
            np.linspace(0, current_length - 1, target_length),
            np.arange(current_length),
            padded_harmonics
        )  # Scale
    return padded_harmonics

# Generate example binary data
def string_to_binary(string):
    """Convert a string to binary representation."""
    return ''.join(format(ord(c), '08b') for c in string)

def hash_to_binary(string):
    """Hash a string using SHA-256 and convert to binary representation."""
    import hashlib
    hashed = hashlib.sha256(string.encode()).hexdigest()
    return ''.join(format(int(c, 16), '04b') for c in hashed)

# Input: padded binary and hashed binary
original_binary = string_to_binary("abc") + '1' + '0' * (448 - len(string_to_binary("abc")) - 1) + format(len("abc") * 8, '064b')
padded_binary = np.array([int(b) for b in original_binary])
hashed_binary = np.array([int(b) for b in hash_to_binary("abc")])

# Generate harmonics
padded_harmonics_truncated = truncate_or_scale_harmonics(harmonic_transform(padded_binary), len(hashed_binary))
hashed_harmonics = harmonic_transform(hashed_binary)

# Apply echo cancellation for 256 iterations
echo_cancelled_harmonics, echo_cancellation_history = cancel_harmonic_echoes(
    padded_harmonics_truncated, hashed_harmonics, iterations=256
)

# Visualization of final alignment
plt.figure(figsize=(12, 6))
plt.plot(padded_harmonics_truncated, label="Padded Harmonics (Macro State)", color='blue')
plt.plot(hashed_harmonics, label="Hashed Harmonics (Initial)", color='orange', linestyle='dashed')
plt.plot(echo_cancelled_harmonics, label="Echo-Cancelled Harmonics (Final)", color='green', linestyle='dotted')
plt.title("Echo Cancellation in Harmonic Space: Padded vs. Hashed vs. Refined")
plt.xlabel("Iteration (n)")
plt.ylabel("Harmonic Amplitude")
plt.legend()
plt.grid()
plt.show()

# Visualization of progression
plt.figure(figsize=(12, 6))
for i in range(0, 256, 32):  # Plot every 32nd step for clarity
    plt.plot(echo_cancellation_history[i], label=f"Step {i+1}", linestyle='dotted')
plt.title("Echo Cancellation Progression")
plt.xlabel("Iteration (n)")
plt.ylabel("Harmonic Amplitude")
plt.legend()
plt.grid()
plt.show()
