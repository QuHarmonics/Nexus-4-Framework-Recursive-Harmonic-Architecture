from tqdm import tqdm
import numpy as np
from scipy.fft import rfft, irfft

# Constants
HARMONIC_CONSTANT = 0.35
TOLERANCE = 0.01
MAX_ITERATIONS = 100

# 1. Read File as Bytes
def read_file_as_bytes(filename):
    """Reads a file as raw bytes."""
    with open(filename, 'rb') as f:
        raw_data = f.read()
    print(f"File size: {len(raw_data)} bytes")
    return raw_data

# 2. Difference Encoding
def difference_encode(audio_data):
    """Encodes differences between consecutive samples."""
    return np.diff(audio_data, prepend=audio_data[0])

def difference_decode(encoded_data):
    """Decodes difference-encoded data back to original form."""
    return np.cumsum(encoded_data)

# 3. Sparse Encoding
def sparse_encode(audio_data):
    """Encodes sparse data by storing non-zero values and their indices."""
    indices = np.nonzero(audio_data)[0]
    values = audio_data[indices]
    return indices, values

def sparse_decode(indices, values, original_length):
    """Decodes sparse data back to the original form."""
    decoded = np.zeros(original_length, dtype=values.dtype)
    decoded[indices] = values
    return decoded

# 4. Fourier Transform Compression
def harmonic_compress(audio_data):
    """Compress audio data using Fourier Transform."""
    transformed = rfft(audio_data)  # Real-to-complex Fourier Transform
    return transformed

def harmonic_expand(compressed_data):
    """Expand compressed data using Inverse Fourier Transform."""
    restored = irfft(compressed_data)
    return np.round(restored).astype(np.int16)  # Round to integers

# Main Workflow
if __name__ == "__main__":
    # Input and Output Files
    input_file = "d:\\test.wav"  # Can be any file, not just .wav
    compressed_file = "d:\\compressed_data.npz"
    output_file = "d:\\restored_output.wav"

    # Step 1: Read Input File as Bytes
    raw_data = read_file_as_bytes(input_file)

    # Step 2: Difference Encode
    encoded_diff = difference_encode(raw_data)

    # Step 3: Apply Fourier Transform Compression
    compressed_data = harmonic_compress(encoded_diff)

    # Save Compressed Data
    np.savez_compressed(compressed_file, compressed_data=compressed_data)
    print(f"Compressed data saved to {compressed_file}")

    # Load and Expand
    loaded_data = np.load(compressed_file)['compressed_data']
    expanded_diff = harmonic_expand(loaded_data)

    # Step 4: Decode Differences
    restored_data = difference_decode(expanded_diff)

    # Step 5: Save Restored File
    save_file_from_bytes(output_file, restored_data)
