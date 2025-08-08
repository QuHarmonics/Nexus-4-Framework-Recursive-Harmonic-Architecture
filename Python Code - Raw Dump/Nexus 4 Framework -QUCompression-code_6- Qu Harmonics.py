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
    return np.frombuffer(raw_data, dtype=np.uint8)  # Convert to NumPy array

# 2. Save File from Raw Bytes
def save_file_from_bytes(filename, raw_bytes):
    """Writes raw bytes to a file."""
    with open(filename, 'wb') as f:
        f.write(raw_bytes)
    print(f"File saved: {filename}, Size: {len(raw_bytes)} bytes")

# 3. Difference Encoding with Progress Meter
def difference_encode(audio_data):
    """Encodes differences between consecutive samples with progress tracking."""
    diff_data = np.empty_like(audio_data, dtype=np.int16)
    diff_data[0] = audio_data[0]
    for i in tqdm(range(1, len(audio_data)), desc="Performing Difference Encoding"):
        diff_data[i] = audio_data[i] - audio_data[i - 1]
    return diff_data

# 4. Difference Decoding with Progress Meter
def difference_decode(encoded_data):
    """Decodes difference-encoded data back to original form with progress tracking."""
    restored_data = np.empty_like(encoded_data, dtype=np.int16)
    restored_data[0] = encoded_data[0]
    for i in tqdm(range(1, len(encoded_data)), desc="Performing Difference Decoding"):
        restored_data[i] = restored_data[i - 1] + encoded_data[i]
    return restored_data

# 5. Fourier Transform Compression
def harmonic_compress(audio_data):
    """Compress audio data using Fourier Transform."""
    print("Performing Fourier Compression...")
    transformed = rfft(audio_data)  # Real-to-complex Fourier Transform
    return transformed

def harmonic_expand(compressed_data):
    """Expand compressed data using Inverse Fourier Transform."""
    print("Performing Fourier Expansion...")
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
    print("Saving compressed data...")
    np.savez_compressed(compressed_file, compressed_data=compressed_data)
    print(f"Compressed data saved to {compressed_file}")

    # Load and Expand
    print("Loading compressed data...")
    loaded_data = np.load(compressed_file)['compressed_data']
    expanded_diff = harmonic_expand(loaded_data)

    # Step 4: Decode Differences
    restored_data = difference_decode(expanded_diff)

    # Step 5: Save Restored File
    save_file_from_bytes(output_file, restored_data.astype(np.uint8))
