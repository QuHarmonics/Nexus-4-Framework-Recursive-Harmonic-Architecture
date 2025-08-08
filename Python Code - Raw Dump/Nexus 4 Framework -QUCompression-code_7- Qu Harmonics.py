from tqdm import tqdm
import numpy as np

# Constants
HARMONIC_CONSTANT = 0.35
TOLERANCE = 0.01
MAX_ITERATIONS = 100

# Harmonic Compression and Expansion
def harmonize_data(data, harmonic_constant, compress=True):
    """Compress or expand data by aligning with the harmonic constant."""
    data = np.copy(data)
    gain = 1.0 if compress else -1.0

    for _ in tqdm(range(MAX_ITERATIONS), desc="Harmonic adjustment"):
        delta = np.mean(data) - harmonic_constant
        adjustment = delta * gain
        data -= adjustment
        data = np.clip(data, 0, 1)  # Ensure binary range
        if abs(delta) < TOLERANCE:
            break

    return data

# File Handling
def read_file_as_bytes(filename):
    """Reads a file as raw bytes."""
    with open(filename, 'rb') as f:
        raw_data = f.read()
    print(f"Read file: {filename}, Size: {len(raw_data)} bytes")
    return np.frombuffer(raw_data, dtype=np.uint8)  # Convert to NumPy array

def save_file_from_bytes(filename, raw_bytes):
    """Writes raw bytes to a file."""
    with open(filename, 'wb') as f:
        f.write(raw_bytes)
    print(f"File saved: {filename}, Size: {len(raw_bytes)} bytes")

# Compression Framework
def hybrid_audio_compression(input_file, compressed_file, output_file):
    """Compresses and decompresses a file using harmonic compression."""
    # Step 1: Read Input File
    raw_data = read_file_as_bytes(input_file)
    binary_data = ''.join(format(byte, '08b') for byte in tqdm(raw_data, desc="Converting bytes to binary"))
    original_binary_length = len(binary_data)

    # Step 2: Prepare Data Matrix
    matrix_size = int(np.ceil(np.sqrt(original_binary_length)))
    data_matrix = np.zeros((matrix_size, matrix_size), dtype=np.float64)
    for i, bit in tqdm(enumerate(binary_data), total=original_binary_length, desc="Filling matrix"):
        data_matrix[i // matrix_size, i % matrix_size] = int(bit)

    # Step 3: Harmonic Compression
    compressed_data = harmonize_data(data_matrix, HARMONIC_CONSTANT, compress=True)
    compressed_size = np.count_nonzero(compressed_data)
    print(f"Compressed Data Non-Zero Elements: {compressed_size}")

    # Step 4: Save Compressed Data
    np.savez_compressed(compressed_file, compressed_matrix=compressed_data)
    print(f"Compressed data saved to {compressed_file}")

    # Step 5: Load and Expand
    loaded_data = np.load(compressed_file)['compressed_matrix']
    expanded_data = harmonize_data(loaded_data, HARMONIC_CONSTANT, compress=False)

    # Step 6: Convert Back to Binary
    flattened_data = ''.join(str(int(round(bit))) for bit in tqdm(expanded_data.flatten()[:original_binary_length], desc="Flattening matrix"))
    restored_bytes = bytearray(int(flattened_data[i:i+8], 2) for i in range(0, len(flattened_data), 8))

    # Step 7: Save Restored File
    save_file_from_bytes(output_file, restored_bytes)

    # Step 8: Compression Statistics
    compression_ratio = compressed_size / original_binary_length
    print(f"Compression Ratio (Compressed/Original): {compression_ratio}")
    print(f"Is Data Lossless? {'Yes' if raw_data.tobytes() == restored_bytes else 'No'}")

# Main Execution
if __name__ == "__main__":
    input_file = "d:\\test.wav"  # Input file path
    compressed_file = "d:\\compressed_data.npz"  # Compressed output
    output_file = "d:\\restored_output.wav"  # Restored file path

    hybrid_audio_compression(input_file, compressed_file, output_file)
