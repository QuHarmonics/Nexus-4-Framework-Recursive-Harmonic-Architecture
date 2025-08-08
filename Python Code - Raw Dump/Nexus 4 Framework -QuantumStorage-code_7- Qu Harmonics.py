from tqdm import tqdm
import numpy as np

# Constants
HARMONIC_CONSTANT = 0.35
TOLERANCE = 0.01
MAX_ITERATIONS = 100

# 1. Read File as Raw Bytes
def read_file_as_bytes(filename):
    """Reads a file as raw bytes."""
    with open(filename, 'rb') as f:
        raw_data = f.read()
    print(f"File size: {len(raw_data)} bytes")
    return raw_data

# 2. Save File from Raw Bytes
def save_file_from_bytes(filename, raw_bytes):
    """Writes raw bytes to a file."""
    with open(filename, 'wb') as f:
        f.write(raw_bytes)
    print(f"File saved: {filename}, Size: {len(raw_bytes)} bytes")

# 3. Convert Bytes to Matrix
def bytes_to_matrix(raw_data):
    """Converts raw bytes into a 2D matrix."""
    binary_length = len(raw_data) * 8  # Total bits
    matrix_size = int(np.ceil(np.sqrt(binary_length)))  # Square matrix size
    data_matrix = np.zeros((matrix_size, matrix_size), dtype=np.float64)

    # Fill the matrix with binary data
    for i in tqdm(range(len(raw_data)), desc="Converting bytes to matrix"):
        byte = raw_data[i]
        for bit_index in range(8):
            bit = (byte >> (7 - bit_index)) & 1
            index = i * 8 + bit_index
            if index < binary_length:
                row = index // matrix_size
                col = index % matrix_size
                data_matrix[row, col] = bit

    return data_matrix, binary_length

# 4. Convert Matrix Back to Bytes
def matrix_to_bytes(expanded_matrix, original_length):
    """Converts a 2D matrix back into raw bytes."""
    flat_bits = []
    for row in expanded_matrix:
        for bit in row:
            flat_bits.append(int(round(bit)))

    # Trim to original bit length
    flat_bits = flat_bits[:original_length]
    raw_bytes = bytearray()

    # Group bits into bytes
    for i in range(0, len(flat_bits), 8):
        byte = 0
        for bit_index in range(8):
            if i + bit_index < len(flat_bits):
                byte |= flat_bits[i + bit_index] << (7 - bit_index)
        raw_bytes.append(byte)

    return bytes(raw_bytes)

# 5. Harmonize Data (Compression/Expansion)
def harmonize_data(data, harmonic_constant, compress=True):
    """Compress or expand data by aligning with the harmonic constant."""
    data = np.copy(data)
    gain = 1.0 if compress else -1.0

    for iteration in tqdm(range(MAX_ITERATIONS), desc="Harmonic adjustment"):
        delta = np.mean(data) - harmonic_constant
        adjustment = delta * gain
        data -= adjustment
        data = np.clip(data, 0, 1)  # Ensure values stay in binary range
        if abs(delta) < TOLERANCE:
            break

    return data

# Main Workflow
if __name__ == "__main__":
    # Input and Output Files
    input_file = "d:\\test.wav"  # Can be any file, not just .wav
    compressed_file = "d:\\compressed_data.npz"
    output_file = "d:\\restored_output.wav"

    # Step 1: Read Input File as Bytes
    raw_data = read_file_as_bytes(input_file)

    # Step 2: Convert Bytes to Matrix
    data_matrix, binary_length = bytes_to_matrix(raw_data)

    # Step 3: Apply Harmonic Compression
    compressed_data = harmonize_data(data_matrix, HARMONIC_CONSTANT, compress=True)
    print(f"Compressed Data Non-Zero Elements: {np.count_nonzero(compressed_data)}")

    # Step 4: Save Compressed Data
    np.savez_compressed(compressed_file, compressed_matrix=compressed_data)
    print(f"Compressed data saved to {compressed_file}")

    # Step 5: Load Compressed Data and Expand
    loaded_data = np.load(compressed_file)['compressed_matrix']
    expanded_data = harmonize_data(loaded_data, HARMONIC_CONSTANT, compress=False)

    # Step 6: Convert Matrix Back to Bytes
    restored_bytes = matrix_to_bytes(expanded_data, binary_length)

    # Step 7: Save Restored File
    save_file_from_bytes(output_file, restored_bytes)
