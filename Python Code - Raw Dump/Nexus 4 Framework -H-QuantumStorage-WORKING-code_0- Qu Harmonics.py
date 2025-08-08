import numpy as np

# Constants
HARMONIC_CONSTANT = 0.35
TOLERANCE = 0.01
MAX_ITERATIONS = 100

# Utility Functions
def text_to_base16(text):
    """Convert text to its base-16 representation."""
    return ''.join(format(ord(char), '02x') for char in text)

def base16_to_text(base16_str):
    """Convert base-16 encoded data back to text."""
    try:
        chars = [chr(int(base16_str[i:i + 2], 16)) for i in range(0, len(base16_str), 2)]
        return ''.join(chars)
    except Exception as e:
        return f"Error decoding base16: {e}"

def binary_to_text(binary_str):
    """Convert binary string back to text."""
    try:
        byte_chunks = [binary_str[i:i + 8] for i in range(0, len(binary_str), 8)]
        text = ''.join(chr(int(chunk, 2)) for chunk in byte_chunks if len(chunk) == 8)
        return text
    except Exception as e:
        return f"Error decoding binary: {e}"

# Harmonic Compression and Expansion
def harmonize_data(data, harmonic_constant, compress=True):
    """
    Harmonize data for compression or expansion.
    If compress=True, align data with harmonic resonance for compression.
    If compress=False, expand harmonized data back to its original form.
    """
    data = np.copy(data)
    gain = 1.0 if compress else -1.0

    for _ in range(MAX_ITERATIONS):
        delta = np.mean(data) - harmonic_constant
        adjustment = delta * gain
        data -= adjustment
        data = np.clip(data, 0, 1)  # Ensure binary range
        if abs(delta) < TOLERANCE:
            break

    return data

# Main Framework
def hybrid_storage(text):
    """Encodes, compresses, and expands text using the hybrid framework."""
    print("Original Text Length:", len(text))
    print("Original Text:", text)

    # Step 1: Base-16 Encoding
    base16_data = text_to_base16(text)
    print("Base-16 Encoded Data Length:", len(base16_data))
    
    # Step 2: Convert Base-16 Data to Binary
    binary_data = ''.join(format(int(char, 16), '04b') for char in base16_data)
    original_binary_length = len(binary_data)
    print("Binary Data Length:", original_binary_length)

    # Step 3: Prepare Data Matrix
    matrix_size = int(np.ceil(np.sqrt(original_binary_length)))
    data_matrix = np.zeros((matrix_size, matrix_size), dtype=np.float64)
    for i, bit in enumerate(binary_data):
        data_matrix[i // matrix_size, i % matrix_size] = int(bit)

    # Step 4: Harmonic Compression
    compressed_data = harmonize_data(data_matrix, HARMONIC_CONSTANT, compress=True)
    compressed_size = np.count_nonzero(compressed_data)
    print("Compressed Data Non-Zero Elements:", compressed_size)

    # Step 5: Harmonic Expansion
    expanded_data = harmonize_data(compressed_data, HARMONIC_CONSTANT, compress=False)
    flattened_data = ''.join(str(int(round(bit))) for bit in expanded_data.flatten()[:original_binary_length])

    # Step 6: Convert Binary Back to Text
    recovered_text = binary_to_text(flattened_data)
    print("Recovered Text:", recovered_text)

    # Compression Statistics
    compression_ratio = compressed_size / original_binary_length
    print("Compression Ratio (Compressed/Original):", compression_ratio)

    return recovered_text, compression_ratio

# Example Usage
if __name__ == "__main__":
    # Large Text Dataset
    input_text = "Quantum Storage Test " * 50  # Scaling input text
    recovered, compression_ratio = hybrid_storage(input_text)
    print("Final Recovered Text:", recovered[:100])  # Display first 100 characters
    print("Compression Ratio:", compression_ratio)
