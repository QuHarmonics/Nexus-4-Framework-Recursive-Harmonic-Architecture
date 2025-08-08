def compress_binary(binary_data):
    """
    Compress binary data using a combination of run-length encoding and base conversion.
    """
    # Run-Length Encoding
    compressed = []
    count = 1
    for i in range(1, len(binary_data)):
        if binary_data[i] == binary_data[i - 1]:
            count += 1
        else:
            compressed.append((binary_data[i - 1], count))
            count = 1
    compressed.append((binary_data[-1], count))  # Add the last group

    # Convert RLE to a more compact base representation
    compressed_string = "".join(f"{val}{count}" for val, count in compressed)
    return compressed_string


def decompress_binary(compressed_data):
    """
    Decompress binary data from a compressed string.
    """
    decompressed = []
    for i in range(0, len(compressed_data), 2):
        val = compressed_data[i]
        count = int(compressed_data[i + 1])
        decompressed.extend([val] * count)
    return "".join(decompressed)


# Example Data
binary_data = "110100101011001010101001" * 10  # Repeated pattern for testing
compressed = compress_binary(binary_data)
decompressed = decompress_binary(compressed)

print("Original Binary:", binary_data)
print("Compressed:", compressed)
print("Decompressed:", decompressed)
print("Matches Original:", decompressed == binary_data)
