def compress_with_pattern_recognition(binary_data):
    """
    Compress binary data by recognizing repeated patterns.
    """
    pattern_map = {}
    compressed_data = []
    index = 0

    while index < len(binary_data):
        pattern = binary_data[index:index + 24]  # Search for 24-bit patterns
        if pattern not in pattern_map:
            pattern_map[pattern] = len(pattern_map)  # Assign a unique ID to the pattern
        compressed_data.append(pattern_map[pattern])
        index += 24

    return compressed_data, pattern_map

def decompress_with_pattern_recognition(compressed_data, pattern_map):
    """
    Decompress binary data using pattern recognition.
    """
    reversed_map = {v: k for k, v in pattern_map.items()}
    binary_data = ''.join(reversed_map[id] for id in compressed_data)
    return binary_data

# Example
original_binary = '110100101011001010101001' * 10
compressed_data, pattern_map = compress_with_pattern_recognition(original_binary)
decompressed_binary = decompress_with_pattern_recognition(compressed_data, pattern_map)

print("Original Binary:", original_binary)
print("Compressed Data:", compressed_data)
print("Decompressed Binary:", decompressed_binary)
print("Matches Original:", decompressed_binary == original_binary)
