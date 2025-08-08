def binary_to_ascii_sliding(binary_string, bit_len=8, step=5):
    """
    Slides through the binary string in steps and extracts ASCII characters from each window.

    Parameters:
    - binary_string (str): Input binary string.
    - bit_len (int): Number of bits to interpret per character (usually 8).
    - step (int): Number of bits to move forward after each window.

    Returns:
    - List of tuples: (position, binary_segment, ASCII char or '.') for each parsed byte.
    """
    results = []
    for i in range(0, len(binary_string) - bit_len + 1, step):
        segment = binary_string[i:i+bit_len]
        try:
            char = chr(int(segment, 2))
            ascii_char = char if 32 <= ord(char) <= 126 else '.'
        except ValueError:
            ascii_char = '.'
        results.append((i, segment, ascii_char))
    
    return results

# Example usage
binary_data = '10111100000010011001001111001011111111101010111001010111001110101111001000010000110110100100101110011111001010010100000011111101101101010000110111100111010110100010101110110100111100011111000010111010000100111011011100111010010110110110010100100011000000010001101011010010110101001011010011111111101101001010100111110011001001101000000011001101011110010100100111001001110111011110010010100011000000001011000001100101101011000'  # binary for 'hello'
result = binary_to_ascii_sliding(binary_data, bit_len=8, step=5)
for position, segment, char in result:
    print(f"Pos {position:>3}: {segment} -> {char}")
