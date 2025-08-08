def binary_to_ascii(binary_string):
    # Ensure the binary string length is a multiple of 8
    length = len(binary_string)
    if length % 8 != 0:
        binary_string = binary_string[:length - (length % 2)]
    
    # Split the binary string into 8-bit chunks
    byte_chunks = [binary_string[i:i+8] for i in range(0, len(binary_string), 8)]
    
    # Convert each chunk to an ASCII character (or '.' if non-printable)
    ascii_output = ''
    for byte in byte_chunks:
        try:
            char = chr(int(byte, 2))
            ascii_output += char if 32 <= ord(char) <= 126 else '.'
        except ValueError:
            ascii_output += '.'
    
    return ascii_output

# Example usage
binary_data = '10111100000010011001001111001011111111101010111001010111001110101111001000010000110110100100101110011111001010010100000011111101101101010000110111100111010110100010101110110100111100011111000010111010000100111011011100111010010110110110010100100011000000010001101011010010110101001011010011111111101101001010100111110011001001101000000011001101011110010100100111001001110111011110010010100011000000001011000001100101101011000'  # binary for 'hello'
result = binary_to_ascii(binary_data)
print(result)
