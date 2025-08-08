def binary_to_ascii_sliding_debug(binary_string, bit_len=8, step=5, watch_chars=None):
    """
    Slides through the binary string in fixed-length windows and decodes ASCII characters.

    Parameters:
    - binary_string (str): The input binary stream.
    - bit_len (int): Bit size per chunk (e.g., 8 for one byte).
    - step (int): Step size to slide the window.
    - watch_chars (set): Optional set of characters to highlight during search.

    Returns:
    - None. Prints each match.
    """
    if watch_chars is None:
        watch_chars = set()

    for i in range(0, len(binary_string) - bit_len + 1, step):
        segment = binary_string[i:i+bit_len]
        try:
            val = int(segment, 2)
            char = chr(val)
            printable = char if 32 <= val <= 126 else '.'
        except:
            printable = '.'

        marker = "<-- match" if printable in watch_chars else ""
        print(f"Pos {i:>4}: {segment} -> {printable} {marker}")

# Example usage
binstream = "10111100000010011001001111001011111111101010111001010111001110101111001000010000110110100100101110011111001010010100000011111101101101010000110111100111010110100010101110110100111100011111000010111010000100111011011100111010010110110110010100100011000000010001101011010010110101001011010011111111101101001010100111110011001001101000000011001101011110010100100111001001110111011110010010100011000000001011000001100101101011000"
binary_to_ascii_sliding_debug(binstream, bit_len=8, step=1, watch_chars={'l'})
