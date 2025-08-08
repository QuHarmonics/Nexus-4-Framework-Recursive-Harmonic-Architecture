def binary_to_ascii_sliding_debug(binary_string, bit_len=8, step=5, watch_chars={'h', 'e', 'l', 'o'}):
    """
    Slides through the binary string in fixed-length windows, converts to ASCII, and highlights matches.

    Parameters:
    - binary_string (str): The input binary stream.
    - bit_len (int): Bit length of each character window (typically 8).
    - step (int): Bitwise step size between windows.
    - watch_chars (set): Characters to highlight when found.

    Output:
    - Prints index, binary segment, ASCII character, and marker if it's in watch_chars.
    """
    print(f"{'Index':>6} {'Binary':>12} {'Char':>6}   Match")
    print("-" * 36)
    for i in range(0, len(binary_string) - bit_len + 1, step):
        segment = binary_string[i:i + bit_len]
        try:
            value = int(segment, 2)
            char = chr(value)
            display_char = char if 32 <= value <= 126 else '.'
        except:
            display_char = '.'

        marker = "<-- match" if display_char in watch_chars else ""
        print(f"{i:6} {segment:>12} {display_char:>6}   {marker}")


# Example usage
if __name__ == "__main__":
    binary_input = (
        "10111100000010011001001111001011111111101010111001010111001110101111001000010000110110100100101110011111001010010100000011111101101101010000110111100111010110100010101110110100111100011111000010111010000100111011011100111010010110110110010100100011000000010001101011010010110101001011010011111111101101001010100111110011001001101000000011001101011110010100100111001001110111011110010010100011000000001011000001100101101011000"
    )

    binary_to_ascii_sliding_debug(
        binary_string=binary_input,
        bit_len=8,
        step=1,
        watch_chars={'h', 'e', 'l', 'o'}
    )
