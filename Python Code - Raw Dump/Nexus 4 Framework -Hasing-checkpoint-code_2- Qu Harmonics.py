def demo_double_compression(seed):
    """
    1) Compress the original seed (a string).
    2) Convert final bits to hex (or text).
    3) Feed that hex (or text) back into the same compression routine.
    4) Compare results.
    """

    print(f"Original Seed: {seed}")

    # Step 1: First compression
    final_bits_1, final_hex_1 = compress_wave_input(seed)  # Suppose it returns (bit_array, hex_string)
    print(f"First Compression - Final Bits (length): {len(final_bits_1)}")
    print(f"First Compression - Final Hex: {final_hex_1[:64]}... (truncated)")

    # Step 2: Convert final bits to a hex or text seed
    # We'll use the final_hex_1 as the "seed" for the next round
    new_seed = final_hex_1  

    print(f"New Seed for Round 2 (Hex): {new_seed[:64]}... (truncated)")

    # Step 3: Second compression
    final_bits_2, final_hex_2 = compress_wave_input(new_seed)
    print(f"Second Compression - Final Bits (length): {len(final_bits_2)}")
    print(f"Second Compression - Final Hex: {final_hex_2[:64]}... (truncated)")

    # (Optional) If you want to see if you can decompress from final_bits_2:
    # decompressed_2 = some_decompress_method(final_bits_2) 
    # But that might just yield 'new_seed', not the original text, 
    # because we changed the seed in the second round.

    return (final_bits_1, final_hex_1, final_bits_2, final_hex_2)

# Now run the demo
if __name__ == "__main__":
    seed = "Hello Turbulent World"
    b1, h1, b2, h2 = demo_double_compression(seed)
