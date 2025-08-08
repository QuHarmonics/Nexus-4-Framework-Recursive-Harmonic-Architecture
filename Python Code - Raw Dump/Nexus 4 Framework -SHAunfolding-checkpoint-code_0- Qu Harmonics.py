def pull_hash_inside_out(hash_bits, debug=False):
    def base_shift(value, current_base, next_base):
        # Convert a number from one base to another
        digits = []
        while value:
            digits.append(value % next_base)
            value //= next_base
        return digits[::-1]  # Preserve order

    # Start with the hash in base 16
    value_base_16 = int(hash_bits, 16)

    # Transition to Base 2 (expanded form)
    expanded_digits = base_shift(value_base_16, 16, 2)
    expanded_value = ''.join(map(str, expanded_digits))

    # Simultaneously compress back to Base 16 (compressed form)
    compressed_digits = base_shift(int(expanded_value, 2), 2, 16)
    compressed_value = ''.join(map(lambda x: hex(x)[2:].upper(), compressed_digits))

    if debug:
        print(f"Base 16 (Compressed): {hash_bits}")
        print(f"Base 2 (Expanded): {expanded_value} (Length: {len(expanded_value)})")
        print(f"Base 16 (Recompressed): {compressed_value}")

    return expanded_value, compressed_value

# Test the inside-out transition
hash_bits = "8f434346648f6b96df89dda901c5176b10a6d83961dd3c1ac88b59b2dc327aa4"
expanded, recompressed = pull_hash_inside_out(hash_bits, debug=True)
