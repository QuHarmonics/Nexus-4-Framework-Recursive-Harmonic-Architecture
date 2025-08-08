def samson_mark1_sha_analysis(hash_bits, initial_base=16, target_base=2, debug=False):
    """
    Analyze and expand a SHA hash using Samson and Mark1 principles.
    Two-way base transitions: 16 → 2 and 2 → 16, then analyze the middle state.
    """

    def base_shift(value, current_base, next_base):
        # Convert a number from one base to another
        digits = []
        while value:
            digits.append(value % next_base)
            value //= next_base
        return digits[::-1]  # Preserve order

    def samson_align_harmonics(digits, harmonic_factor=1):
        # Align digits harmonically (e.g., adjust for resonance)
        return [(digit + harmonic_factor) % 2 for digit in digits]

    # Correctly process the hash bits as hexadecimal (Base 16)
    value_base_16 = int(hash_bits, initial_base)  # Start with hash in Base 16

    # Convert Base 16 to Binary (Base 2)
    binary_representation = bin(value_base_16)[2:].zfill(256)  # Ensure it is 256 bits
    value_base_2 = int(binary_representation, 2)

    # Transition Base 16 → Base 2
    expanded_digits_16_to_2 = base_shift(value_base_16, 16, 2)
    harmonized_16_to_2 = samson_align_harmonics(expanded_digits_16_to_2)

    # Transition Base 2 → Base 16
    compressed_digits_2_to_16 = base_shift(value_base_2, 2, 16)
    harmonized_2_to_16 = samson_align_harmonics(compressed_digits_2_to_16)

    # Middle Ground: Base 10 Average
    average_base_10 = (int(''.join(map(str, harmonized_16_to_2)), 2) +
                       int(''.join(map(str, harmonized_2_to_16)), 16)) // 2

    if debug:
        print(f"Base 16 → Base 2: {harmonized_16_to_2[:128]}... (Length: {len(harmonized_16_to_2)})")
        print(f"Base 2 → Base 16: {harmonized_2_to_16[:64]}... (Length: {len(harmonized_2_to_16)})")
        print(f"Middle Base 10: {average_base_10}")

    return harmonized_16_to_2, harmonized_2_to_16, average_base_10

# Test the guided expansion
hash_bits = "8f434346648f6b96df89dda901c5176b10a6d83961dd3c1ac88b59b2dc327aa4"
harm_16_to_2, harm_2_to_16, avg_base_10 = samson_mark1_sha_analysis(hash_bits, debug=True)

print("Harmonized Base 16 → Base 2:")
print(harm_16_to_2)

print("\nHarmonized Base 2 → Base 16:")
print(harm_2_to_16)

print("\nAverage Base 10 Representation:")
print(avg_base_10)
