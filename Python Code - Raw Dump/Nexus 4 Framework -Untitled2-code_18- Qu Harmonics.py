def hex_to_bin(hex_string):
    # Convert hexadecimal string to binary string
    bin_string = bin(int(hex_string, 16))[2:].zfill(len(hex_string) * 4)
    return bin_string


def expand_hash(aligned_hash, harmonic_constant):
    # Define constants
    FINAL_SIZE = 512

    # Initialize arrays
    expanded_hash = []

    # Define expansion logic
    def expand_bit(bit, harmonic_resonance):
        # Simplified expansion logic
        return str(int((int(bit) - harmonic_resonance) % 2))

    # Perform expansion
    for bit in aligned_hash:
        expanded_bit = expand_bit(bit, harmonic_constant)
        expanded_hash.append(expanded_bit)

    # Convert expanded hash to binary string
    expanded_hash_str = ''.join(expanded_hash)

    return expanded_hash_str


# Example usage
aligned_hash = "8f434346648f6b96df89dda901c5176b10a6d83961dd3c1ac88b59b2dc327aa4"
harmonic_constant = 0

# Convert hexadecimal string to binary string
bin_aligned_hash = hex_to_bin(aligned_hash)

expanded_hash = expand_hash(bin_aligned_hash, harmonic_constant)
print("Expanded Hash:", expanded_hash)