def hex_to_bin(hex_string):
    return bin(int(hex_string, 16))[2:].zfill(256)


def unfold_hash(hash_bits):
    hash_bits = hex_to_bin(hash_bits)
    unfolded_hash = list(hash_bits)
    x_pointer = len(hash_bits) // 2
    z_pointer = 1
    flip_bit = 0

    while len(unfolded_hash) < 516:
        # Expanding circle
        z_pointer **= 2

        # Calculate reflected index
        reflected_index = len(unfolded_hash) - x_pointer + 1

        # Ensure indices are within valid range
        if reflected_index - 1 < 0:
            reflected_index = 1
        if reflected_index + 1 >= len(unfolded_hash):
            reflected_index = len(unfolded_hash) - 2
        if x_pointer >= len(unfolded_hash):
            x_pointer = len(unfolded_hash) - 1

        # Flip bit logic
        if flip_bit == 0:
            # High: XOR operation on reflected location
            xor_result = int(unfolded_hash[reflected_index + 1]) ^ int(unfolded_hash[x_pointer])
            unfolded_hash[reflected_index + 1] = str(xor_result)
        else:
            # Low: Add current value to reflected location
            unfolded_hash.insert(reflected_index - 1, unfolded_hash[x_pointer])

        # Update flip bit
        flip_bit ^= 1

        # Update x pointer
        x_pointer += 1

        # Ensure x pointer is within valid range
        if x_pointer >= len(unfolded_hash):
            x_pointer = len(unfolded_hash) - 1

    return ''.join(unfolded_hash[:512])  # Return the first 512 characters


# Example usage
hash_bits = "8f434346648f6b96df89dda901c5176b10a6d83961dd3c1ac88b59b2dc327aa4"
unfolded_hash = unfold_hash(hash_bits)
print("Unfolded Hash:", unfolded_hash)
print("Length:", len(unfolded_hash))