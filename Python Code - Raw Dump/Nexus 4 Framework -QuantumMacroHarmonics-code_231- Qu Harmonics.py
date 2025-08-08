def unfold_hash_with_phases(hash_bits, max_length=256):  # Larger max_length for better testing
    def hex_to_bin(hex_string):
        return bin(int(hex_string, 16))[2:].zfill(256)

    def calculate_center_point(array_length):
        # Calculate the center point; fractional if even length
        return array_length / 2

    def reflect_position(array_length, position):
        # Reflect a position across the center
        center = array_length / 2
        return int(2 * center - position)

    # Initialize hash array
    hash_array = list(hex_to_bin(hash_bits))
    linear_pointer = 1  # Forward progress
    flip_bit = 0

    while len(hash_array) < max_length:
        # Phase 1: Determine center and reflection
        center_point = calculate_center_point(len(hash_array))
        is_fractional = center_point != int(center_point)
        center_index = int(center_point)
        reflected_index = reflect_position(len(hash_array), center_index)

        # Phase 2: Apply operation based on center type and flip_bit
        if is_fractional:
            # Insert a space (inverted value of the current linear pointer)
            inverted_value = str(1 - (linear_pointer % 2))
            if flip_bit == 0:
                # Insert on the left side (before center)
                hash_array.insert(max(reflected_index - 1, 0), inverted_value)
            else:
                # Insert on the right side (after center)
                hash_array.insert(min(reflected_index + 1, len(hash_array)), inverted_value)
        else:
            # Perform XOR at the reflected position
            if 0 <= reflected_index < len(hash_array):
                xor_value = str(int(hash_array[center_index]) ^ int(hash_array[reflected_index]))
                hash_array[reflected_index] = xor_value

        # Phase 3: Update pointers and states
        linear_pointer += 1  # Forward progress
        flip_bit = (flip_bit + 1) % 2  # Alternate between 0 and 1

    return ''.join(hash_array[:max_length])

# Test with adjusted logic
hash_bits = "8f434346648f6b96df89dda901c5176b10a6d83961dd3c1ac88b59b2dc327aa4"
unfolded_hash = unfold_hash_with_phases(hash_bits)
unfolded_hash, len(unfolded_hash)
