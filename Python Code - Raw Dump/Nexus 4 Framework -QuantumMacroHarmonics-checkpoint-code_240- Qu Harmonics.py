def unfold_hash_with_xor_and_swing(hash_bits, target_length=512, debug=False):
    def hex_to_bin(hex_string):
        return bin(int(hex_string, 16))[2:].zfill(256)

    # Initialize hash array
    hash_array = list(hex_to_bin(hash_bits))
    linear_pointer = 1  # Forward progress
    flip_bit = 0

    while len(hash_array) < target_length:
        array_length = len(hash_array)
        center_point = array_length / 2  # Calculate center point
        is_fractional = center_point != int(center_point)
        center_index = int(center_point)

        if debug:
            print(f"Length: {array_length}, Center: {center_point}, Fractional: {is_fractional}, Flip Bit: {flip_bit}")

        if is_fractional:
            # Perform XOR operation for the center
            if 0 < center_index < len(hash_array) - 1:
                xor_value = str(int(hash_array[center_index - 1]) ^ int(hash_array[center_index + 1]))
                hash_array.insert(center_index, xor_value)
        else:
            # Insert a space with the XOR of values before and after
            if 1 <= center_index < len(hash_array) - 1:
                xor_value = str(int(hash_array[center_index - 1]) ^ int(hash_array[center_index]))
                if flip_bit == 0:
                    # Insert on the left side
                    hash_array.insert(center_index - 1, xor_value)
                else:
                    # Insert on the right side
                    hash_array.insert(center_index + 1, xor_value)

        # Update the flip bit and the linear pointer
        flip_bit ^= 1
        linear_pointer += 1

    return ''.join(hash_array[:target_length]), len(hash_array)

# Test the refined logic
hash_bits = "8f434346648f6b96df89dda901c5176b10a6d83961dd3c1ac88b59b2dc327aa4"
unfolded_hash, final_length = unfold_hash_with_xor_and_swing(hash_bits, debug=True)
unfolded_hash, final_length
