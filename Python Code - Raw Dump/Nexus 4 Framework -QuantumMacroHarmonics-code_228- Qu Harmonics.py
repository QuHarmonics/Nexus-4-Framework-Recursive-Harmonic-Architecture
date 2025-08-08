def unfold_hash_based_on_center(hash_bits, max_length=512):
    def hex_to_bin(hex_string):
        return bin(int(hex_string, 16))[2:].zfill(256)

    def calculate_center_point(array_length):
        # Calculate the center point; if even, it's fractional; if odd, it's whole
        return array_length / 2

    def insert_at_position(hash_array, position, value):
        # Insert value at the specified position
        return hash_array[:position] + [value] + hash_array[position:]

    # Initialize hash array
    hash_array = list(hex_to_bin(hash_bits))
    flip_bit = 0

    while len(hash_array) < max_length:
        # Calculate center point
        center_point = calculate_center_point(len(hash_array))
        is_fractional = center_point != int(center_point)
        center_index = int(center_point)

        if is_fractional:
            # Expand by adding new bits
            new_value = '0'  # Example value; adjust as needed
            if flip_bit == 0:
                # Insert to the left of the center
                hash_array = insert_at_position(hash_array, center_index, new_value)
            else:
                # Insert to the right of the center
                hash_array = insert_at_position(hash_array, center_index + 1, new_value)
        else:
            # Apply XOR operation to change state
            if flip_bit == 0:
                # XOR with the left neighbor
                if center_index > 0:
                    hash_array[center_index] = str(int(hash_array[center_index]) ^ int(hash_array[center_index - 1]))
            else:
                # XOR with the right neighbor
                if center_index < len(hash_array) - 1:
                    hash_array[center_index] = str(int(hash_array[center_index]) ^ int(hash_array[center_index + 1]))

        # Update flip_bit
        flip_bit ^= 1

    return ''.join(hash_array[:max_length])

# Example usage
hash_bits = "8f434346648f6b96df89dda901c5176b10a6d83961dd3c1ac88b59b2dc327aa4"
unfolded_hash = unfold_hash_based_on_center(hash_bits)
print("Unfolded Hash:", unfolded_hash)
print("Length:", len(unfolded_hash))
