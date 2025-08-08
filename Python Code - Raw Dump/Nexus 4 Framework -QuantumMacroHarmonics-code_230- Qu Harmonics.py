def unfold_hash_with_axes(hash_bits, max_length=511):
    def hex_to_bin(hex_string):
        return bin(int(hex_string, 16))[2:].zfill(256)

    def calculate_center_point(array_length):
        # Calculate the center point; fractional if even length
        return array_length / 2

    def reflect_position(array_length, position):
        # Reflect a position across the center
        center = array_length / 2
        return int(2 * center - position)

    def insert_at_position(hash_array, position, value):
        # Insert value at the specified position
        return hash_array[:position] + [value] + hash_array[position:]

    # Initialize hash array
    hash_array = list(hex_to_bin(hash_bits))
    linear_pointer = 1  # Forward progress
    exponential_pointer = 1  # Starts at 1 and grows exponentially
    flip_bit = 0

    while len(hash_array) < max_length:
        # Calculate center point and determine if whole or fractional
        center_point = calculate_center_point(len(hash_array))
        is_fractional = center_point != int(center_point)
        center_index = int(center_point)
        reflected_index = reflect_position(len(hash_array), center_index)

        if is_fractional:
            # Insert a space (inverted value of the current linear pointer)
            inverted_value = str(1 - (linear_pointer % 2))
            if flip_bit == 0:
                # Insert on the left side (reflected position)
                hash_array = insert_at_position(hash_array, reflected_index, inverted_value)
            else:
                # Insert on the right side (center + 1)
                hash_array = insert_at_position(hash_array, reflected_index + 1, inverted_value)
        else:
            # Perform XOR at the reflected position
            if reflected_index < len(hash_array):
                xor_value = str(int(hash_array[center_index]) ^ int(hash_array[reflected_index]))
                hash_array[reflected_index] = xor_value

        # Update pointers
        linear_pointer += 1  # Arrow of time
        exponential_pointer *= 2  # Exponential growth
        flip_bit ^= 1  # Toggle flip bit for alternating sides

    return ''.join(hash_array[:max_length])

# Example usage
hash_bits = "8f434346648f6b96df89dda901c5176b10a6d83961dd3c1ac88b59b2dc327aa4"
unfolded_hash = unfold_hash_with_axes(hash_bits)
print("Unfolded Hash:", unfolded_hash)
print("Length:", len(unfolded_hash))
