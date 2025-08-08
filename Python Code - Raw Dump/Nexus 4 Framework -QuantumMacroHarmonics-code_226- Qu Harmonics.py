def unfold_hash_with_middle_insert(hash_bits, max_length=512):
    def hex_to_bin(hex_string):
        return bin(int(hex_string, 16))[2:].zfill(256)

    def calculate_midpoint(array_length):
        # Calculate the dynamic midpoint
        return array_length // 2

    def insert_at_middle(hash_array, value, flip_bit):
        # Calculate midpoint
        midpoint = calculate_midpoint(len(hash_array))

        if flip_bit == 1:
            # Insert on the right side of the midpoint
            return hash_array[:midpoint + 1] + [value] + hash_array[midpoint + 1:]
        else:
            # Insert on the left side of the midpoint
            return hash_array[:midpoint] + [value] + hash_array[midpoint:]

    # Initialize hash array
    hash_array = list(hex_to_bin(hash_bits))
    flip_bit = 0

    while len(hash_array) < max_length:
        # Calculate current position and value
        midpoint = calculate_midpoint(len(hash_array))
        current_value = hash_array[midpoint]

        # Determine value to insert
        if flip_bit == 0:
            # XOR current value with its reflection
            reflected_index = max(0, len(hash_array) - midpoint - 1)
            new_value = str(int(current_value) ^ int(hash_array[reflected_index]))
        else:
            # Clone the current value
            new_value = current_value

        # Insert new value symmetrically
        hash_array = insert_at_middle(hash_array, new_value, flip_bit)

        # Update flip_bit
        flip_bit ^= 1

    return ''.join(hash_array[:max_length])

# Example usage
hash_bits = "8f434346648f6b96df89dda901c5176b10a6d83961dd3c1ac88b59b2dc327aa4"
unfolded_hash = unfold_hash_with_middle_insert(hash_bits)
print("Unfolded Hash:", unfolded_hash)
print("Length:", len(unfolded_hash))
