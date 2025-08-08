def unfold_hash_with_pointers(hash_bits, max_length=512):
    def initialize_pointers(array_length):
        # current_position, hash_pointer, exp_pointer, flip_bit, center_pointer
        return array_length // 2 + 1, 1, 1, 0, 0

    def calculate_insert_position(current_position, hash_pointer):
        return current_position - (hash_pointer * hash_pointer)

    def insert_and_shift(hash_array, insert_spot, value):
        # Insert the value at the specified position and shift the array
        return hash_array[:insert_spot] + [value] + hash_array[insert_spot:]

    # Initialize hash array
    hash_array = list(hex_to_bin(hash_bits))
    current_position, hash_pointer, exp_pointer, flip_bit, center_pointer = initialize_pointers(len(hash_array))

    while len(hash_array) < max_length:
        # Calculate insertion spot
        insert_spot = calculate_insert_position(current_position, hash_pointer)

        # Determine value to insert
        if (center_pointer + flip_bit) % 2 == 0:
            # XOR value at reflected position
            value = str(int(hash_array[insert_spot % len(hash_array)]) ^ int(hash_array[current_position]))
        else:
            # Clone value at the current position
            value = hash_array[current_position]

        # Insert and shift array
        hash_array = insert_and_shift(hash_array, max(0, insert_spot), value)

        # Update pointers
        hash_pointer += 1
        exp_pointer *= 2
        current_position = (current_position + 1) % len(hash_array)
        flip_bit ^= 1
        center_pointer -= 2

    return ''.join(hash_array[:max_length])

# Example usage
hash_bits = "8f434346648f6b96df89dda901c5176b10a6d83961dd3c1ac88b59b2dc327aa4"
unfolded_hash = unfold_hash_with_pointers(hash_bits)
print("Unfolded Hash:", unfolded_hash)
print("Length:", len(unfolded_hash))
