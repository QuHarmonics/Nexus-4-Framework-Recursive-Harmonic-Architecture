def unfold_hash_with_alternating_expansion(hash_bits, max_length=512):
    def hex_to_bin(hex_string):
        return bin(int(hex_string, 16))[2:].zfill(256)

    def calculate_midpoint(array_length):
        # Calculate the dynamic midpoint
        return array_length // 2

    def insert_alternating(hash_array, pointer_value, flip_bit):
        # Calculate midpoint
        midpoint = calculate_midpoint(len(hash_array))

        if flip_bit == 1:
            # Insert on the right side of the midpoint
            return hash_array[:midpoint + 1] + [pointer_value] + hash_array[midpoint + 1:]
        else:
            # Insert on the left side of the midpoint
            return hash_array[:midpoint] + [pointer_value] + hash_array[midpoint:]

    # Initialize hash array
    hash_array = list(hex_to_bin(hash_bits))
    current_pointer = 1  # Start the pointer at 1
    flip_bit = 0

    while len(hash_array) < max_length:
        # Calculate pointer value (as a binary string)
        pointer_value = str(current_pointer % 2)

        # Insert new bit alternately based on the flip bit
        hash_array = insert_alternating(hash_array, pointer_value, flip_bit)

        # Update pointers
        current_pointer += 1
        flip_bit ^= 1

    return ''.join(hash_array[:max_length])

# Example usage
hash_bits = "8f434346648f6b96df89dda901c5176b10a6d83961dd3c1ac88b59b2dc327aa4"
unfolded_hash = unfold_hash_with_alternating_expansion(hash_bits)
print("Unfolded Hash:", unfolded_hash)
print("Length:", len(unfolded_hash))
