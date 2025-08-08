def unfold_hash(hash_bits):
    unfolded_hash = list(hash_bits)
    x_pointer = len(hash_bits) // 2
    y_pointer = 1
    z_pointer = 1

    for _ in range(len(hash_bits)):
        # X-pointer movement
        if z_pointer % 2 == 0:
            x_pointer += 0.5
        else:
            x_pointer -= 0.5

        # Y-pointer movement
        y_pointer += 1

        # Z-pointer movement
        z_pointer **= 2

        # Calculate insert index
        insert_index = int(x_pointer)

        # Check if insert index is within bounds
        if insert_index < len(unfolded_hash):
            # Perform XOR operation
            if z_pointer % 2 == 0:
                unfolded_hash[insert_index] = str(int(unfolded_hash[insert_index]) ^ int(unfolded_hash[(insert_index - 1) % len(unfolded_hash)]))
            # Insert space
            else:
                unfolded_hash.insert(insert_index, '0')
        else:
            # Append space
            unfolded_hash.append('0')

    return unfolded_hash


# Example usage
hash_bits = "8f434346648f6b96df89dda901c5176b10a6d83961dd3c1ac88b59b2dc327aa4"  # Replace with your 256-bit hash
unfolded_hash = unfold_hash(hash_bits)
print("Unfolded Hash:", unfolded_hash)
print("Length:", len(unfolded_hash))