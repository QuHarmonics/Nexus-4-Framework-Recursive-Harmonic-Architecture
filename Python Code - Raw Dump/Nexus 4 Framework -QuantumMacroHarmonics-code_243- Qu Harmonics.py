def unfold_hash_harmonized(hash_bits, target_length=512, debug=False):
    def hex_to_bin(hex_string):
        return bin(int(hex_string, 16))[2:].zfill(256)

    def rotate_data(data, steps):
        # Rotate the data by a fixed number of steps
        return data[-steps:] + data[:-steps]

    def modular_addition(a, b):
        # Perform modular addition
        return (a + b) % (2 ** len(bin(a)[2:]))

    # Initialize hash array
    hash_array = list(hex_to_bin(hash_bits))
    linear_pointer = 1  # Forward progress
    wheel_speeds = [1]  # Speeds for each wheel, starting with the center wheel

    while len(hash_array) < target_length:
        array_length = len(hash_array)
        center_point = array_length // 2  # Calculate center point

        if debug:
            print(f"Length: {array_length}, Center: {center_point}, Wheel Speeds: {wheel_speeds}")

        # Rotate previous rings
        for i, speed in enumerate(wheel_speeds):
            start = max(0, center_point - (len(wheel_speeds) - i))
            end = min(len(hash_array), center_point + (len(wheel_speeds) - i))
            if end > start:
                segment = hash_array[start:end]
                hash_array[start:end] = rotate_data(segment, speed)

        # Perform SHA-like transformations
        if 0 < center_point < len(hash_array) - 1:
            left = int(hash_array[center_point - 1])
            center = int(hash_array[center_point])
            right = int(hash_array[center_point + 1])

            # Three-way operation: modular addition + XOR + AND
            transformed_value = modular_addition(left, right) ^ (center & right)
            hash_array.insert(center_point, str(transformed_value % 2))  # Ensure bit-level insertion

        # Increment wheel speeds and pointer
        if len(wheel_speeds) < len(hash_array) // 2:
            wheel_speeds.append(len(wheel_speeds) + 1)
        linear_pointer += 1

    return ''.join(hash_array[:target_length]), len(hash_array)

# Test the updated harmonized logic
hash_bits = "8f434346648f6b96df89dda901c5176b10a6d83961dd3c1ac88b59b2dc327aa4"
harmonized_hash, harmonized_length = unfold_hash_harmonized(hash_bits, target_length=512, debug=True)
harmonized_hash, harmonized_length
