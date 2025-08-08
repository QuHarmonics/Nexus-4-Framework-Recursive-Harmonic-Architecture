def hex_to_bin(hex_string):
    return bin(int(hex_string, 16))[2:].zfill(256)

def unfold_hash_with_harmonics(hash_bits, harmonic_constant=0.35, max_length=512):
    hash_bits = hex_to_bin(hash_bits)
    unfolded_hash = list(hash_bits)
    x_pointer = len(hash_bits) // 2
    flip_bit = 0
    feedback_error = 0

    while len(unfolded_hash) < max_length:
        # Calculate reflection index with harmonic correction
        reflected_index = max(1, len(unfolded_hash) - x_pointer + 1)
        reflected_index = min(reflected_index, len(unfolded_hash) - 2)
        
        # Recursive adjustment using harmonic constant
        feedback_error = (harmonic_constant - len(unfolded_hash) / max_length) * 0.5

        if flip_bit == 0:
            xor_result = (int(unfolded_hash[reflected_index]) ^ int(unfolded_hash[x_pointer])) ^ int(feedback_error > 0)
            unfolded_hash[reflected_index] = str(xor_result)
        else:
            unfolded_hash.insert(reflected_index, unfolded_hash[x_pointer])

        # Update flip bit and pointers
        flip_bit ^= 1
        x_pointer = (x_pointer + 1) % len(unfolded_hash)

    return ''.join(unfolded_hash[:max_length])

# Example usage
hash_bits = "8f434346648f6b96df89dda901c5176b10a6d83961dd3c1ac88b59b2dc327aa4"
unfolded_hash = unfold_hash_with_harmonics(hash_bits)
print("Unfolded Hash:", unfolded_hash)
print("Length:", len(unfolded_hash))
