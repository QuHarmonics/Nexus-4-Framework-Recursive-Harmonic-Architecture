import numpy as np

def harmonic_feedback_unfolding(hash_bits, target_length=512, alpha=1.5, target=0.5, debug=False):
    def hex_to_bin(hex_string):
        return bin(int(hex_string, 16))[2:].zfill(256)

    def harmonic_feedback(n, previous, target):
        # Recursive harmonic feedback formula
        correction = alpha * (target - previous) / (n + 1)
        return previous * (-0.5) * np.cos(n * np.pi) + correction

    # Initialize unfolding process
    hash_array = list(hex_to_bin(hash_bits))
    linear_pointer = 1  # Forward progress
    harmonic_values = [target]  # Initialize with target for feedback

    while len(hash_array) < target_length:
        array_length = len(hash_array)
        n = len(harmonic_values) - 1
        previous_value = harmonic_values[-1]

        # Compute harmonic feedback
        harmonic_value = harmonic_feedback(n, previous_value, target)
        harmonic_values.append(harmonic_value)

        # Insert feedback value into the center
        center_index = array_length // 2
        hash_array.insert(center_index, str(int(harmonic_value > 0)))

        if debug:
            print(f"Length: {len(hash_array)}, Harmonic Value: {harmonic_value}")

        linear_pointer += 1

    return ''.join(hash_array[:target_length]), harmonic_values

# Test the harmonic feedback unfolding
hash_bits = "8f434346648f6b96df89dda901c5176b10a6d83961dd3c1ac88b59b2dc327aa4"
harmonic_unfolded_hash, harmonic_values = harmonic_feedback_unfolding(hash_bits, debug=True)

harmonic_unfolded_hash[:256], harmonic_unfolded_hash[-256:]
