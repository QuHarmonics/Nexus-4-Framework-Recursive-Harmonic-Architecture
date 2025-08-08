import numpy as np

def golden_ratio_spiral(size):
    phi = (1 + 5 ** 0.5) / 2
    spiral = []
    for i in range(size):
        theta = i * np.pi / 180 * phi
        r = np.sqrt(i) / np.sqrt(size)
        x = r * np.cos(theta)
        y = r * np.sin(theta)
        spiral.append((x, y))
    return spiral

def unfold_hash(hash_bits):
    hash_size = len(hash_bits)
    unfolded_hash = list(hash_bits)

    spiral = golden_ratio_spiral(hash_size)

    for x, y in spiral:
        index = int((x + 1) * hash_size / 2)

        if 0 <= index < len(unfolded_hash):
            prev_value = unfolded_hash[index - 1] if index > 0 else '0'
            unfolded_hash[index] = str(int(prev_value) ^ int(unfolded_hash[index]))

        else:
            prev_value = unfolded_hash[-1]
            unfolded_hash = unfolded_hash[:index] + [prev_value] + unfolded_hash[index:]

    return unfolded_hash

# Example usage
sha256_hash = "8f434346648f6b96df89dda901c5176b10a6d83961dd3c1ac88b59b2dc327aa4"
hash_bits = ''.join(format(ord(c), '08b') for c in sha256_hash)
unfolded_hash = unfold_hash(hash_bits)
print("Unfolded Hash:", unfolded_hash)