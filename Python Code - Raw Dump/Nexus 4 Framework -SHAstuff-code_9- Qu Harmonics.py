import hashlib

def sha256_hex(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()

def mirror_hex_string(hex_str: str) -> str:
    return hex_str[::-1]

def reverse_nibbles(hex_str: str) -> str:
    reversed_nibbles = ''
    for char in hex_str:
        bin_nibble = format(int(char, 16), '04b')
        reversed_bin = bin_nibble[::-1]
        reversed_hex = format(int(reversed_bin, 2), 'x')
        reversed_nibbles += reversed_hex
    return reversed_nibbles

def hamming_distance(hex1: str, hex2: str) -> int:
    bin1 = bin(int(hex1, 16))[2:].zfill(256)
    bin2 = bin(int(hex2, 16))[2:].zfill(256)
    return sum(c1 != c2 for c1, c2 in zip(bin1, bin2))

def trailing_zeros(hex_str: str) -> int:
    binary = bin(int(hex_str, 16))[2:].zfill(256)
    return len(binary) - len(binary.rstrip('0'))

def harmonic_mirrorwalker(seed: str, generations: int = 20):
    current = seed.encode()

    print(f"{'Gen':<4} {'Forward Hash':<64} {'T-Zeros':<8} {'ΔMirror':<9} {'ΔNibble':<9}")
    print("-" * 100)

    for i in range(generations):
        forward_hash = sha256_hex(current)
        mirror_input = mirror_hex_string(forward_hash)
        mirror_hash = sha256_hex(mirror_input.encode())
        reversed_nibbles_input = reverse_nibbles(forward_hash)
        nibble_hash = sha256_hex(reversed_nibbles_input.encode())

        # Bitwise analytics
        tzeros = trailing_zeros(forward_hash)
        delta_mirror = hamming_distance(forward_hash, mirror_hash)
        delta_nibble = hamming_distance(forward_hash, nibble_hash)

        # Print row
        print(f"{i+1:<4} {forward_hash:<64} {tzeros:<8} {delta_mirror:<9} {delta_nibble:<9}")

        # Move to next gen
        current = bytes.fromhex(forward_hash)

# 🔑 Run It
seed = "01"  # Minimal symbolic seed
harmonic_mirrorwalker(seed, generations=600)
