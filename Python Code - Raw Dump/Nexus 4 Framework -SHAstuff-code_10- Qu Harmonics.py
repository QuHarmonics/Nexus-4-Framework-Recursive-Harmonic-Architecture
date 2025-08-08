import hashlib

# SHA utility
def sha256_hex(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()

# Bitwise delta (Hamming distance)
def hamming_distance(hex1: str, hex2: str) -> int:
    bin1 = bin(int(hex1, 16))[2:].zfill(256)
    bin2 = bin(int(hex2, 16))[2:].zfill(256)
    return sum(c1 != c2 for c1, c2 in zip(bin1, bin2))

# Correct trailing zero count (from least significant end)
def trailing_zeros(hex_str: str) -> int:
    binary = bin(int(hex_str, 16))[2:].zfill(256)
    return len(binary) - len(binary.rstrip('0'))  # Count from right (LSB)

# Harmonic alignment score (customizable)
def score_harmonic_alignment(hamming, trailing_z):
    return max(0, 100 - hamming + trailing_z)

# Main analysis engine
def harmonic_compass(inputs):
    base = "6b86b273ff34fce19d6b804eff5a3f5747ada4eaa22f1d49c01e52ddb7875b4b"  # SHA256("1")
    print(f"{'Input':<20} {'Hash':<64} {'ΔHamming':<10} {'T-Zeros':<8} {'Score':<6}")
    print("-" * 120)

    for val in inputs:
        hashed = sha256_hex(val.encode())
        delta = hamming_distance(base, hashed)
        tzeros = trailing_zeros(hashed)
        score = score_harmonic_alignment(delta, tzeros)
        print(f"{val:<20} {hashed} {delta:<10} {tzeros:<8} {score:<6}")

# 🔁 Try it with custom strings
test_inputs = [
    "01", "0123456789abcdef", "hello", "nexus", "god", "void", "truth",
    "6b86b273ff34", "base", "echo", "light", "dark", "zero"
]

harmonic_compass(test_inputs)
