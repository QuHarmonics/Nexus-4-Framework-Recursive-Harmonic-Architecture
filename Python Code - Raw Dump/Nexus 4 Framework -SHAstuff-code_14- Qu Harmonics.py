import hashlib

# === Config ===
ITERATIONS = 256
SEED = "1"

# === Utilities ===

def sha256_hex(text: str) -> str:
    """SHA-256 of a string."""
    return hashlib.sha256(text.encode()).hexdigest()

def text_to_ascii_hex(text: str) -> str:
    """Convert each character to its ASCII hex representation."""
    return ''.join(format(ord(c), '02x') for c in text)

def trailing_zeros(hex_str: str) -> int:
    """Count number of trailing zeros in the binary form of a SHA hash."""
    binary = bin(int(hex_str, 16))[2:].zfill(256)
    return len(binary) - len(binary.rstrip('0'))

def hamming_distance(hex1: str, hex2: str) -> int:
    """Hamming distance between two SHA256 hex strings."""
    bin1 = bin(int(hex1, 16))[2:].zfill(256)
    bin2 = bin(int(hex2, 16))[2:].zfill(256)
    return sum(c1 != c2 for c1, c2 in zip(bin1, bin2))

# === SHA Harmonic Series: ASCII-Hex Reflection ===

def recursive_sha_ascii_hex(seed: str, steps: int):
    results = []
    current = seed
    for i in range(steps):
        digest = sha256_hex(current)
        zeros = trailing_zeros(digest)
        results.append((i+1, current, digest, zeros))
        current = text_to_ascii_hex(digest)  # Trick SHA: treat digest as raw text and convert to hex
    return results

# === Display ===

def display_series(series):
    print(f"{'Step':>4} | {'Input (truncated)':<32} → {'Digest'} | T-Zeros")
    print("-" * 96)
    for step, seed, digest, zeros in series:
        print(f"{step:>4} | {seed[:32]:<32} → {digest} | {zeros}")

# === Run ===

if __name__ == "__main__":
    sha_series = recursive_sha_ascii_hex(SEED, ITERATIONS)
    display_series(sha_series)
