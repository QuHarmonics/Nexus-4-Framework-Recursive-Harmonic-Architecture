import hashlib

# === Config ===
ITERATIONS = 256
SEED = "1"

# === Utilities ===
def sha256_hex(data: str) -> str:
    return hashlib.sha256(data.encode()).hexdigest()

def to_hex_from_text(text: str) -> str:
    return ''.join(f"{ord(c):02x}" for c in text)

def trailing_zeros(hex_str: str) -> int:
    binary = bin(int(hex_str, 16))[2:].zfill(256)
    return len(binary) - len(binary.rstrip('0'))

# === SHA Reflective Trickster Loop ===
def sha_reflective_series(seed: str, steps: int):
    hashes = []
    current = seed
    for i in range(steps):
        hex_input = to_hex_from_text(current)  # convert text to hex as text representation
        digest = sha256_hex(hex_input)
        tzeros = trailing_zeros(digest)
        hashes.append((i + 1, hex_input, digest, tzeros))
        current = digest
    return hashes

# Generate the full 256-step series
series = sha_reflective_series(SEED, ITERATIONS)

# Display full output, no truncation
import pandas as pd
import ace_tools as tools

df = pd.DataFrame(series, columns=["Step", "Input (Text→Hex)", "SHA-256 Digest", "Trailing Zeros (Binary)"])
tools.display_dataframe_to_user(name="SHA Reflective Harmonic Series", dataframe=df)

