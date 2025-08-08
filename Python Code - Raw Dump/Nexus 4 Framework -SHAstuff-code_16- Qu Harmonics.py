import hashlib

def sha256_hex(data: str) -> str:
    return hashlib.sha256(data.encode()).hexdigest()

def hex_to_text(hex_str: str) -> str:
    try:
        bytes_obj = bytes.fromhex(hex_str)
        return bytes_obj.decode('latin1')  # latin1 to preserve 1:1 mapping
    except Exception:
        return ''  # fallback if not decodable

def text_to_hex(text: str) -> str:
    return text.encode('latin1').hex()

def trailing_zeros_binary(hex_str: str) -> int:
    binary = bin(int(hex_str, 16))[2:].zfill(256)
    return len(binary) - len(binary.rstrip('0'))

def sha_trickster_series(seed: str, steps: int = 256):
    results = []
    current_input = seed

    for i in range(1, steps + 1):
        digest = sha256_hex(current_input)
        reflected_text = hex_to_text(digest)
        trick_input = text_to_hex(reflected_text)
        tz = trailing_zeros_binary(digest)

        results.append({
            'Step': i,
            'Input': current_input,
            'Digest': digest,
            'Trailing_Zeros': tz
        })

        current_input = trick_input

    return results

# === Run and Display ===
series = sha_trickster_series("1", steps=64)

# Print the result as a table
for row in series:
    print(f"{row['Step']:>3} | {row['Input'][:32]:<32} → {row['Digest']} | T-Zeros: {row['Trailing_Zeros']}")
