def recursive_xor_decode(last, now):
    # Convert hex strings to byte arrays for XOR operations
    last_bytes = bytes.fromhex(last)
    now_bytes = bytes.fromhex(now)

    # Extract the 20-byte header (first 40 characters in hex)
    header = now_bytes[:20]

    # Work with the rest of the data after the header
    last_remainder = last_bytes[20:]
    now_remainder = now_bytes[20:]

    # Initialize a list to keep track of intermediate results
    results = []

    while now_remainder:
        # XOR the non-matching data
        new_value = bytes(a ^ b for a, b in zip(last_remainder, now_remainder))

        # Store the result
        results.append(new_value.hex())

        # If we find data matching the original size, assume we reached the original Merkle data
        if len(new_value) == 32:  # Standard SHA-256 hash length in bytes
            break

        # Update last and now for the next recursive step
        last_remainder = now_remainder
        now_remainder = new_value

    # Return the reconstructed data
    return {
        "original_merkle": results[-1] if results else None,
        "steps": results,
        "header": header.hex()
    }


# Input strings
last = "9f0a1b2c3d4e5f6a7b8c9d0e1f2a3b4c5d6e7f8a9b0c1d2e3f4a5b6c7d8e9f0a"
now = "9F0A1B2C3D4E5F6A7B8C9D0E1F2A3B4C5D6E0F8F9AFFFFFF9B0C1D2E3F4A5B6C0F8DA8FFFFFF9F0A"

# Perform the recursive XOR decode process
result = recursive_xor_decode(last, now)
print(result)
