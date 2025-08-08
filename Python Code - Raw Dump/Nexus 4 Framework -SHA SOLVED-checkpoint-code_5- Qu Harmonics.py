# Define the two sequences (hash and anti-hash) as hexadecimal strings.
hash_hex = '185F8DB32271FE25F561A6FC938B2E264306EC304EDA518007D1764826381969'
anti_hash_hex = '25F561A6FC938B2E264306EC304EDA518007D10F8660'

# Convert to byte sequences.
hash_bytes = bytes.fromhex(hash_hex)
anti_hash_bytes = bytes.fromhex(anti_hash_hex)

# XOR the first 5 bytes of each.
xor_result = bytes(x ^ y for x, y in zip(hash_bytes[:5], anti_hash_bytes[:5]))
print("First 5 bytes XOR:", xor_result.hex())

# Optionally, incorporate the 6th byte of anti_hash as a finishing factor.
finishing_byte = anti_hash_bytes[5]
print("Finishing byte:", hex(finishing_byte))
