import hashlib

# Seed data
pre_sha_input = "33333332333333353334333633333335333333363333333133343331"

# Generate SHA-512 hash
def sha512_hash(input_value):
    return hashlib.sha512(input_value.encode('utf-8')).hexdigest()

# Calculate the hash
generated_hash = sha512_hash(pre_sha_input)

# Compare to the provided hash
provided_hash = "ee1f36c84231330aa9ea6a842a36ff2eb09d5d004c0f449e11df865f309a3ac9d71f115b981eda6d4cfaa531434e4e5f4a73a61720430e8763700e14a1387a05"

print("Generated Hash:", generated_hash)
print("Matches Provided Hash:", generated_hash == provided_hash)
