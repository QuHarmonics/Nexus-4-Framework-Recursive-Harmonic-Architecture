import hashlib, math, random

# Block header components (Genesis block example)
version = bytes.fromhex('01000000')            # Version 1 (little-endian)
prev_hash = bytes.fromhex('00'*32)             # Previous hash (no parent for genesis)
merkle_root = bytes.fromhex('3BA3EDFD7A7B12B27AC72C3E67768F617FC81BC3888A51323A9FB8AA4B1E5E4A')
merkle_root = merkle_root[::-1]                # Convert to internal little-endian
time = bytes.fromhex('29ab5f49')               # Timestamp (0x495FAB29 little-endian)
bits = bytes.fromhex('ffff001d')               # Difficulty bits (0x1d00ffff little-endian)
# (Using actual genesis block values above for context, though difficulty is adjusted in concept)

# Difficulty target for demo: require 16 leading zero *bits* (i.e., hash < 2^(256-16))
TARGET_ZERO_BITS = 20

# Initialize nonce using digits of π (3.141592653 -> 3141592653)
initial_nonce = 3141592653 % (1<<32)
nonce = initial_nonce

# Function to count leading zero bits in a 256-bit hash
def count_leading_zero_bits(h_bytes):
    h_int = int.from_bytes(h_bytes, 'big')
    if h_int == 0:
        return 256  # hash is all zeros (unlikely in practice)
    return 256 - h_int.bit_length()

# Feedback mining loop
best_zero_bits = 0
history_best = []  # track best alignment over time
for iteration in range(1000000):  # cap iterations to avoid infinite loop
    # Construct block header with current nonce
    nonce_bytes = nonce.to_bytes(4, 'little')
    header = version + prev_hash + merkle_root + time + bits + nonce_bytes
    # Double SHA-256
    h1 = hashlib.sha256(header).digest()
    h2 = hashlib.sha256(h1).digest()
    # Check alignment (leading zeros)
    zero_bits = count_leading_zero_bits(h2)
    if zero_bits > best_zero_bits:
        best_zero_bits = zero_bits
    history_best.append(best_zero_bits)
    if zero_bits >= TARGET_ZERO_BITS:
        print(f"Success: Nonce {nonce} yields hash {h2.hex()} with {zero_bits} leading zero bits!")
        break
    # Harmonic feedback: XOR a portion of the hash into the nonce for next try
    feedback = int.from_bytes(h2[:4], 'little')  # take first 32 bits of hash as feedback
    nonce = nonce ^ feedback  # adjust nonce by XORing feedback (phase reflection)
