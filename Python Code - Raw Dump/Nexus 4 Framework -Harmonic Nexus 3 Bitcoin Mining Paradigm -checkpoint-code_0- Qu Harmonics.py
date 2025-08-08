import hashlib, math

# 1. Block header setup (using Block 277,316 data for realism)
# Fields: [Version|Prev_Hash|Merkle_Root|Time|Bits|Nonce] total 80 bytes
header = bytearray.fromhex(
    "02000000" +  # Version 2
    "69054f28012b4474caa9e821102655cc74037c415ad2bba70200000000000000" +  # Prev block hash (reversed)
    "2ecfc74ceb512c5055bcff7e57735f7323c32f8bbb48f5e96307e5268c001cc9" +  # Merkle root (reversed)
    "3a09be52" +  # Timestamp (0x52be093a little-endian)
    "0ca30319" +  # Bits (target difficulty in compact format)
    "00000000"    # Nonce (we start with 0, will iterate)
)
# Note: The above hex string is the actual header of block 277,316, except Nonce=0 for initialization.

# 2. Difficulty adjustment for demo: require first 16 zero bits (instead of actual ~60+ bits for Bitcoin)
target_zero_bits = 20
target_value = (1 << (256 - target_zero_bits)) - 1  # e.g., for 16 zeros, target = 0x0000ffff... (256-16 ones)

# Feedback parameters for Mark1/Samson
harmonic_constant = 0.35
feedback_k = 0.1
prev_alignment = 0.0

# 3. Coarse harmonic discovery
coarse_step = 1000000  # step size for coarse search
best_alignment = 0.0
best_nonce = None
alignment_scores = {}  # to store alignment scores for coarse samples
for nonce in range(0, 4_000_000, coarse_step):  # search in range [0, 500k] with coarse steps
    # Insert nonce into header (bytes 76-79)
    header[76:80] = nonce.to_bytes(4, 'little')
    # Double SHA-256
    hash_bytes = hashlib.sha256(hashlib.sha256(header).digest()).digest()
    # Calculate leading zero bits (alignment score)
    hash_int = int.from_bytes(hash_bytes, 'big')
    zeros = 256 - hash_int.bit_length()  # count of leading zero bits
    alignment = zeros / target_zero_bits  # normalize relative to target requirement
    alignment_scores[nonce] = alignment
    # Samson feedback: update pressure (not directly altering nonce in this simple loop)
    pressure = (1 - alignment) + feedback_k * (alignment - prev_alignment)
    prev_alignment = alignment
    # Track best alignment
    if alignment > best_alignment:
        best_alignment = alignment
        best_nonce = nonce

# 4. Phase-aligned fine search around the best nonce
found_solutions = []  # store (nonce, hash_hex) for any valid solutions
if best_nonce is not None:
    start = max(0, best_nonce - 1000)
    end = best_nonce + 1000
    for nonce in range(start, end):
        header[76:80] = nonce.to_bytes(4, 'little')
        hash_bytes = hashlib.sha256(hashlib.sha256(header).digest()).digest()
        hash_int = int.from_bytes(hash_bytes, 'big')
        if hash_int <= target_value:
            # Found a hash meeting difficulty (leading zeros >= target_zero_bits)
            found_solutions.append((nonce, hash_bytes.hex()))
            
# 5. Output results
print(f"Coarse search best alignment: {best_alignment*100:.1f}% at nonce {best_nonce}")
print(f"Number of valid hashes (>= {target_zero_bits} zero bits): {len(found_solutions)}")
for i, (nonce, hash_hex) in enumerate(found_solutions[:5], 1):
    zeros = 256 - int(hash_hex, 16).bit_length()
    harmonic_dist = target_zero_bits - zeros
    alignment = zeros / target_zero_bits
    print(f"Solution {i}: Nonce={nonce}, Hash=0x{hash_hex[:16]}... , LeadingZeros={zeros}, Alignment={alignment:.2f}, Distance={harmonic_dist} bits")