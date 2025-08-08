import hashlib, struct, random

# Parameters for demonstration
difficulty_bits = 16  # target: 16 leading zero bits (for faster finding in demo)
target_zero_count = difficulty_bits
# Example block header components (simplified)
version = 1
prev_hash = b'\x00' * 32  # previous block hash (all zeros for genesis-like scenario)
merkle_root = hashlib.sha256(b'tx1').digest()  # dummy Merkle root of one tx (would normally double SHA, but not critical here)
timestamp = 1637184000  # example timestamp
bits = 0x1f00ffff        # example difficulty bits (not used directly in this demo)
# Construct 76-byte block header (version, prev_hash, merkle_root, time, bits) in little-endian where appropriate
header_no_nonce = struct.pack("<L", version) + prev_hash[::-1] + merkle_root[::-1] + struct.pack("<LL", timestamp, bits)

def double_sha256(data: bytes) -> bytes:
    """Compute double SHA-256 hash of given data."""
    return hashlib.sha256(hashlib.sha256(data).digest()).digest()

def count_leading_zero_bits(h: bytes) -> int:
    """Count number of leading zero bits in a 256-bit hash."""
    count = 0
    for byte in h:
        if byte == 0:
            count += 8
            continue
        # If not a 0 byte, count leading zeros in this byte (0x08 -> 4, 0x0F -> 4, 0x1F -> 3, etc.)
        # Find position of first 1-bit from MSB:
        leading_zeros = 7 - (byte.bit_length() - 1)
        count += leading_zeros
        break  # stop at the first non-zero byte
    return count

def harmonic_mine(prefix: bytes, target_zero_count: int, max_attempts: int = 500000, range_factor: float = 0.35):
    """
    Search for a nonce that yields at least target_zero_count leading zero bits.
    Uses a harmonic alignment strategy: narrows the search range recursively around promising nonces.
    """
    best_nonce = None
    best_zcount = -1
    # Initial search bounds (0 to 2^32-1 for a 32-bit nonce)
    low, high = 0, 2**32 - 1
    current_center = random.randrange(low, high+1)  # start at a random nonce
    current_range = high - low
    attempts = 0
    found = False
    
    while attempts < max_attempts and not found:
        # Select a candidate nonce around the current center within the current range
        half_range = int(current_range * range_factor / 2)
        start = max(low, current_center - half_range)
        end   = min(high, current_center + half_range)
        if start > end:
            # If range has collapsed too much, reset to full range (or break)
            start, end = low, high
        nonce = random.randrange(start, end+1)
        
        # Compute the double SHA-256 hash for prefix+nonce
        full_header = prefix + struct.pack("<L", nonce)
        hash_val = double_sha256(full_header)
        zcount = count_leading_zero_bits(hash_val)
        attempts += 1
        
        # Check if this attempt is our best so far or meets the target
        if zcount > best_zcount:
            best_zcount = zcount
            best_nonce = nonce
            # Re-center search around this promising nonce
            current_center = nonce
            # Narrow the search range multiplicatively (steer towards this region)
            current_range = max(1, int(current_range * range_factor))
        # If we meet or exceed target, we can stop
        if zcount >= target_zero_count:
            found = True
    
    # Prepare metrics to return
    result_hash = double_sha256(prefix + struct.pack("<L", best_nonce)) if best_nonce is not None else None
    return {
        "found": found,
        "attempts": attempts,
        "best_leading_zeros": best_zcount,
        "best_nonce": best_nonce,
        "best_hash": result_hash.hex() if result_hash else None
    }

# Run the harmonic mining simulation
result = harmonic_mine(header_no_nonce, target_zero_count=difficulty_bits)
print("Harmonic mining result:", result)