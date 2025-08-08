import hashlib
import math

def byte1_nonce_generator(seed_a=1, seed_b=4, limit=10000):
    """Generate nonces based on recursive Byte1 folding logic."""
    stack = [seed_a, seed_b]
    nonces = []
    
    for _ in range(limit):
        if len(stack) < 2:
            break
        A = stack[-2]
        B = stack[-1]
        C = len(bin(B - A)[2:]) if B - A > 0 else 1
        stack.append(C)
        stack.append(C)
        Z = A + B
        stack.append(Z)
        D = Z - B
        stack[-3] = D
        Y = Z + B
        stack.append(Y)
        X = len(stack)
        stack.append(X)
        total = sum(stack)
        compressed = len(bin(total)[2:])
        stack.append(compressed)
        close_byte = A + B
        stack.append(close_byte)
        # Harvest the most recent value as a nonce candidate
        nonces.append(stack[-1])
    return nonces

def try_nonces_with_block_header(nonces, version, merkle_root, timestamp):
    results = []
    for nonce_int in nonces:
        nonce = nonce_int.to_bytes(4, byteorder='little', signed=False)
        block_header = version + merkle_root + timestamp + nonce
        first_hash = hashlib.sha256(block_header).digest()
        second_hash = hashlib.sha256(first_hash).hexdigest()
        results.append({
            "nonce": nonce_int,
            "second_hash": second_hash,
            "leading_zeros": len(second_hash) - len(second_hash.lstrip('0'))
        })
    return sorted(results, key=lambda x: x["leading_zeros"], reverse=True)

# Setup block data
version = bytes.fromhex('20000000')
merkle_root = bytes.fromhex('4d5e6f7a8b9c0d1e2f3a4b5c6d7e8f90123456789abcdef0fedcba9876543210')
timestamp = bytes.fromhex('5f5e1000')

# Generate Byte1 nonces and test them
byte1_nonces = byte1_nonce_generator(limit=20)
test_results = try_nonces_with_block_header(byte1_nonces, version, merkle_root, timestamp)
test_results[:5]  # Show top 5 results with most leading zeros

