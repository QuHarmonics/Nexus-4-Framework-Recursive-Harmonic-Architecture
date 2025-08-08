import hashlib

def byte1_nonce_generator(seed_a=0, seed_b=9, limit=1000, memory_stack=None):
    """Generate nonces with recursive folding logic and optional memory stack."""
    stack = memory_stack[:] if memory_stack else [seed_a, seed_b]
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
        nonces.append(stack[-1])
    return nonces, stack

def count_leading_hex_zeros(hex_digest):
    """Count leading '0' characters in a hex digest string (nibble-level)."""
    return len(hex_digest) - len(hex_digest.lstrip('0'))

def try_nonces_with_block_header(nonces, version, merkle_root, timestamp):
    results = []
    for nonce_int in nonces:
        try:
            nonce = nonce_int.to_bytes(4, byteorder='little', signed=False)
        except OverflowError:
            continue
        
        nonce_hex = nonce_int.to_bytes(4, byteorder='little').hex()
        
        block_header = bytearray.fromhex(
            "02000000" +  # Version 2
            "69054f28012b4474caa9e821102655cc74037c415ad2bba70200000000000000" +  # Prev block hash (reversed)
            "2ecfc74ceb512c5055bcff7e57735f7323c32f8bbb48f5e96307e5268c001cc9" +  # Merkle root (reversed)
            "3a09be52" +  # Timestamp (0x52be093a little-endian)
            "0ca30319" +  # Bits (target difficulty in compact format)
             nonce_hex    # Nonce (we start with 0, will iterate)
        )    
        first_hash = hashlib.sha256(block_header).digest()
        second_hash = hashlib.sha256(first_hash).hexdigest()
        leading_zeros = count_leading_hex_zeros(second_hash)
        results.append({
            "nonce": nonce_int,
            "second_hash": second_hash,
            "leading_zeros": leading_zeros
        })
    return sorted(results, key=lambda x: x["leading_zeros"], reverse=True)

# Initial block data


# Recursive harmonic loop
generations = 256
nonce_limit_per_gen = 256
memory_stack = None
best_hashes = []
harmonic_nonces = []

for gen in range(generations):
    byte1_nonces, memory_stack = byte1_nonce_generator(limit=nonce_limit_per_gen, memory_stack=memory_stack)
    test_results = try_nonces_with_block_header(byte1_nonces, version, merkle_root, timestamp)
    if not test_results:
        continue
    best_result = test_results[0]
    best_hashes.append(best_result)

    if best_result["leading_zeros"] >= 4:
        print(f"[!] 🚀 Gen {gen} | Nonce: {best_result['nonce']} | Zeros: {best_result['leading_zeros']} | Hash: {best_result['second_hash']}")
        harmonic_nonces.append(best_result)
        memory_stack = [best_result["nonce"], best_result["nonce"]]
    else:
        seed_a = len(bin(best_result["nonce"])[2:])
        seed_b = best_result["leading_zeros"] + 1 or 4
        memory_stack = [seed_a, seed_b] + memory_stack

# Show harmonic zero matches
print("\n🔮 Harmonic Nonces (>= 5 hex-zeros):")
for h in harmonic_nonces:
    if h["leading_zeros"] >= 5:
        print(h)
