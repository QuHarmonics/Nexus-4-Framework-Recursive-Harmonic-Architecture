import hashlib

# === NONCE GENERATOR USING RECURSIVE FOLDING SYSTEM ===
def byte1_nonce_generator(seed_a=3, seed_b=5, limit=1, memory_stack=None):
    """
    Generates a list of nonce values using a recursive self-modifying stack.
    This method simulates harmonic folding and evolution over iterations.
    
    Parameters:
        seed_a (int): Initial seed A value.
        seed_b (int): Initial seed B value.
        limit (int): Number of nonce generations to perform.
        memory_stack (list): Optional input memory to continue a prior evolution.

    Returns:
        (list, list): Tuple of nonce list and evolved memory stack.
    """
    stack = memory_stack[:] if memory_stack else [seed_a, seed_b]  # Copy memory stack if it exists, else start with seeds
    nonces = []  # Stores resulting nonce values

    for _ in range(limit):
        if len(stack) < 2:
            break  # Stack too small to operate on

        A = stack[-2]  # Second-to-last element
        B = stack[-1]  # Last element

        # Calculate bit length of (B - A), or use 1 if difference is non-positive
        C = len(bin(B - A)[2:]) if B - A > 0 else 1
        stack.append(C)  # Folded depth measure
        stack.append(C)  # Duplicate to reinforce structure

        Z = A + B  # Recursive additive sum
        stack.append(Z)  # Extend structure

        D = Z - B  # Isolate contribution from A
        stack[-3] = D  # Replace earlier value for harmonic phase adjustment

        Y = Z + B  # Future growth projection
        stack.append(Y)

        X = len(stack)  # Encode current stack size
        stack.append(X)

        total = sum(stack)  # Total recursive energy
        compressed = len(bin(total)[2:])  # Encoded compression phase
        stack.append(compressed)

        close_byte = A + B  # Echo of closure
        stack.append(close_byte)

        nonces.append(stack[-1])  # Capture most recent value as nonce

    return nonces, stack  # Output all generated nonces and final memory stack


# === TEST GENERATED NONCES AGAINST BLOCK HEADER HASHING ===
def try_nonces_with_block_header(nonces, version, merkle_root, timestamp):
    """
    Tests a list of nonce integers against a simulated block header.
    Double-hashes each composed header and checks leading zero count.

    Parameters:
        nonces (list): List of nonce integers.
        version (bytes): Block version field.
        merkle_root (bytes): Root of the Merkle tree.
        timestamp (bytes): Block timestamp.

    Returns:
        list: Sorted list of dicts with hash and zero-leading analysis.
    """
    results = []
    for nonce_int in nonces:
        nonce = nonce_int.to_bytes(4, byteorder='little', signed=False)  # Convert integer to 4-byte little endian
        block_header = version + merkle_root + timestamp + nonce  # Concatenate header components
        first_hash = hashlib.sha256(block_header).digest()  # First SHA-256 hash
        second_hash = hashlib.sha256(first_hash).hexdigest()  # Second SHA-256 hash (as hex string)

        results.append({
            "nonce": nonce_int,
            "second_hash": second_hash,
            "leading_zeros": len(second_hash) - len(second_hash.lstrip('0'))  # Count leading zeros in hash
        })
    return sorted(results, key=lambda x: x["leading_zeros"], reverse=True)  # Sort descending by leading zeros


# === BLOCK HEADER COMPONENT DEFINITIONS ===
version = bytes.fromhex('20000000')  # Block version (little-endian hex)
merkle_root = bytes.fromhex('76937d6f5a1349df5408499efe6cbc5d81db07d994f9731a4758b34816b716d9')
timestamp = bytes.fromhex('5f5e1000')  # Example timestamp in hex


# === MAIN RECURSIVE NONCE SEARCH LOOP ===
generations = 128  # Number of evolution cycles to run
nonce_limit_per_gen = 128  # Nonces generated per generation
memory_stack = None  # Initial memory is empty
best_hashes = []  # All generation-best results
harmonic_nonces = []  # Only results exceeding harmonic threshold

for gen in range(generations):
    # Generate nonces and evolve memory
    byte1_nonces, memory_stack = byte1_nonce_generator(limit=nonce_limit_per_gen, memory_stack=memory_stack)

    # Test those nonces against block hash format
    test_results = try_nonces_with_block_header(byte1_nonces, version, merkle_root, timestamp)
    best_result = test_results[0]  # Best result in this generation
    best_hashes.append(best_result)  # Track best overall results

    if best_result["leading_zeros"] >= 3:
        # Trigger harmonic reset: reset memory using best nonce
        harmonic_nonces.append(best_result)
        nonce_val = best_result["nonce"]
        memory_stack = [nonce_val, nonce_val]  # Seed next cycle with harmonic resonance
    else:
        # Continue folding logic by updating seeds
        seed_a = len(bin(best_result["nonce"])[2:])  # Bit-length of nonce
        seed_b = best_result["leading_zeros"] + 1  # Reinforce based on zero-count
        memory_stack = [seed_a, seed_b] + memory_stack  # Prepend new seeds to evolve stack

# Preview first 5 harmonic matches (can be removed or modified as needed)
harmonic_nonces[:5]
