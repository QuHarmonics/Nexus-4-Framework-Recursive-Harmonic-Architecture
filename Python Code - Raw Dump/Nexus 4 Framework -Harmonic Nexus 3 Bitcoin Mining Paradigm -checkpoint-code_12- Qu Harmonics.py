import hashlib

def byte1_nonce_generator(seed_a=3, seed_b=5, limit=10000, memory_stack=None):
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

# --- Initial block data ---
version = bytes.fromhex('20000000')
merkle_root = bytes.fromhex('4d5e6f7a8b9c0d1e2f3a4b5c6d7e8f90123456789abcdef0fedcba9876543210')
timestamp = bytes.fromhex('5f5e1000')

# --- Grid parameters (unchanged) ---
generations = 128           # more depth for longer phase space
nonce_limit_per_gen = 128

# === PHASE 1: RUN THE GRID (INNER LOOP) ===
# Start with no memory stack. Run until we get our first good nonce (>= 5 leading zeros).
memory_stack = None
harmonic_nonces = []  # to hold good nonce(s)
good_found = False

print("=== Running Grid Phase: Searching for initial good state (>=5 leading zeros) ===")
for gen in range(generations):
    byte1_nonces, memory_stack = byte1_nonce_generator(limit=nonce_limit_per_gen, memory_stack=memory_stack)
    test_results = try_nonces_with_block_header(byte1_nonces, version, merkle_root, timestamp)
    best_result = test_results[0]
    print(f"[Grid] Generation {gen} | Nonce: {best_result['nonce']} | Zeros: {best_result['leading_zeros']}")
    
    if best_result["leading_zeros"] >= 5:
        harmonic_nonces.append(best_result)
        good_found = True
        print(f"*** Good state found at generation {gen}: {best_result} ***")
        # Break out of the grid loop once a good nonce is found.
        break
    else:
        # Otherwise update the memory stack as in your original code.
        seed_a = len(bin(best_result["nonce"])[2:])
        seed_b = best_result["leading_zeros"] + 1 or 4
        memory_stack = [seed_a, seed_b] + memory_stack

if not good_found:
    print("No nonce with >= 5 leading zeros was found in the grid phase.")
    exit(1)

# === PHASE 2: FINE-TUNING AROUND THE GOOD NONCE ===
# Use the nonce found above as the starting value. For each candidate,
# test nonce+1 and nonce-1 for 10,000 iterations.

base_nonce = harmonic_nonces[-1]["nonce"]
fine_tune_iterations = 10000
fine_results = []

print("\n=== Fine-Tuning Phase ===")
print(f"Using base nonce: {base_nonce}")

for offset in range(1, fine_tune_iterations + 1):
    for candidate in [base_nonce + offset, base_nonce - offset]:
        nonce = candidate.to_bytes(4, byteorder='little', signed=False)
        block_header = version + merkle_root + timestamp + nonce
        first_hash = hashlib.sha256(block_header).digest()
        second_hash = hashlib.sha256(first_hash).hexdigest()
        result = {
            "nonce": candidate,
            "second_hash": second_hash,
            "leading_zeros": len(second_hash) - len(second_hash.lstrip('0'))
        }
        fine_results.append(result)

# Sort fine-tuning results by number of leading zeros (highest first)
fine_results = sorted(fine_results, key=lambda x: x["leading_zeros"], reverse=True)

print("\n=== Top 10 Fine-Tuning Results ===")
for res in fine_results[:10]:
    print(res)
