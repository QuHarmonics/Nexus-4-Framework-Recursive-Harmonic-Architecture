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

# --- Grid Parameters (the unchanged "grid") ---
generations = 128           # more depth for longer phase space
nonce_limit_per_gen = 128

# === PHASE 1: INNER LOOP ===
# Run the grid until we find our first "good state" (nonce with >=5 leading zeros).
memory_stack = None
good_state = None  # Will store: {'nonce': int, 'leading_zeros': int, 'memory_stack': [...]}
harmonic_nonces = []  # To hold all good states we capture

print("=== Inner Loop Phase: Searching for initial good state (>=5 leading zeros) ===")
for gen in range(generations):
    byte1_nonces, memory_stack = byte1_nonce_generator(limit=nonce_limit_per_gen, memory_stack=memory_stack)
    test_results = try_nonces_with_block_header(byte1_nonces, version, merkle_root, timestamp)
    best_result = test_results[0]
    print(f"[Inner] Generation {gen} | Nonce: {best_result['nonce']} | Zeros: {best_result['leading_zeros']}")
    
    # If we get a nonce with at least 5 leading zeros, capture the good state and seed the grid with it.
    if best_result["leading_zeros"] >= 5:
        good_state = {
            "nonce": best_result["nonce"],
            "leading_zeros": best_result["leading_zeros"],
            "memory_stack": memory_stack[:]  # Preserve the stack exactly at this moment
        }
        harmonic_nonces.append(good_state)
        # Reset the memory stack using the nonce (as in your original logic)
        memory_stack = [best_result["nonce"], best_result["nonce"]]
        print(f"*** Initial good state found: {good_state} ***")
        break
    else:
        # Otherwise, update the memory stack as before
        seed_a = len(bin(best_result["nonce"])[2:])
        seed_b = best_result["leading_zeros"] + 1 or 4
        memory_stack = [seed_a, seed_b] + memory_stack

if good_state is None:
    print("No good state (nonce with >=5 zeros) found in the inner loop.")
    exit(1)

# === PHASE 2: OUTER LOOP ===
# Using the saved good state's memory stack, run additional iterations to refine the result.
outer_iterations = 128
print("\n=== Outer Loop Phase: Refining good state ===")
for gen in range(outer_iterations):
    # Start each outer iteration by seeding with the preserved good state's memory stack.
    memory_stack = good_state["memory_stack"][:]
    byte1_nonces, memory_stack = byte1_nonce_generator(limit=nonce_limit_per_gen, memory_stack=memory_stack)
    test_results = try_nonces_with_block_header(byte1_nonces, version, merkle_root, timestamp)
    best_result = test_results[0]
    print(f"[Outer] Generation {gen} | Nonce: {best_result['nonce']} | Zeros: {best_result['leading_zeros']}")
    
    # If the new result is improved—that is, it still has >=5 zeros and a lower nonce value—
    # then update the good state.
    if best_result["leading_zeros"] >= 5 and best_result["nonce"] < good_state["nonce"]:
        good_state = {
            "nonce": best_result["nonce"],
            "leading_zeros": best_result["leading_zeros"],
            "memory_stack": memory_stack[:]  # save this new good state's memory
        }
        harmonic_nonces.append(good_state)
        print(f"*** Updated good state found: {good_state} ***")
    else:
        # Otherwise, if the new nonce is worse (i.e. a higher value) then revert to the good state's memory stack.
        memory_stack = good_state["memory_stack"][:]
        print("No improvement; reverting to good state's memory stack.")

print("\n=== Final Harmonic Nonces ===")
for h in harmonic_nonces:
    print(h)
