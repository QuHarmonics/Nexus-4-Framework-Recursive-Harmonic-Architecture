import hashlib

def byte1_nonce_generator(seed_a=3, seed_b=5, limit=10000, memory_stack=None, direction=1):
    """
    Generate nonces using recursive folding logic and an optional memory stack.
    The `direction` parameter lets us “jump” forward (direction=1, default) or
    backward (direction=-1) by inverting key operations.
    """
    # Start with a copy of the provided memory stack, or initialize it.
    stack = memory_stack[:] if memory_stack is not None else [seed_a, seed_b]
    nonces = []
    
    for _ in range(limit):
        if len(stack) < 2:
            break
        A = stack[-2]
        B = stack[-1]
        # Compute C; invert the sign in reverse mode.
        if direction == 1:
            C = len(bin(B - A)[2:]) if (B - A) > 0 else 1
        else:
            C = - (len(bin(B - A)[2:]) if (B - A) > 0 else 1)
        stack.append(C)
        stack.append(C)
        
        # Compute Z: add in forward, subtract in reverse.
        if direction == 1:
            Z = A + B
        else:
            Z = A - B
        stack.append(Z)
        
        # Adjust the element three from the end.
        if direction == 1:
            D = Z - B
            stack[-3] = D
        else:
            D = Z + B
            stack[-3] = D
        
        # Compute Y.
        if direction == 1:
            Y = Z + B
        else:
            Y = Z - B
        stack.append(Y)
        
        X = len(stack)
        stack.append(X)
        total = sum(stack)
        compressed = len(bin(total)[2:])
        stack.append(compressed)
        
        # Compute a closing byte.
        if direction == 1:
            close_byte = A + B
        else:
            close_byte = A - B
        stack.append(close_byte)
        
        nonces.append(stack[-1])
    return nonces, stack

def try_nonces_with_block_header(nonces, version, merkle_root, timestamp):
    """
    For each nonce, assemble a block header and compute two SHA256 hashes.
    Returns a list of dicts sorted by the number of leading zeros (highest first).
    """
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

# --- Block Header Data ---
version = bytes.fromhex('20000000')
merkle_root = bytes.fromhex('4d5e6f7a8b9c0d1e2f3a4b5c6d7e8f90123456789abcdef0fedcba9876543210')
timestamp = bytes.fromhex('5f5e1000')

# === STEP 2: OUTER LOOP ===
# Now that we have our initial good state, we run the outer loop.
outer_iterations = 10
good_states = [good_state]  # We can store multiple good states if needed.

for outer in range(outer_iterations):
    print(f"\n=== Outer Loop Iteration {outer} ===")
    # Use the latest good state's memory stack as the starting point.
    current_state = good_states[-1]
    memory_stack = current_state["memory_stack"][:]  # Use the preserved stack
    
    # Run a forward iteration.
    nonces, forward_stack = byte1_nonce_generator(limit=nonce_limit_per_gen, memory_stack=memory_stack, direction=1)
    forward_results = try_nonces_with_block_header(nonces, version, merkle_root, timestamp)
    forward_best = forward_results[0]
    print(f"Forward: Nonce: {forward_best['nonce']} | Zeros: {forward_best['leading_zeros']}")
    
    # Check: if the forward nonce is worse (i.e. a higher nonce value) than our current good state's nonce,
    # we initiate a backward jump.
    if forward_best["nonce"] > current_state["nonce"]:
        print("Forward nonce is higher than current good state. Initiating backward jump.")
        nonces_b, backward_stack = byte1_nonce_generator(limit=nonce_limit_per_gen,
                                                          memory_stack=current_state["memory_stack"],
                                                          direction=-1)
        backward_results = try_nonces_with_block_header(nonces_b, version, merkle_root, timestamp)
        backward_best = backward_results[0]
        print(f"Backward: Nonce: {backward_best['nonce']} | Zeros: {backward_best['leading_zeros']}")
        # If the backward jump yields an improved (lower) nonce with >=5 zeros, update the good state.
        if backward_best["leading_zeros"] >= 5 and backward_best["nonce"] < current_state["nonce"]:
            new_state = {
                "nonce": backward_best["nonce"],
                "leading_zeros": backward_best["leading_zeros"],
                "memory_stack": backward_stack[:]
            }
            good_states.append(new_state)
            print("Updated good state from backward jump:", new_state)
        else:
            # No improvement; keep the current state.
            print("No improvement from backward jump; retaining current good state.")
    else:
        # If forward yields an improved state (lower nonce) and has >=5 zeros, update.
        if forward_best["leading_zeros"] >= 5 and forward_best["nonce"] < current_state["nonce"]:
            new_state = {
                "nonce": forward_best["nonce"],
                "leading_zeros": forward_best["leading_zeros"],
                "memory_stack": forward_stack[:]
            }
            good_states.append(new_state)
            print("Updated good state from forward iteration:", new_state)
        else:
            print("Forward iteration did not improve; maintaining current good state.")

# === FINAL REPORT ===
print("\n=== Final Good States ===")
for state in good_states:
    print(state)
