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
        block_header = bytearray.fromhex(
            "02000000" +  # Version 2
            "69054f28012b4474caa9e821102655cc74037c415ad2bba70200000000000000" +  # Prev block hash (reversed)
            "2ecfc74ceb512c5055bcff7e57735f7323c32f8bbb48f5e96307e5268c001cc9" +  # Merkle root (reversed)
            "3a09be52" +  # Timestamp (0x52be093a little-endian)
            "0ca30319"   # Bits (target difficulty in compact format)
        ) + nonce  # Append nonce bytes here

        first_hash = hashlib.sha256(block_header).digest()
        second_hash = hashlib.sha256(first_hash).hexdigest()
        results.append({
            "nonce": nonce_int,
            "second_hash": second_hash,
            "leading_zeros": len(second_hash) - len(second_hash.lstrip('0'))
        })
    return sorted(results, key=lambda x: x["leading_zeros"], reverse=True)

# Initial block data


# Recursive harmonic loop
generations = 256  # more depth for longer phase space
nonce_limit_per_gen = 256
memory_stack = None
best_hashes = []
harmonic_nonces = []

for gen in range(generations):
    byte1_nonces, memory_stack = byte1_nonce_generator(limit=nonce_limit_per_gen, memory_stack=memory_stack)
    test_results = try_nonces_with_block_header(byte1_nonces, version, merkle_root, timestamp)
    best_result = test_results[0]
    best_hashes.append(best_result)

    if best_result["leading_zeros"] >= 3:
        print(f"[!] 🚀 Generation {gen} | Nonce: {best_result['nonce']} | Zeros: {best_result['leading_zeros']}")
        harmonic_nonces.append(best_result)
        memory_stack = [best_result["nonce"], best_result["nonce"]]
    else:
        seed_a = len(bin(best_result["nonce"])[2:])
        seed_b = best_result["leading_zeros"] + 1 or 4
        memory_stack = [seed_a, seed_b] + memory_stack

print("\n🔮 Harmonic Nonces (5+ zeros):")
for h in harmonic_nonces:
    print(h)

# ------------------------------------------------------------------
# Incremental Fine-Tuning Phase:
#
# 1. Extract the best nonce based on the highest leading_zeros.
# 2. Divide that nonce by 2.
# 3. Loop incrementally from: 
#       start = best_nonce - (best_nonce // 2)
#       to   end   = best_nonce + (best_nonce // 2)
#    and for each candidate nonce, call try_nonces_with_block_header
# ------------------------------------------------------------------

if harmonic_nonces:
    # Extract the best entry based on leading_zeros (if tie, one of them is chosen)
    best_entry = max(harmonic_nonces, key=lambda x: x["leading_zeros"])
    best_nonce = best_entry["nonce"]
    print("\n[+] Best nonce from harmonic_nonces (highest leading zeros):", best_nonce)
    
    # Calculate half the value (using integer division)
    half_range = best_nonce/2
    
    fine_tune_results = []
    
    # Incrementally loop from (best_nonce - half_range) to (best_nonce + half_range)
    for candidate in range(best_nonce - half_range, best_nonce + half_range + 1):
        candidate_result = try_nonces_with_block_header([candidate], version, merkle_root, timestamp)
        # candidate_result is a list with one element.
        fine_tune_results.append(candidate_result[0])
    
    # Sort the fine-tune results by leading zeros (highest first)
    fine_tune_results = sorted(fine_tune_results, key=lambda x: x["leading_zeros"], reverse=True)
    
    print("\n[+] Fine-Tune Results (Top 10):")
    for res in fine_tune_results[:1000]:
        print(res)
else:
    print("No harmonic nonces available for fine-tuning.")
