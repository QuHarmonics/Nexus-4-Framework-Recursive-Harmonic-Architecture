import numpy as np

def int_to_binary_array(n, bits=256):
    return np.array([int(b) for b in bin(n)[2:].zfill(bits)])

def binary_array_to_int(arr):
    return int(''.join(map(str, arr)), 2)

def rotl(x, n, bits=256):
    return np.roll(x, n)

def sha256_round(hash_state, round_constants):
    hash_state = np.bitwise_xor(hash_state, round_constants)
    hash_state = np.bitwise_and(hash_state, np.roll(hash_state, 2))
    hash_state = rotl(hash_state, 7)
    return hash_state

def simulate_sha256(initial_hash_state, num_rounds):
    hash_state = initial_hash_state
    for i in range(num_rounds):
        round_constants = np.random.randint(0, 2, size=256)
        hash_state = sha256_round(hash_state, round_constants)
    return hash_state

# Initialize random 256-bit hash state
initial_hash_state = np.random.randint(0, 2, size=256)

# Simulate 64 rounds of SHA-256
final_hash_state = simulate_sha256(initial_hash_state, 64)

# Print the final hash state
print(final_hash_state)