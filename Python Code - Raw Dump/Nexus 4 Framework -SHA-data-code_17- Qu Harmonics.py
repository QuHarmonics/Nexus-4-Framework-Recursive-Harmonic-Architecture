import hashlib

def sha256_hash(state):
    message = bytes([state % 256])
    return hashlib.sha256(message).digest()

def normalize_hash(hash_bytes):
    hash_int = int.from_bytes(hash_bytes, 'big')
    return hash_int / (2 ** 256)

def delta_phase(state, hash_bytes):
    a = state % 256
    b = int.from_bytes(hash_bytes, 'big') % 256
    return b - a

def next_state(state, delta):
    return (state + delta) % 256

state = 0
for t in range(5):
    hash_bytes = sha256_hash(state)
    norm_hash = normalize_hash(hash_bytes)
    delta = delta_phase(state, hash_bytes)
    next_s = next_state(state, delta)
    fold = "Yes" if norm_hash >= 0.35 else "No"
    print(f"Iteration {t}: State = {state}, Norm Hash = {norm_hash:.4f}, Delta = {delta}, Next State = {next_s}, Fold = {fold}")
    state = next_s