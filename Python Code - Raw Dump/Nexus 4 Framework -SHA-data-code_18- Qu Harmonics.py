import hashlib

# Toy target: very low difficulty for demonstration
TARGET = 0x00000FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF

def double_sha256(header_bytes: bytes) -> bytes:
    return hashlib.sha256(hashlib.sha256(header_bytes).digest()).digest()

# Fixed header fields (except nonce)
version = (1).to_bytes(4, 'little')
prev_hash = bytes.fromhex('00' * 32)
merkle_root = bytes.fromhex('3ba3edfd7a7b12b27ac72c3e67768f617fc81bc3888a51323a9fb8aa4b1e5e4a')
timestamp = int(1609459200).to_bytes(4, 'little')  # e.g., 2021-01-01 UTC
bits = bytes.fromhex('1d00ffff')

for nonce in range(0, 2**32):
    header = version + prev_hash + merkle_root + timestamp + bits + nonce.to_bytes(4, 'little')
    h = int.from_bytes(double_sha256(header), 'little')
    if h <= TARGET:
        print(f"Found valid nonce: {nonce}, hash: {h:064x}")
        break