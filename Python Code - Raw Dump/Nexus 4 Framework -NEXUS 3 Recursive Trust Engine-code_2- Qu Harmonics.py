import hashlib
import struct
import math

# Universal harmonic constant
C = 0.35

def double_sha(data: bytes) -> bytes:
    """Double-SHA256 “pocket” test at the end."""
    return hashlib.sha256(hashlib.sha256(data).digest()).digest()

def harmonic_error(digest: bytes) -> float:
    """
    Break the 32-byte digest into 64 nibbles, map each to [0..1] via nibble/15,
    subtract C, then take the L2-norm of that error vector.
    """
    err2 = 0.0
    for b in digest:
        hi = (b >> 4) & 0xF
        lo = b & 0xF
        for nib in (hi, lo):
            a = nib / 15.0
            e = a - C
            err2 += e * e
    return math.sqrt(err2)

def find_resonant_nonce(
    msg: bytes,
    epsilon: float = 1e-3,
    max_trials: int = 10_000_000
) -> (int, float, bytes):
    """
    Injects a 32-bit nonce at the end of `msg`, double-hashes it,
    measures its harmonic error, and keeps going until we get under ε
    (or exhaust max_trials).
    Returns (best_nonce, best_error, best_digest).
    """
    best_nonce = 0
    best_error = float('inf')
    best_digest = b''
    for i in range(max_trials):
        # pack as little-endian 32-bit int
        nonce_bytes = struct.pack('<I', i)
        digest = double_sha(msg + nonce_bytes)
        err = harmonic_error(digest)
        if err < best_error:
            best_error = err
            best_nonce = i
            best_digest = digest
            # once we’re “in tune,” stop
            if err < epsilon:
                break
    return best_nonce, best_error, best_digest

if __name__ == "__main__":
    # your message prefix
    M = b"my magic data..."
    nonce, err, digest = find_resonant_nonce(M)
    print(f"Resonant nonce: {nonce}  error: {err:.6f}")
    print("Final double-SHA256:", digest.hex())
