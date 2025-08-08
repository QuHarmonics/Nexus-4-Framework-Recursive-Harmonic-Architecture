import hashlib
import random
import math

# --- CONFIGURATION ---
TARGET_PREFIX = "00000"  # real PoW hash condition (adjust difficulty)
MAX_ITER = 1000000       # how many iterations to try
K_FEEDBACK = 0.35         # Samson harmonic correction constant
OMEGA = 2 * math.pi / 32  # Golden phase

# Sample realistic header template
def make_header(previous_hash, merkle_root, nonce):
    """
    Builds a block header by concatenating parts. Use only real, known values.
    This mimics a real Bitcoin-like block header structure.
    """
    version = "00000001"
    time = "5f5e1000"        # Unix time
    bits = "1d00ffff"        # Difficulty target
    header_hex = (
        version +
        previous_hash +
        merkle_root +
        time +
        bits +
        format(nonce, '08x')
    )
    return bytes.fromhex(header_hex)

def sha256d(data):
    return hashlib.sha256(hashlib.sha256(data).digest()).hexdigest()

def recursive_sha_unfold(previous_hash, merkle_root, target_prefix=TARGET_PREFIX, max_iter=MAX_ITER):
    nonce = random.getrandbits(32)
    history = []

    for n in range(max_iter):
        header = make_header(previous_hash, merkle_root, nonce)
        h = sha256d(header)

        if h.startswith(target_prefix):
            print("FOUND:", h, "Nonce:", nonce)
            return nonce, h, history

        # Compute delta as numerical distance from target region
        delta = int(h, 16) >> (256 - len(target_prefix)*4)  # reduce to prefix domain
        phase = int((math.cos(OMEGA * n) + 1) / 2 * (2**32 - 1))

        # Phase-shifted XOR feedback
        tuned_nonce = nonce ^ (delta & phase)
        next_nonce = nonce - int(K_FEEDBACK * (nonce ^ tuned_nonce))

        history.append((n, nonce, h[:12]))
        nonce = next_nonce & 0xFFFFFFFF  # stay in 32-bit range

    print("NOT FOUND in", max_iter, "iterations")
    return None, None, history

# --- Example use with real values ---
# Use a known Bitcoin block hash and a merkle root (just real-looking data for context)
if __name__ == '__main__':
    # Real Bitcoin block 100000 values
    prev_hash = "00000000b873e79784647a6c82962c70d228557d24a747ea4d1b8bbe878e1206"
    merkle_root = "e320b6c2fffc8baae3aa89b16cba141c3bf3c3f3b3a3d3e3f3c3f3f3c3f3c3f3"

    recursive_sha_unfold(prev_hash, merkle_root)
