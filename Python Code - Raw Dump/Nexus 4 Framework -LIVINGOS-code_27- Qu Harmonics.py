import hashlib
import binascii
from mpmath import mp

def sha256_digest(data: bytes) -> bytes:
    return hashlib.sha256(data).digest()

def hamming_distance(a: bytes, b: bytes) -> int:
    return sum(bin(x ^ y).count("1") for x, y in zip(a, b))

def get_pi_bytes(n_bytes, start_digit=0):
    """Return n_bytes of pi, starting at start_digit after the decimal."""
    mp.dps = n_bytes * 3  # digits per byte (≈2.41, round up)
    pi_str = str(mp.pi)[2:]  # skip "3."
    pi_digits = pi_str[start_digit:start_digit + n_bytes * 3]
    # pack three decimal digits into a byte (mod 256)
    return bytes([int(pi_digits[i:i+3]) % 256 for i in range(0, len(pi_digits), 3)])

def pi_injector(step, state, hashed, history):
    # Get a fresh vector of pi bytes each resonance, unique per step
    pi_bytes = get_pi_bytes(len(hashed), start_digit=step)
    return bytes([b ^ p for b, p in zip(hashed, pi_bytes)])

def resonance_validator(hamming, mark1_ratio=0.35, tol=0.08):
    return abs((hamming / 256) - mark1_ratio) < tol

def recursive_harmonic_identity(seed: bytes, rounds=512, inject_external_on_resonance=None, verbose=True):
    state = seed
    history = [state]
    resonance_hits = []

    for step in range(rounds):
        history_digest = sha256_digest(b"".join(history))
        next_state = bytes([s ^ h for s, h in zip(state, history_digest)])
        hashed = sha256_digest(sha256_digest(next_state))
        ham = hamming_distance(state, hashed)
        resonance = resonance_validator(ham)

        if verbose:
            print(f"Step {step:3d} | Hamming: {ham:3d} | {'✔️' if resonance else '❌'}")

        if resonance and inject_external_on_resonance is not None:
            injected = inject_external_on_resonance(step, state, hashed, history)
            hashed = sha256_digest(injected)
            print("    [Injected real π bytes at resonance]")

        state = hashed
        history.append(state)

    print("\nFinal SHA state:", binascii.hexlify(state).decode())

if __name__ == "__main__":
    seed = bytes([0] * 32)
    recursive_harmonic_identity(
        seed,
        rounds=256,
        inject_external_on_resonance=pi_injector,
        verbose=True
    )
