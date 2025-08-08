import numpy as np
import hashlib

def sha_curvature_resonance(seed, nonces=1000, k=0.1, attractor=0.35):
    """Compute curvature resonance of SHA-256 digests of (seed+nonce) strings."""
    resonances = []
    prev_digest = None
    for n in range(nonces):
        s = f"{seed}{n}"
        h = hashlib.sha256(s.encode('utf-8')).digest()
        val = int.from_bytes(h[:4], "big") / 2**32  # Normalize to [0,1)
        if prev_digest is not None:
            diff = val - prev_digest
            curvature = abs(diff)
            # Resonance feedback
            delta = abs(curvature - attractor)
            resonance = curvature / (1 + k * delta)
            resonances.append(resonance)
        prev_digest = val
    return np.array(resonances)

def pi_chunk_resonance(pi_digits, chunk_size=8, mode='sum', attractor=0.35):
    """Folding/aggregation of pi digits (as string or array) for resonance calculation."""
    chunks = [pi_digits[i:i+chunk_size] for i in range(0, len(pi_digits)-chunk_size+1, chunk_size)]
    vals = []
    for chunk in chunks:
        digits = np.array([int(d) for d in chunk])
        if mode == 'sum':
            val = digits.sum() / (chunk_size * 9)  # Normalize sum to [0,1]
        elif mode == 'mirror_sum':
            val = (digits + digits[::-1]).sum() / (2 * chunk_size * 9)
        elif mode == 'product':
            prod = np.prod(digits + 1)  # Avoid zeros
            val = np.log(prod) / (chunk_size * np.log(10))  # Normalize log-product
        else:
            raise ValueError("Unknown mode")
        delta = abs(val - attractor)
        resonance = val / (1 + 0.1 * delta)  # Feedback as above
        vals.append(resonance)
    return np.array(vals)

def harmonic_window_stats(resonances, attractor=0.35, tol=0.05):
    """Compute % of values inside harmonic window, mean, std, entropy."""
    in_win = np.abs(resonances - attractor) <= tol
    mean = np.mean(resonances)
    std = np.std(resonances)
    p = np.histogram(resonances, bins=32, range=(0,1))[0]
    p = p / p.sum() + 1e-12
    entropy = -np.sum(p * np.log2(p))
    return {
        "count": len(resonances),
        "hits_in_window": in_win.sum(),
        "window_percent": 100 * in_win.sum() / len(resonances),
        "mean": mean,
        "std": std,
        "entropy": entropy
    }

# Example usage
# (Replace with true pi digits, e.g., via mpmath or preloaded string)
pi_digits = "14159265358979323846264338327950288419716939937510" * 100  # toy sample

sha_results = sha_curvature_resonance("14159265", nonces=10000, k=0.1)
pi_results = pi_chunk_resonance(pi_digits, chunk_size=8, mode='sum')

sha_stats = harmonic_window_stats(sha_results)
pi_stats = harmonic_window_stats(pi_results)

print("SHA-256 Curvature:", sha_stats)
print("Pi Digits Folded:", pi_stats)
