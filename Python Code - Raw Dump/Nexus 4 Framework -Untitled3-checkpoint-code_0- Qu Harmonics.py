###############################################################################
# π-CARRIER  SHA-256  BREATHING  DEMODULATOR  – SKELETON
###############################################################################

import math
from hashlib import sha256

# ─────────────────────────────────────────────────────────────────────────────
# 1. Helpers
# ─────────────────────────────────────────────────────────────────────────────
def sha256_bits(message: bytes) -> list[int]:
    """Return the 256-bit digest as a list of ±1 values."""
    digest = sha256(message).digest()
    bits = [(bit ^ 1) * -1 + 1        # 0 → –1, 1 → +1
            for byte in digest
            for bit in [(byte >> i) & 1 for i in range(7, -1, -1)]]
    return bits                      # len == 256, values in {-1, +1}

def carrier(phi: float) -> list[float]:
    """Pre-compute sin(π·i + φ) for i = 1 … 256."""
    return [math.sin(math.pi * i + phi) for i in range(1, 257)]

def segment_indices(seg_len: int) -> list[tuple[int, int]]:
    """Return (start, end) index pairs for non-overlapping segments."""
    return [(k, k + seg_len) for k in range(0, 256, seg_len)]

# ─────────────────────────────────────────────────────────────────────────────
# 2. Global phase-lock
# ─────────────────────────────────────────────────────────────────────────────
def best_phase(bits: list[int], sweep_step: float = 0.001) -> float:
    """Find φ ∈ [0, 2π) maximising Σ bits[i] · sin(π·i + φ)."""
    best_phi, best_r = 0.0, -1.0
    phi = 0.0
    while phi < 2 * math.pi:
        r = sum(b * s for b, s in zip(bits, carrier(phi)))
        if r > best_r:
            best_phi, best_r = phi, r
        phi += sweep_step
    return best_phi

# ─────────────────────────────────────────────────────────────────────────────
# 3. First-pass drift map  (16-bit slices by default)
# ─────────────────────────────────────────────────────────────────────────────
def drift_map(bits: list[int], phi: float, seg_len: int = 16) -> list[float]:
    """Phase offset (deg) per segment relative to aligned carrier."""
    c      = carrier(phi)
    phase  = []
    for a, b in segment_indices(seg_len):
        r   = sum(bits[i] * c[i] for i in range(a, b))
        # Normalise: r_max = seg_len (all bits match carrier sign)
        norm = max(min(r / seg_len, 1.0), -1.0)
        phase.append(math.degrees(math.asin(norm)))   # -90° … +90°
    return phase

# ─────────────────────────────────────────────────────────────────────────────
# 4. Recursive refinement  (simple two-level demo)
# ─────────────────────────────────────────────────────────────────────────────
def refine(bits: list[int], phi: float,
           coarse_len: int = 16, fine_len: int = 4) -> dict:
    """
    Return a nested dict:
        {segment_id: {'phase': θ_deg,
                      'sub': {sub_id: θ_deg, …}}, …}
    """
    c_phase = drift_map(bits, phi, coarse_len)
    result  = {}
    for idx, θ in enumerate(c_phase):
        a, b  = segment_indices(coarse_len)[idx]
        subφ  = phi + math.radians(θ)                 # local re-align
        subθs = drift_map(bits[a:b], subφ, fine_len)
        result[idx] = {
            'phase' : θ,
            'sub'   : {j: θ_sub for j, θ_sub in enumerate(subθs)}
        }
    return result

# ─────────────────────────────────────────────────────────────────────────────
# 5. Top-level callable
# ─────────────────────────────────────────────────────────────────────────────
def decode(message: bytes,
           coarse_len: int = 16,
           fine_len:   int = 4) -> dict:
    bits = sha256_bits(message)
    φ    = best_phase(bits)                        # global lock
    return refine(bits, φ, coarse_len, fine_len)   # one refinement level

# ─────────────────────────────────────────────────────────────────────────────
# EXAMPLE  (prints drift map for the string “Hello”)
# ─────────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    import json
    drift = decode(b"Hello")
    print(json.dumps(drift, indent=2))
