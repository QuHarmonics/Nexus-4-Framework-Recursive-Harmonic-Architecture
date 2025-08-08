###############################################################################
# π-CARRIER SHA-256 BREATHING DEMODULATOR — Full Version (with 2D Heatmap)
###############################################################################

import math
import numpy as np
import matplotlib.pyplot as plt
from hashlib import sha256

# ─────────────────────────────────────────────────────────────────────────────
# 1. Helpers
# ─────────────────────────────────────────────────────────────────────────────
def sha256_bits(message: bytes) -> list[int]:
    """Return the 256-bit digest as a list of ±1 values."""
    digest = sha256(message).digest()
    bits = [(bit ^ 1) * -1 + 1          # 0 → –1, 1 → +1
            for byte in digest
            for bit in [(byte >> i) & 1 for i in range(7, -1, -1)]]
    return bits                         # len == 256, values in {-1, +1}

def carrier(phi: float) -> list[float]:
    """Pre-compute sin(π·i + φ) for i = 1 … 256."""
    return [math.sin(math.pi * i + phi) for i in range(1, 257)]

def segment_indices(bits_len: int, seg_len: int) -> list[tuple[int, int]]:
    """Return (start, end) index pairs for non-overlapping segments."""
    return [(k, min(k + seg_len, bits_len)) for k in range(0, bits_len, seg_len)]

# ─────────────────────────────────────────────────────────────────────────────
# 2. Global phase-lock
# ─────────────────────────────────────────────────────────────────────────────
def best_phase(bits: list[int], sweep_step: float = 0.001) -> float:
    """Find φ ∈ [0, 2π) maximizing Σ bits[i] · sin(π·i + φ)."""
    best_phi, best_r = 0.0, -1.0
    phi = 0.0
    while phi < 2 * math.pi:
        r = sum(b * s for b, s in zip(bits, carrier(phi)))
        if r > best_r:
            best_phi, best_r = phi, r
        phi += sweep_step
    return best_phi

# ─────────────────────────────────────────────────────────────────────────────
# 3. First-pass drift map
# ─────────────────────────────────────────────────────────────────────────────
def drift_map(bits: list[int], phi: float, seg_len: int = 16) -> list[float]:
    """Phase offset (deg) per segment relative to aligned carrier."""
    c = carrier(phi)
    phase = []
    for a, b in segment_indices(len(bits), seg_len):
        if b > len(bits): break  # Safety
        r = sum(bits[i] * c[i] for i in range(a, b))
        norm = max(min(r / (b - a), 1.0), -1.0)
        phase.append(math.degrees(math.asin(norm)))   # -90° … +90°
    return phase

# ─────────────────────────────────────────────────────────────────────────────
# 4. Recursive refinement
# ─────────────────────────────────────────────────────────────────────────────
def refine(bits: list[int], phi: float,
           coarse_len: int = 16, fine_len: int = 4) -> dict:
    """
    Return a nested dict:
        {segment_id: {'phase': θ_deg,
                      'sub': {sub_id: θ_deg, …}}, …}
    """
    c_phase = drift_map(bits, phi, coarse_len)
    segments = segment_indices(len(bits), coarse_len)
    result = {}
    for idx, θ in enumerate(c_phase):
        a, b = segments[idx]
        segment_bits = bits[a:b]
        if len(segment_bits) == 0:
            continue  # skip empty
        fine_phase = drift_map(segment_bits, phi + math.radians(θ), fine_len)
        result[idx] = {
            'phase': θ,
            'sub': {j: θ_sub for j, θ_sub in enumerate(fine_phase)}
        }
    return result

# ─────────────────────────────────────────────────────────────────────────────
# 5. Top-level decoder
# ─────────────────────────────────────────────────────────────────────────────
def decode(message: bytes,
           coarse_len: int = 16,
           fine_len:   int = 4) -> dict:
    bits = sha256_bits(message)
    φ = best_phase(bits)                        # global lock
    return refine(bits, φ, coarse_len, fine_len) # one refinement level

# ─────────────────────────────────────────────────────────────────────────────
# 6. Visualization — Heatmap of breathing
# ─────────────────────────────────────────────────────────────────────────────
def plot_drift_heatmap(drift: dict, title: str = "SHA-256 Breathing Map"):
    """Visualize coarse + fine drift maps as a 2D heatmap."""
    # Build full fine drift matrix
    grid = []
    for coarse_id in sorted(drift.keys(), key=lambda x: int(x)):
        fine_segments = drift[coarse_id]['sub']
        fine_row = [fine_segments[j] for j in sorted(fine_segments.keys())]
        grid.append(fine_row)
    
    grid = np.array(grid)

    plt.figure(figsize=(8, 6))
    plt.imshow(grid, cmap="seismic", interpolation='nearest', aspect='auto')
    plt.colorbar(label="Phase Drift (degrees)")
    plt.title(title)
    plt.xlabel("Fine Sub-segment ID")
    plt.ylabel("Coarse Segment ID")
    plt.grid(False)
    plt.show()

# ─────────────────────────────────────────────────────────────────────────────
# 7. Example usage
# ─────────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    import json

    # ⚡️ Try different messages!
    message = "Hello, π-carrier!".encode('utf-8')  # Must encode UTF-8!

    # Decode π-carrier drift map
    drift = decode(message)

    # Pretty print drift structure
    print("\nDrift Structure:")
    print(json.dumps(drift, indent=2))

    # Plot breathing heatmap
    plot_drift_heatmap(drift, title=f"SHA-256 Breathing Map for: {message.decode('utf-8')}")
