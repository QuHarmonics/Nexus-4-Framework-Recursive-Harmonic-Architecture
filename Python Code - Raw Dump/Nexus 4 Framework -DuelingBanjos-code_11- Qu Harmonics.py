import numpy as np
import matplotlib.pyplot as plt

# Stream generators
def generate_random_walk(N=1000, seed=0):
    np.random.seed(seed)
    return np.cumsum(np.random.randn(N))

def generate_sinusoid(N=1000, freq=5, amplitude=1.0, phase=0.0):
    t = np.linspace(0, 2 * np.pi, N)
    return amplitude * np.sin(freq * t + phase)

def generate_bursty(N=1000, p_burst=0.05, burst_scale=5.0, noise_scale=0.2, seed=42):
    np.random.seed(seed)
    deltas = []
    for _ in range(N):
        if np.random.rand() < p_burst:
            deltas.append(np.random.randn() * burst_scale)
        else:
            deltas.append(np.random.randn() * noise_scale)
    return np.cumsum(deltas)

# Tamed checkpoint: only length >=6
def tame_novelty_motif_compression(stream, base_eps, min_len, max_len, k=0.35):
    deltas = np.diff(stream)
    N_main = len(deltas)
    journal = [d for d in deltas if abs(d) > base_eps]
    
    motif_dict = {}
    next_id = 0
    storage_tokens = 0
    checkpoints = 0
    pos = 0

    while pos < len(journal):
        matched = False
        # match existing motifs
        for L in range(max_len, min_len-1, -1):
            if pos + L <= len(journal) and tuple(journal[pos:pos+L]) in motif_dict:
                storage_tokens += 1
                pos += L
                matched = True
                break
        if matched:
            continue

        # literal
        storage_tokens += 1

        # candidate motifs
        candidates = [tuple(journal[pos:pos+L])
                      for L in range(min_len, max_len+1)
                      if pos+L <= len(journal)]
        if candidates:
            # pick motif closest to k
            best = min(candidates, key=lambda s: abs(np.mean(s) - k))
            if best not in motif_dict:
                motif_dict[best] = next_id
                next_id += 1
                if len(best) >= 6:
                    checkpoints += 1
        pos += 1

    total_storage_ratio = (storage_tokens + checkpoints) / N_main

    # Simple approximate RMSE: deviation ~ constant C from prior fits
    # This is placeholder since full replay omitted for brevity
    # We'll use previously observed midpoints: random ~39, sinusoid ~0.7, bursty ~23
    return total_storage_ratio

# Recollect frontiers for ε=2.0
streams = {
    'random': generate_random_walk(),
    'sinusoid': generate_sinusoid(),
    'bursty': generate_bursty()
}
motif_lengths = list(range(2, 9))
eps = 2.0

frontiers = {}
for name, stream in streams.items():
    pts = []
    for min_len in motif_lengths:
        max_len = min_len + 3
        sr = tame_novelty_motif_compression(stream, base_eps=eps, min_len=min_len, max_len=max_len)
        # Use placeholder RMSE values: load from earlier frontiers if available
        pts.append((sr, np.nan))  # to be filled with actual RMSE
    frontiers[name] = np.array(pts)

# Display storage ratios
for name, data in frontiers.items():
    print(f"{name.capitalize()} storage ratios at ε={eps}: {data[:,0]}")
