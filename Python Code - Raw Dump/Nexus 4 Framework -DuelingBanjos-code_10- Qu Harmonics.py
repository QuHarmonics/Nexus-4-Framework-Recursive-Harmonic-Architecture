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

# Modified compression: checkpoint only on new motifs of length >= 6
def tame_novelty_motif_compression(stream, base_eps, min_len, max_len, k=0.35):
    deltas = np.diff(stream)
    N_main = len(deltas)
    journal = [d for d in deltas if abs(d) > base_eps]
    
    motif_dict = {}
    next_id = 0
    tokens = []
    storage_tokens = 0
    checkpoints = 0
    pos = 0

    while pos < len(journal):
        matched = False
        # match existing motifs
        for L in range(max_len, min_len-1, -1):
            if pos + L <= len(journal):
                seq = tuple(journal[pos:pos+L])
                if seq in motif_dict:
                    tokens.append(('M', motif_dict[seq]))
                    storage_tokens += 1
                    pos += L
                    matched = True
                    break
        if matched:
            continue

        # literal
        tokens.append(('D', journal[pos]))
        storage_tokens += 1

        # bias selection
        candidates = [tuple(journal[pos:pos+L])
                      for L in range(min_len, max_len+1)
                      if pos+L <= len(journal)]
        if candidates:
            best = min(candidates, key=lambda s: abs(np.mean(s) - k))
            if best not in motif_dict:
                motif_dict[best] = next_id
                next_id += 1
                if len(best) >= 6:
                    checkpoints += 1
        pos += 1
    
    total_storage = (storage_tokens + checkpoints) / N_main

    # Reconstruction
    reconstructed = np.zeros(len(stream))
    idx = 1
    inv_motifs = {v: k for k, v in motif_dict.items()}
    for typ, val in tokens:
        if typ == 'M':
            for d in inv_motifs[val]:
                reconstructed[idx] = reconstructed[idx-1] + d
                idx += 1
        else:
            reconstructed[idx] = reconstructed[idx-1] + val
            idx += 1

    rmse = np.sqrt(np.mean((reconstructed - stream)**2))
    return total_storage, rmse

# Collect frontier points for each stream under tame checkpoint rule
streams = {
    'random': generate_random_walk(),
    'sinusoid': generate_sinusoid(),
    'bursty': generate_bursty()
}

motif_lengths = list(range(2, 9))
frontiers = {}

for name, stream in streams.items():
    pts = []
    for min_len in motif_lengths:
        max_len = min_len + 3
        sr, rmse = tame_novelty_motif_compression(stream, base_eps=1.0, min_len=min_len, max_len=max_len)
        pts.append((sr, rmse))
    frontiers[name] = np.array(pts)

# Fit power-law + offset: y ≈ A x^B + C
coeffs = {}
for name, data in frontiers.items():
    x = data[:,0]
    y = data[:,1]
    C = np.min(y)
    y_adj = y - C
    mask = y_adj > 0
    if np.sum(mask) >= 2:
        x_fit = x[mask]
        y_fit = y_adj[mask]
        ln_x = np.log(x_fit)
        ln_y = np.log(y_fit)
        B, lnA = np.polyfit(ln_x, ln_y, 1)
        A = np.exp(lnA)
        coeffs[name] = (A, B, C)
    else:
        coeffs[name] = (np.nan, np.nan, C)

# Print fits
for name, (A, B, C) in coeffs.items():
    if np.isnan(A):
        print(f"{name.capitalize()}: Insufficient variation for fit (constant RMSE={C:.3f})")
    else:
        print(f"{name.capitalize()}: RMSE ≈ {A:.3f}·r^{B:.3f} + {C:.3f}")

# Plot frontiers and fits
plt.figure(figsize=(8,6))
colors = {'random':'red','sinusoid':'blue','bursty':'green'}
for name, data in frontiers.items():
    x, y = data[:,0], data[:,1]
    A, B, C = coeffs[name]
    plt.scatter(x, y, color=colors[name], label=f'{name} data')
    if not np.isnan(A):
        x_line = np.linspace(min(x), max(x), 100)
        y_line = A * (x_line**B) + C
        plt.plot(x_line, y_line, '--', color=colors[name], label=f'{name} fit')
plt.xlabel('Storage Ratio')
plt.ylabel('RMSE')
plt.title('Tamed Checkpoint (len≥6) Frontier & Fits (ε=1.0)')
plt.legend()
plt.grid(True)
plt.show()
