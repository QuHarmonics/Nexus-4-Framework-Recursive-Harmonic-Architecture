# Full Prototype: Spaced-Checkpoint Motif Compressor with Power-Law Fits

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# -----------------------------------------------------------------------------
# 1) STREAM GENERATORS
# -----------------------------------------------------------------------------
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

# -----------------------------------------------------------------------------
# 2) SPACED-CHECKPOINT COMPRESSION FUNCTION
# -----------------------------------------------------------------------------
def spaced_checkpoint_compression(stream, base_eps, m, min_len, max_len, k=0.35):
    """
    - base_eps: journal threshold (fixed)
    - m: minimum number of journal-entries between full checkpoints
    - min_len, max_len: motif window bounds
    - k: harmonic bias target
    Returns (storage_ratio, rmse)
    """
    deltas = np.diff(stream)
    N = len(deltas)

    # 2a) Journal only “significant” deltas
    journal = [d for d in deltas if abs(d) > base_eps]

    motif_dict = {}   # seq -> id
    next_id = 0
    tokens = []       # list of ("M",id) or ("D",delta)
    storage_tokens = 0
    checkpoints = 0
    since_last_cp = m   # allow immediate checkpoint

    pos = 0
    while pos < len(journal):
        # 2b) Try matching existing motifs (longest first)
        matched = False
        for L in range(max_len, min_len - 1, -1):
            if pos + L <= len(journal):
                seq = tuple(journal[pos:pos+L])
                if seq in motif_dict:
                    tokens.append(("M", motif_dict[seq]))
                    storage_tokens += 1
                    since_last_cp += 1
                    pos += L
                    matched = True
                    break
        if matched:
            continue

        # 2c) Emit literal + potential new motif
        tokens.append(("D", journal[pos]))
        storage_tokens += 1
        since_last_cp += 1

        # Candidate motifs
        candidates = [
            tuple(journal[pos:pos+L])
            for L in range(min_len, max_len + 1)
            if pos + L <= len(journal)
        ]
        if candidates:
            # Pick motif whose mean Δ is closest to k
            best = min(candidates, key=lambda s: abs(np.mean(s) - k))
            if best not in motif_dict:
                motif_dict[best] = next_id
                next_id += 1
                # Spaced checkpoint
                if since_last_cp >= m:
                    checkpoints += 1
                    since_last_cp = 0

        pos += 1

    storage_ratio = (storage_tokens + checkpoints) / N

    # 2d) Reconstruct and compute RMSE
    reconstructed = np.zeros(len(stream))
    idx = 1
    inv = {v: k for k, v in motif_dict.items()}
    for typ, val in tokens:
        if typ == "M":
            for d in inv[val]:
                reconstructed[idx] = reconstructed[idx-1] + d
                idx += 1
        else:
            reconstructed[idx] = reconstructed[idx-1] + val
            idx += 1

    rmse = np.sqrt(np.mean((reconstructed - stream)**2))
    return storage_ratio, rmse

# -----------------------------------------------------------------------------
# 3) SWEEP PARAMETERS AND COLLECT RESULTS
# -----------------------------------------------------------------------------
epsilon = 1.0
motif_lengths = list(range(2, 9))
m_values = [25, 50, 100]

streams = {
    'random': generate_random_walk(),
    'sinusoid': generate_sinusoid(),
    'bursty': generate_bursty()
}

rows = []
for name, stream in streams.items():
    for m in m_values:
        for min_len in motif_lengths:
            max_len = min_len + 3
            sr, rmse = spaced_checkpoint_compression(
                stream, base_eps=epsilon,
                m=m, min_len=min_len, max_len=max_len
            )
            rows.append({
                'stream': name,
                'm': m,
                'min_len': min_len,
                'storage_ratio': sr,
                'rmse': rmse
            })
df = pd.DataFrame(rows)

# -----------------------------------------------------------------------------
# 4) FIT POWER-LAW + OFFSET PER (stream, m)
# -----------------------------------------------------------------------------
fit_rows = []
for (name, m), group in df.groupby(['stream','m']):
    x = group['storage_ratio'].values
    y = group['rmse'].values
    C = y.min()
    y_adj = y - C
    mask = y_adj > 0
    if mask.sum() >= 2:
        ln_x = np.log(x[mask])
        ln_y = np.log(y_adj[mask])
        B, lnA = np.polyfit(ln_x, ln_y, 1)
        A = np.exp(lnA)
        y_pred = A * (x**B) + C
        ss_res = np.sum((y - y_pred)**2)
        ss_tot = np.sum((y - y.mean())**2)
        R2 = 1 - ss_res/ss_tot
    else:
        A = B = R2 = np.nan
    fit_rows.append({'stream': name, 'm': m, 'A': A, 'B': B, 'C': C, 'R2': R2})

fits = pd.DataFrame(fit_rows)
print("Power-law + offset fits:")
print(fits)

# -----------------------------------------------------------------------------
# 5) VISUALIZE FRONTIERS & FITS
# -----------------------------------------------------------------------------
plt.figure(figsize=(8,6))
colors = {'random':'red','sinusoid':'blue','bursty':'green'}
markers = {25:'o',50:'s',100:'^'}

for _, row in fits.iterrows():
    name, m = row['stream'], row['m']
    A, B, C = row['A'], row['B'], row['C']
    subset = df[(df['stream']==name)&(df['m']==m)]
    x, y = subset['storage_ratio'], subset['rmse']
    plt.scatter(x, y, color=colors[name], marker=markers[m],
                label=f'{name}, m={m}')
    if not np.isnan(A):
        x_line = np.linspace(x.min(), x.max(), 100)
        y_line = A * x_line**B + C
        plt.plot(x_line, y_line, '--', color=colors[name])

plt.xlabel('Storage Ratio')
plt.ylabel('RMSE')
plt.title('Spaced-Checkpoint Frontiers & Power-Law Fits (ε=1.0)')
plt.legend(loc='best')
plt.grid(True)
plt.tight_layout()
plt.show()

