import numpy as np
import matplotlib.pyplot as plt

# Stream generators
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

# Spaced-checkpoint compressor with wide motif window
def spaced_checkpoint_compression(stream, base_eps, min_len, max_len, bias, m, k=0.35):
    deltas = np.diff(stream)
    N_main = len(deltas)
    journal = [d for d in deltas if abs(d) > base_eps]
    motif_dict = {}
    next_id = 0
    tokens = []
    storage_tokens = 0
    checkpoints = 0
    since_cp = m  # allow immediate checkpoint

    pos = 0
    while pos < len(journal):
        matched = False
        for L in range(max_len, min_len - 1, -1):
            if pos + L <= len(journal):
                seq = tuple(journal[pos:pos+L])
                if seq in motif_dict:
                    tokens.append(('M', motif_dict[seq]))
                    storage_tokens += 1
                    since_cp += 1
                    pos += L
                    matched = True
                    break
        if matched:
            if since_cp >= m:
                checkpoints += 1
                since_cp = 0
            continue

        # literal
        tokens.append(('D', journal[pos]))
        storage_tokens += 1
        since_cp += 1

        # new motif insertion
        candidates = [
            tuple(journal[pos:pos+L])
            for L in range(min_len, max_len+1) if pos+L <= len(journal)
        ]
        if candidates and bias:
            best = min(candidates, key=lambda s: abs(np.mean(s) - k))
            if best not in motif_dict:
                motif_dict[best] = next_id
                next_id += 1

        if since_cp >= m:
            checkpoints += 1
            since_cp = 0

        pos += 1

    storage_ratio = (storage_tokens + checkpoints) / N_main

    # Reconstruction for RMSE
    reconstructed = np.zeros(len(stream))
    idx = 1
    inv = {v: k for k, v in motif_dict.items()}
    for typ, val in tokens:
        if typ == 'M':
            for d in inv[val]:
                reconstructed[idx] = reconstructed[idx-1] + d
                idx += 1
        else:
            reconstructed[idx] = reconstructed[idx-1] + val
            idx += 1

    rmse = np.sqrt(np.mean((reconstructed - stream)**2))
    return storage_ratio, rmse

# Parameters
epsilons = [0.1, 0.2, 0.5]
motif_min_lengths = range(2, 9)
m = 50
bias = True

# Streams
streams = {
    'Sinusoid': generate_sinusoid(),
    'Bursty': generate_bursty()
}

# Run simulation and collect
results = {name: {eps: [] for eps in epsilons} for name in streams}
for name, stream in streams.items():
    for eps in epsilons:
        for min_len in motif_min_lengths:
            max_len = min_len + 8
            sr, rmse = spaced_checkpoint_compression(stream, eps, min_len, max_len, bias, m)
            results[name][eps].append((sr, rmse))

# Plotting
fig, axs = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
for ax, (name, eps_data) in zip(axs, results.items()):
    for eps in epsilons:
        pts = np.array(eps_data[eps])
        ax.plot(pts[:,0], pts[:,1], marker='o', linestyle='-', label=f'ε={eps}')
    ax.set_title(name)
    ax.set_xlabel('Storage Ratio')
    ax.set_ylabel('RMSE')
    ax.grid(True)
    ax.legend()
plt.suptitle('Storage vs RMSE for ε ∈ {0.1,0.2,0.5} (m=50, wide motifs)')
plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.show()
