import numpy as np
import matplotlib.pyplot as plt

# Simulator parameters
length = 1000
epsilons = [0.5, 1, 2, 4]
window_vol = 10  # for adaptive epsilon
k = 0.35  # harmonic constant for replay drift
max_motif = 5
min_motif = 2

def generate_random_walk(n):
    np.random.seed(0)
    return np.cumsum(np.random.randn(n))

def novelty_checkpoint_compression(stream, base_eps, adaptive=False):
    # Compute deltas
    deltas = np.diff(stream)
    N = len(deltas)
    
    # Volatility history for adaptive epsilon
    vol_hist = [np.std(deltas[max(0, i-window_vol):i+1]) for i in range(N)]
    
    # Journal entries
    journal = []
    epsilons_n = []
    for i, d in enumerate(deltas):
        eps = (base_eps * vol_hist[i]) if adaptive else base_eps
        if abs(d) > eps:
            journal.append(d)
            epsilons_n.append(eps)
    
    # Motif compression with novelty-checkpoint
    motif_dict = set()
    storage_tokens = 0
    checkpoints = 0
    tokens = []
    pos = 0
    # Flatten journal list for processing
    journal = list(journal)
    while pos < len(journal):
        matched = False
        # Try to match existing motifs (longest first)
        for L in range(max_motif, min_motif - 1, -1):
            if pos + L <= len(journal):
                seq = tuple(journal[pos:pos+L])
                if seq in motif_dict:
                    tokens.append(('motif', seq))
                    storage_tokens += 1
                    pos += L
                    matched = True
                    break
        if not matched:
            # Emit literal
            d = journal[pos]
            tokens.append(('lit', d))
            storage_tokens += 1
            # Insert new motifs (seqs length 2..5) and checkpoint on first insertion
            inserted = False
            for L in range(min_motif, max_motif+1):
                if pos + L <= len(journal):
                    seq = tuple(journal[pos:pos+L])
                    if seq not in motif_dict:
                        motif_dict.add(seq)
                        if not inserted:
                            checkpoints += 1
                            inserted = True
            pos += 1
    
    total_storage = storage_tokens + checkpoints
    
    # Reconstruction for RMSE
    reconstructed = np.zeros(len(stream))
    # We'll reconstruct context from last checkpoint splittings
    # For simplicity, we replay from start, resetting at checkpoints:
    idx = 1  # stream index for delta application
    jpos = 0
    for token in tokens:
        typ, val = token
        if typ == 'motif':
            for d in val:
                reconstructed[idx] = reconstructed[idx-1] + d
                idx += 1
        else:  # literal
            d = val
            reconstructed[idx] = reconstructed[idx-1] + d
            idx += 1
    # Compute RMSE over deltas in stream vs reconstructed
    true_vals = stream
    rmse = np.sqrt(np.mean((reconstructed - true_vals)**2))
    
    return total_storage / len(deltas), rmse

# Run simulations for random walk
stream = generate_random_walk(length)
results = {'fixed': [], 'adaptive': []}

for adaptive in [False, True]:
    key = 'adaptive' if adaptive else 'fixed'
    for eps in epsilons:
        storage_ratio, rmse = novelty_checkpoint_compression(stream, eps, adaptive)
        results[key].append((storage_ratio, rmse))

# Plotting
plt.figure(figsize=(6,5))
for tag, marker in [('fixed', 'o'), ('adaptive', 'x')]:
    ratios, rmses = zip(*results[tag])
    plt.scatter(ratios, rmses, label=tag, marker=marker)
plt.title('Novelty-Checkpoint: RMSE vs Total Storage Ratio (Random Walk)')
plt.xlabel('Total Storage Ratio (#checkpoints + #tokens) / stream length')
plt.ylabel('RMSE over entire stream')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
