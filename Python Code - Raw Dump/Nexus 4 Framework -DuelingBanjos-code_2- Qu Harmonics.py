import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from IPython.display import display

# 1. Delta-stream generators
def generate_sinusoid(N=1000, period=100, noise_std=0.1):
    t = np.arange(N)
    x = np.sin(2 * np.pi * t / period) + np.random.normal(scale=noise_std, size=N)
    return np.diff(x, prepend=x[0])

def generate_random_walk(N=1000, step_std=1.0):
    return np.random.normal(loc=0.0, scale=step_std, size=N)

def generate_bursty(N=1000, noise_std=0.5, burst_prob=0.01, burst_mag=10.0):
    x = np.zeros(N)
    for i in range(1, N):
        if np.random.rand() < burst_prob:
            x[i] = x[i-1] + np.random.choice([-burst_mag, burst_mag])
        else:
            x[i] = x[i-1] + np.random.normal(scale=noise_std)
    return np.diff(x, prepend=x[0])

# 2. Journaling function returning stored deltas and flags
def journal_stream(deltas, scheme='fixed', base_eps=1.0, window=20):
    N = len(deltas)
    recent = []
    stored_flags = []
    journal = []
    for d in deltas:
        if scheme == 'fixed':
            eps = base_eps
        else:  # adaptive
            local_vol = np.std(recent[-window:]) if len(recent) >= window else (np.std(recent) if recent else 1.0)
            eps = base_eps * local_vol
        if abs(d) > eps:
            journal.append(d)
            stored_flags.append(True)
        else:
            stored_flags.append(False)
        recent.append(d)
    return journal, stored_flags

# 3. Motif compressor (LZ-style sliding window)
def compress_journal(journal, min_len=2, max_len=5):
    motifs = {}  # sequence -> id
    motif_id = 0
    tokens = []
    i = 0
    while i < len(journal):
        matched = False
        for L in range(max_len, min_len - 1, -1):
            if i + L <= len(journal):
                seq = tuple(journal[i:i+L])
                if seq in motifs:
                    tokens.append(('M', motifs[seq]))
                    i += L
                    matched = True
                    break
        if not matched:
            # Add new motif of length min_len if possible
            if i + min_len <= len(journal):
                seq = tuple(journal[i:i+min_len])
                if seq not in motifs:
                    motifs[seq] = motif_id
                    motif_id += 1
            tokens.append(('D', journal[i]))
            i += 1
    return tokens, motifs

# 4. Decompress tokens back into original journal list
def decompress_tokens(tokens, motifs):
    inv_motifs = {idx: seq for seq, idx in motifs.items()}
    decompressed = []
    for typ, val in tokens:
        if typ == 'D':
            decompressed.append(val)
        else:  # 'M'
            decompressed.extend(inv_motifs[val])
    return decompressed

# 5. Reconstruction RMSE from stored flags and decompressed journal
def reconstruction_rmse(deltas, stored_flags, decompressed):
    true_cumsum = np.cumsum(deltas)
    recon_sum = 0.0
    recon = []
    idx = 0
    for flag in stored_flags:
        if flag:
            recon_sum += decompressed[idx]
            idx += 1
        recon.append(recon_sum)
    recon = np.array(recon)
    rmse = np.sqrt(np.mean((recon - true_cumsum)**2))
    return rmse

# 6. Run simulation with motif compression
stream_funcs = {
    'sinusoid': generate_sinusoid,
    'random_walk': generate_random_walk,
    'bursty': generate_bursty
}
schemes = ['fixed', 'adaptive']
epsilon_bases = [0.5, 1.0, 2.0, 4.0]

results = []
for name, func in stream_funcs.items():
    deltas = func()
    for scheme in schemes:
        for base_eps in epsilon_bases:
            journal, flags = journal_stream(deltas, scheme, base_eps)
            tokens, motifs = compress_journal(journal)
            decompressed = decompress_tokens(tokens, motifs)
            rmse = reconstruction_rmse(deltas, flags, decompressed)
            results.append({
                'stream': name,
                'scheme': scheme,
                'base_eps': base_eps,
                'journal_ratio': len(journal)/len(deltas),
                'tokens_ratio': len(tokens)/len(deltas),
                'rmse': rmse
            })

df = pd.DataFrame(results)

# 7. Display results
display(df)

# 8. Plot RMSE vs tokens_ratio for each stream
for stream_name in df['stream'].unique():
    df_stream = df[df['stream'] == stream_name]
    plt.figure()
    for scheme in schemes:
        df_s = df_stream[df_stream['scheme'] == scheme]
        marker = 'o' if scheme == 'fixed' else 'x'
        plt.scatter(df_s['tokens_ratio'], df_s['rmse'], label=scheme, marker=marker)
    plt.xlabel("Compressed Journal Ratio")
    plt.ylabel("RMSE")
    plt.title(f"Motif Compression: RMSE vs Compressed Ratio ({stream_name.replace('_',' ').title()})")
    plt.legend()
    plt.show()
