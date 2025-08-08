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

# 2. Journaling & reconstruction function
def journal_and_reconstruct(deltas, scheme='fixed', base_eps=1.0, window=20):
    N = len(deltas)
    journal = []
    recon = np.zeros(N)
    recon_sum = 0.0
    recent = []

    for i, d in enumerate(deltas):
        # Determine epsilon
        if scheme == 'fixed':
            eps = base_eps
        elif scheme == 'adaptive':
            if len(recent) >= window:
                local_vol = np.std(recent[-window:])
            else:
                local_vol = np.std(recent) if recent else 1.0
            eps = base_eps * local_vol
        else:
            raise ValueError("Unknown scheme")

        # Journal decision & reconstruction
        if abs(d) > eps:
            journal.append(d)
            recon_sum += d
        recent.append(d)
        recon[i] = recon_sum

    true_cumsum = np.cumsum(deltas)
    rmse = np.sqrt(np.mean((recon - true_cumsum)**2))
    journal_ratio = len(journal) / N
    return len(journal), journal_ratio, rmse

# 3. Simulation parameters & run
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
            size, ratio, error = journal_and_reconstruct(deltas, scheme, base_eps)
            results.append({
                'stream': name,
                'scheme': scheme,
                'base_eps': base_eps,
                'journal_size': size,
                'journal_ratio': ratio,
                'rmse': error
            })

df = pd.DataFrame(results)

# 4. Display results table
display(df)

# 5. Plot RMSE vs Journal Ratio for each stream
for stream_name in df['stream'].unique():
    df_stream = df[df['stream'] == stream_name]
    plt.figure()
    for scheme in schemes:
        df_s = df_stream[df_stream['scheme'] == scheme]
        marker = 'o' if scheme == 'fixed' else 'x'
        plt.scatter(df_s['journal_ratio'], df_s['rmse'], label=scheme, marker=marker)
    plt.xlabel("Journal Ratio")
    plt.ylabel("RMSE")
    plt.title(f"RMSE vs Journal Ratio ({stream_name.replace('_',' ').title()})")
    plt.legend()
    plt.show()
