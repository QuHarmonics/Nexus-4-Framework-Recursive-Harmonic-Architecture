import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# For interactive display
from ace_tools import display_dataframe_to_user

# Simulation functions
def generate_sinusoid(N=1000, period=100, noise_std=0.1):
    t = np.arange(N)
    x = np.sin(2 * np.pi * t / period) + np.random.normal(scale=noise_std, size=N)
    return np.diff(x, prepend=x[0])

def generate_random_walk(N=1000, step_std=1.0):
    steps = np.random.normal(loc=0.0, scale=step_std, size=N)
    return steps  # Δ = step

def generate_bursty(N=1000, noise_std=0.5, burst_prob=0.01, burst_magnitude=10.0):
    x = np.zeros(N)
    for i in range(1, N):
        if np.random.rand() < burst_prob:
            x[i] = x[i-1] + np.random.choice([-burst_magnitude, burst_magnitude])
        else:
            x[i] = x[i-1] + np.random.normal(scale=noise_std)
    return np.diff(x, prepend=x[0])

# Journaling and reconstruction
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

        # Journal decision
        if abs(d) > eps:
            journal.append(d)
            recon_sum += d  # apply stored delta
        recent.append(d)
        recon[i] = recon_sum

    true = np.cumsum(deltas)
    rmse = np.sqrt(np.mean((recon - true)**2))
    journal_size = len(journal)

    return journal_size, rmse

# Experiment setup
stream_funcs = {
    'sinusoid': generate_sinusoid,
    'random_walk': generate_random_walk,
    'bursty': generate_bursty
}
schemes = ['fixed', 'adaptive']
epsilon_bases = [0.5, 1.0, 2.0, 4.0]

# Run simulations
results = []
for stream_name, func in stream_funcs.items():
    deltas = func()
    for scheme in schemes:
        for base in epsilon_bases:
            size, error = journal_and_reconstruct(deltas, scheme, base)
            results.append({
                'stream': stream_name,
                'scheme': scheme,
                'base_eps': base,
                'journal_size': size,
                'journal_ratio': size / len(deltas),
                'rmse': error
            })

df = pd.DataFrame(results)

# Display the results table
display_dataframe_to_user("Delta Memory Simulator Results", df)

# Plot for one example: random_walk
df_rw = df[df['stream'] == 'random_walk']
plt.figure()
plt.scatter(df_rw['journal_ratio'], df_rw['rmse'], label='fixed', marker='o')
plt.scatter(df_rw['journal_ratio'][df_rw['scheme']=='adaptive'], df_rw['rmse'][df_rw['scheme']=='adaptive'], label='adaptive', marker='x')
plt.xlabel("Journal Ratio")
plt.ylabel("RMSE")
plt.title("RMSE vs Journal Ratio (Random Walk)")
plt.legend()
plt.show()
