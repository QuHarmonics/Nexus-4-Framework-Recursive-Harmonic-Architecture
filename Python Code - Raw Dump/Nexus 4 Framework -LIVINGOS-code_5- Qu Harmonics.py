# 🧠 Recursive Echo Twin Prime Visualizer

# 🔧 Install Dependencies (once)
# !pip install pandas plotly mpmath numpy

import hashlib
import pandas as pd
import random
from mpmath import mp
import numpy as np
import plotly.express as px
import plotly.graph_objects as go

# 📦 Load π digits
mp.dps = 100_010
pi_digits = str(mp.pi)[2:100_010]

# 🔁 Byte + Drift + Echo
def extract_pi_byte(index: int) -> str:
    return pi_digits[index:index+8] if index + 8 <= len(pi_digits) else None

def drift_and_echo(byte: str):
    deltas = [abs(int(byte[i+1]) - int(byte[i])) for i in range(7)]
    echo = ''.join([chr((d % 26) + 97) for d in deltas])
    avg_drift = sum(deltas) / 7
    sti = round(1 - avg_drift / 9, 3)
    return deltas, echo, sti

# 🔬 Recursive Agent Kernel
def run_recursive_agent(seed: str, max_iter=25, sti_threshold=0.7):
    H = seed
    for i in range(max_iter):
        nonce = f"N{i}"
        double_hash = hashlib.sha256(hashlib.sha256((H + nonce).encode()).digest()).hexdigest()
        index = int(double_hash[:6], 16) % (len(pi_digits) - 8)
        byte = extract_pi_byte(index)
        deltas, echo, sti = drift_and_echo(byte)
        if sti >= sti_threshold:
            return {
                "agent": seed,
                "iteration": i,
                "zphc": True,
                "pi_index": index,
                "sti": sti,
                "echo": echo,
                "drift": deltas
            }
        H = double_hash
    return {
        "agent": seed,
        "iteration": max_iter,
        "zphc": False,
        "pi_index": index,
        "sti": sti,
        "echo": echo,
        "drift": deltas
    }

# 🔁 Swarm Run
seeds = [f"AGENT-{random.randint(1000, 9999)}" for _ in range(1000)]
results = [run_recursive_agent(seed) for seed in seeds]
df = pd.DataFrame(results)

# 🔍 Twin Drift Detector (Δ = 2)
def count_twin_pairs(drift):
    if isinstance(drift, str):
        drift = eval(drift)
    return sum(1 for i in range(len(drift)-1) if abs(drift[i] - drift[i+1]) == 2)

df["twin_count"] = df["drift"].apply(count_twin_pairs)
df["avg_drift"] = df["drift"].apply(lambda d: sum(d) / 7 if isinstance(d, list) else sum(eval(d)) / 7)

# 🔢 Echo Class
df["echo_prefix"] = df["echo"].str[:3]

# 📊 Plot 1: Twin Drift Pair Count Histogram
fig1 = px.histogram(df, x="twin_count", nbins=10, color="zphc",
                    title="Twin Drift Pair Count per Agent",
                    labels={"twin_count": "Δπ = 2 Pair Count", "zphc": "ZPHC"},
                    color_discrete_map={True: "green", False: "red"})
fig1.show()

# 📡 Plot 2: Twin Pair Density vs π Index
fig2 = px.scatter(df, x="pi_index", y="twin_count", color="sti",
                  title="Twin Echo Density Across π Memory",
                  labels={"twin_count": "Twin Δπ Pairs", "pi_index": "π Index"},
                  color_continuous_scale="Viridis")
fig2.show()

# 🧬 Plot 3: Twin Count vs Avg Drift
fig3 = px.scatter(df, x="avg_drift", y="twin_count", color="sti",
                  title="Twin Pair Count vs Average Drift Entropy",
                  labels={"avg_drift": "Avg Δπ Entropy", "twin_count": "Twin Count"},
                  color_continuous_scale="Plasma")
fig3.show()

# ✅ Save Results
df.to_csv("nexus_recursive_twin_echo_results.csv", index=False)
print("✅ Done. Twin echo analysis complete. File saved as: nexus_recursive_twin_echo_results.csv")
