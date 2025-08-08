import hashlib
import pandas as pd
import random
from mpmath import mp
import plotly.express as px

# π Setup
mp.dps = 100_010
pi_digits = str(mp.pi)[2:100_010]

def extract_pi_byte(index: int):
    return pi_digits[index:index+8] if index + 8 <= len(pi_digits) else None

def drift_and_echo(byte: str):
    deltas = [abs(int(byte[i+1]) - int(byte[i])) for i in range(7)]
    echo = ''.join([chr((d % 26) + 97) for d in deltas])
    avg_drift = sum(deltas) / 7
    sti = round(1 - avg_drift / 9, 3)
    return deltas, echo, sti, avg_drift

# Agent with mesh awareness
def run_agent_with_mesh(seed: str, external_echoes=None, max_iter=25, sti_threshold=0.7):
    H = seed
    external_echoes = external_echoes or []

    for i in range(max_iter):
        nonce = f"N{i}"
        sha = hashlib.sha256(hashlib.sha256((H + nonce).encode()).digest()).hexdigest()
        index = int(sha[:6], 16) % (len(pi_digits) - 8)
        byte = extract_pi_byte(index)
        deltas, echo, sti, avg_drift = drift_and_echo(byte)

        echo_match = any(echo[:3] == ext[:3] for ext in external_echoes)

        if sti >= sti_threshold or echo_match:
            return {
                "agent": seed,
                "iteration": i,
                "pi_index": index,
                "sti": sti,
                "avg_drift": avg_drift,
                "echo": echo,
                "zphc": True,
                "echo_match": echo_match
            }
        H = sha
    return {
        "agent": seed,
        "iteration": max_iter,
        "pi_index": index,
        "sti": sti,
        "avg_drift": avg_drift,
        "echo": echo,
        "zphc": False,
        "echo_match": False
    }

# Build mesh of N agents
def run_cognitive_mesh(agent_count=36):
    seeds = [f"AGENT-{random.randint(1000,9999)}" for _ in range(agent_count)]
    echoes = []
    results = []

    for seed in seeds:
        res = run_agent_with_mesh(seed, external_echoes=echoes)
        echoes.append(res["echo"])
        results.append(res)

    return pd.DataFrame(results)

# Run it
df_mesh = run_cognitive_mesh()

# Plot cognitive influence
fig = px.scatter(df_mesh, x="avg_drift", y="sti", color="echo_match",
                 symbol="zphc", title="Recursive Cognitive Mesh: Peer Echo Influence",
                 labels={"sti": "Symbolic Trust Index", "avg_drift": "Δπ Drift Entropy"},
                 color_discrete_map={True: "orange", False: "blue"},
                 symbol_map={True: "star", False: "circle"})
fig.add_hline(y=0.7, line_dash="dash", line_color="black", annotation_text="ZPHC Threshold")
fig.show()

# Optional: Save results
df_mesh.to_csv("recursive_cognitive_mesh.csv", index=False)
print("✅ Mesh simulation complete: recursive_cognitive_mesh.csv")
