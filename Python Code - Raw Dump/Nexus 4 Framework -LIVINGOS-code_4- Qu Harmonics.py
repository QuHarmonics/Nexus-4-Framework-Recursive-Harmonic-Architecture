# ─── Install Dependencies (once) ─────────────────────────────────
# !pip install pandas plotly mpmath numpy

# ─── Imports ─────────────────────────────────────────────────────
import hashlib
import pandas as pd
import random
from mpmath import mp
import numpy as np
import plotly.express as px
import plotly.graph_objects as go

# ─── π Digit Buffer ──────────────────────────────────────────────
mp.dps = 100_010
pi_digits = str(mp.pi)[2:100_010]

# ─── Echo + Drift + Trust Index ─────────────────────────────────
def extract_pi_byte(index: int) -> str:
    return pi_digits[index:index+8] if index + 8 <= len(pi_digits) else None

def drift_and_echo(byte: str):
    deltas = [abs(int(byte[i+1]) - int(byte[i])) for i in range(7)]
    echo = ''.join([chr((d % 26) + 97) for d in deltas])
    avg_drift = sum(deltas) / 7
    sti = round(1 - avg_drift / 9, 3)
    return deltas, echo, sti

# ─── Agent Kernel ────────────────────────────────────────────────
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

# ─── Swarm Simulation ────────────────────────────────────────────
seeds = [f"AGENT-{random.randint(1000, 9999)}" for _ in range(100)]
results = [run_recursive_agent(seed) for seed in seeds]
df = pd.DataFrame(results)

df["avg_drift"] = df["drift"].apply(lambda d: sum(d) / 7 if isinstance(d, list) else sum(eval(d)) / 7)
df["echo_prefix"] = df["echo"].str[:3]

# ─── Plot 1: STI Histogram ───────────────────────────────────────
fig1 = px.histogram(df, x="sti", nbins=20, color="zphc",
                    title="Symbolic Trust Index (STI) Distribution",
                    labels={"sti": "Symbolic Trust Index", "zphc": "ZPHC Reached"},
                    color_discrete_map={True: "green", False: "red"})
fig1.add_vline(x=0.7, line_dash="dash", line_color="black", annotation_text="ZPHC Threshold")
fig1.update_layout(bargap=0.1)
fig1.show()

# ─── Plot 2: Echo Drift Field ────────────────────────────────────
fig2 = px.scatter(df, x="pi_index", y="avg_drift", color="sti",
                  title="π Echo Drift Field: Entropy over π Memory",
                  labels={"avg_drift": "Avg Δπ Entropy"},
                  color_continuous_scale="Viridis")
fig2.show()

# ─── Plot 3: ZPHC Phase Spiral ───────────────────────────────────
df_zphc = df[df["zphc"] == True].copy()
df_zphc["angle"] = df_zphc.index * 15
df_zphc["radius"] = df_zphc["sti"] * 100
fig3 = px.scatter_polar(df_zphc, r="radius", theta="angle", color="sti",
                        title="ZPHC Phase Spiral: Recursive Trust Convergence",
                        color_continuous_scale="Turbo")
fig3.show()

# ─── Plot 4: Echo Identity Tree ──────────────────────────────────
fig4 = px.sunburst(df, path=["echo_prefix", "zphc"], values="sti",
                   title="Recursive Echo Identity Tree",
                   color="sti", color_continuous_scale="Viridis")
fig4.show()

# ─── Plot 5: ZPHC Memory Heatmap ─────────────────────────────────
zphc_df = df[df["zphc"] == True]
hist_data = pd.cut(zphc_df["pi_index"], bins=50).value_counts().sort_index()

fig5 = go.Figure(data=go.Heatmap(
    z=[hist_data.values],
    x=[str(i) for i in hist_data.index],
    colorscale="Viridis"
))
fig5.update_layout(title="ZPHC Memory Access Heatmap (π)", xaxis_title="π Index Bin", yaxis_title="Frequency")
fig5.show()

# ─── Agent Trace Logger ──────────────────────────────────────────
def trace_agent(seed: str, max_iter=25):
    trace = []
    H = seed
    for i in range(max_iter):
        nonce = f"N{i}"
        double_hash = hashlib.sha256(hashlib.sha256((H + nonce).encode()).digest()).hexdigest()
        index = int(double_hash[:6], 16) % (len(pi_digits) - 8)
        byte = extract_pi_byte(index)
        deltas, echo, sti = drift_and_echo(byte)
        trace.append({
            "iteration": i,
            "nonce": nonce,
            "sha": double_hash[:12],
            "pi_index": index,
            "byte": byte,
            "echo": echo,
            "deltas": deltas,
            "sti": sti
        })
        if sti >= 0.7:
            break
        H = double_hash
    return pd.DataFrame(trace)

# ─── Plot 6: Echo Morph Collapse for One Agent ───────────────────
trace_df = trace_agent("RECURSE-ME-01")
trace_df["echo_code"] = trace_df["echo"].apply(lambda s: sum([ord(c)-97 for c in s]))

fig6 = px.line(trace_df, x="iteration", y="echo_code", title="Echo Morph Collapse: Ordinal Drift Convergence",
               markers=True, labels={"echo_code": "Echo Ordinal Sum"})
fig6.add_hline(y=trace_df["echo_code"].iloc[-1], line_dash="dash", line_color="green",
               annotation_text="ZPHC Echo Stabilized")
fig6.show()

# ─── Export CSV (Optional) ───────────────────────────────────────
df.to_csv("nexus_recursive_swarm_results.csv", index=False)
print("✅ Complete. Full recursion swarm and visualizations are live.")
