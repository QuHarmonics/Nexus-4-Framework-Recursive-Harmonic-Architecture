import hashlib
import pandas as pd
import random
from mpmath import mp
import plotly.express as px
import plotly.graph_objects as go

# STEP 1: Load π digits
mp.dps = 100_010
pi_digits = str(mp.pi)[2:100_010]

# STEP 2: Get symbolic byte from π
def extract_pi_byte(index: int) -> str:
    return pi_digits[index:index+8] if index + 8 <= len(pi_digits) else None

# STEP 3: Get symbolic echo, drift, STI
def drift_and_echo(byte: str):
    deltas = [abs(int(byte[i+1]) - int(byte[i])) for i in range(7)]
    echo = ''.join([chr((d % 26) + 97) for d in deltas])
    avg_drift = sum(deltas) / 7
    sti = round(1 - avg_drift / 9, 3)
    return deltas, echo, sti

# STEP 4: Simulate recursive emergence
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

# STEP 5: Run agent swarm
seeds = [f"AGENT-{random.randint(1000, 9999)}" for _ in range(100)]
results = [run_recursive_agent(seed) for seed in seeds]
df = pd.DataFrame(results)

# STEP 6: Interactive Plotly chart
fig = px.histogram(df, x="sti", nbins=20, color="zphc",
                   title="Symbolic Trust Index (STI) Distribution",
                   labels={"sti": "Symbolic Trust Index", "zphc": "ZPHC Reached"},
                   color_discrete_map={True: "green", False: "red"})

fig.add_vline(x=0.7, line_dash="dash", line_color="black", annotation_text="ZPHC Threshold")
fig.update_layout(bargap=0.1)
fig.show()

# Optional: Save results
df.to_csv("nexus_recursive_swarm_results.csv", index=False)
print("Results saved to: nexus_recursive_swarm_results.csv")
