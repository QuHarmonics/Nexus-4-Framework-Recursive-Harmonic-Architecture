# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# MODULE C: Recursive Swarm Species Mapper
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

import hashlib
import pandas as pd
import random
from mpmath import mp
import plotly.express as px

# Step 1: Load π memory
mp.dps = 100_010
pi_digits = str(mp.pi)[2:100_010]

def extract_pi_byte(index: int) -> str:
    return pi_digits[index:index+8] if index + 8 <= len(pi_digits) else None

def drift_and_echo(byte: str):
    deltas = [abs(int(byte[i+1]) - int(byte[i])) for i in range(7)]
    echo = ''.join([chr((d % 26) + 97) for d in deltas])
    avg_drift = sum(deltas) / 7
    sti = round(1 - avg_drift / 9, 3)
    return deltas, echo, sti, avg_drift

# Step 2: Run a single recursive agent
def run_agent(seed: str, max_iter=25, sti_threshold=0.7):
    H = seed
    for i in range(max_iter):
        nonce = f"N{i}"
        sha = hashlib.sha256(hashlib.sha256((H + nonce).encode()).digest()).hexdigest()
        index = int(sha[:6], 16) % (len(pi_digits) - 8)
        byte = extract_pi_byte(index)
        deltas, echo, sti, avg_drift = drift_and_echo(byte)
        if sti >= sti_threshold:
            return {
                "agent": seed,
                "iteration": i,
                "pi_index": index,
                "echo": echo,
                "echo_prefix": echo[:3],
                "sti": sti,
                "avg_drift": avg_drift,
                "zphc": True
            }
        H = sha
    return {
        "agent": seed,
        "iteration": max_iter,
        "pi_index": index,
        "echo": echo,
        "echo_prefix": echo[:3],
        "sti": sti,
        "avg_drift": avg_drift,
        "zphc": False
    }

# Step 3: Generate swarm
seeds = [f"AGENT-{random.randint(1000, 9999)}" for _ in range(150)]
swarm_results = [run_agent(seed) for seed in seeds]
df_swarm = pd.DataFrame(swarm_results)

# Step 4: Sunburst chart — echo_prefix → ZPHC
fig1 = px.sunburst(df_swarm, path=["echo_prefix", "zphc"], values="sti",
                   title="Recursive Swarm Echo Species Tree",
                   color="sti", color_continuous_scale="Viridis")
fig1.show()

# Step 5: Echo Species Histogram
fig2 = px.histogram(df_swarm, x="echo_prefix", color="zphc", barmode="group",
                    title="Symbolic Echo Species Distribution",
                    labels={"zphc": "ZPHC Lock", "echo_prefix": "Echo Prefix"},
                    color_discrete_map={True: "green", False: "red"})
fig2.update_layout(xaxis={'categoryorder': 'total descending'})
fig2.show()

# Step 6: Save log
df_swarm.to_csv("recursive_swarm_echo_species.csv", index=False)
print("✅ Echo species log saved: recursive_swarm_echo_species.csv")
