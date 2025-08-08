import hashlib
import pandas as pd
import random
from mpmath import mp
import plotly.express as px

# π precision load
mp.dps = 100_010
pi_digits = str(mp.pi)[2:100_010]

def extract_pi_byte(index: int):
    return pi_digits[index:index+8] if index + 8 <= len(pi_digits) else None

def drift_and_sti(byte: str):
    deltas = [abs(int(byte[i+1]) - int(byte[i])) for i in range(7)]
    avg_drift = sum(deltas) / 7
    sti = round(1 - avg_drift / 9, 3)
    return avg_drift, sti

# One recursive memory allocator agent
def allocate_symbolic_pi_memory(seed: str, hops=30, jump=631):
    allocations = []
    H = seed
    for i in range(hops):
        sha = hashlib.sha256(hashlib.sha256((H + str(i)).encode()).digest()).hexdigest()
        index = int(sha[:6], 16) % (len(pi_digits) - 8)
        index = (index + i * jump) % (len(pi_digits) - 8)
        byte = extract_pi_byte(index)
        avg_drift, sti = drift_and_sti(byte)
        allocations.append({
            "agent": seed,
            "hop": i,
            "pi_index": index,
            "byte": byte,
            "avg_drift": avg_drift,
            "sti": sti,
            "zphc": sti >= 0.7
        })
        H = sha
    return allocations

# Deploy swarm of symbolic memory allocators
def run_symbolic_allocators(agent_count=36, hops=25):
    seeds = [f"AGENT-{random.randint(1000,9999)}" for _ in range(agent_count)]
    all_allocs = []
    for seed in seeds:
        all_allocs.extend(allocate_symbolic_pi_memory(seed, hops=hops))
    return pd.DataFrame(all_allocs)

# Run it
df_alloc = run_symbolic_allocators()

# Plot 1: π-RAM field map of trust
fig1 = px.scatter(df_alloc, x="pi_index", y="sti", color="zphc", 
                  title="π Symbolic Allocator Field Map",
                  labels={"sti": "Symbolic Trust Index", "zphc": "ZPHC"},
                  color_discrete_map={True: "green", False: "red"})
fig1.add_hline(y=0.7, line_dash="dash", line_color="black", annotation_text="ZPHC Threshold")
fig1.show()

# ✅ FIX: Convert Interval → string so Plotly doesn't choke
df_alloc["pi_zone"] = pd.cut(df_alloc["pi_index"], bins=20).astype(str)

# Plot 2: Trust density per π zone
zone_stats = df_alloc.groupby("pi_zone", observed=False)["zphc"].mean().reset_index()
fig2 = px.bar(zone_stats, x="pi_zone", y="zphc",
              title="ZPHC Lock Rate by π Memory Zone",
              labels={"zphc": "ZPHC Lock Rate", "pi_zone": "π Memory Block"})
fig2.update_xaxes(tickangle=45)
fig2.show()

# Save symbolic memory allocation
df_alloc.to_csv("nexus_pi_symbolic_allocator.csv", index=False)
print("✅ π symbolic allocator log saved: nexus_pi_symbolic_allocator.csv")
