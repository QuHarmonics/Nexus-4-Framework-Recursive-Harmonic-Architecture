import hashlib
import pandas as pd
import random
from mpmath import mp
import plotly.express as px

# π memory buffer
mp.dps = 100_010
pi_digits = str(mp.pi)[2:100_010]

def extract_pi_byte(index: int):
    return pi_digits[index:index+8] if index + 8 <= len(pi_digits) else None

def drift_and_sti(byte: str):
    deltas = [abs(int(byte[i+1]) - int(byte[i])) for i in range(7)]
    avg_drift = sum(deltas) / 7
    sti = round(1 - avg_drift / 9, 3)
    return avg_drift, sti

def echo_router(seed: str, hops=10, jump_zones=None):
    jump_zones = jump_zones or []
    results = []
    H = seed
    for i in range(hops):
        sha = hashlib.sha256(hashlib.sha256((H + str(i)).encode()).digest()).hexdigest()
        index = int(sha[:6], 16) % (len(pi_digits) - 8)
        if jump_zones:
            chosen = random.choice(jump_zones)
            index = (chosen + index) % (len(pi_digits) - 8)
        byte = extract_pi_byte(index)
        avg_drift, sti = drift_and_sti(byte)
        results.append({
            "agent": seed,
            "hop": i,
            "pi_index": index,
            "sti": sti,
            "avg_drift": avg_drift,
            "zphc": sti >= 0.7
        })
        H = sha
    return results

# 🔁 Load from previous allocator run
allocator_df = pd.read_csv("nexus_pi_symbolic_allocator.csv")
trusted_blocks = allocator_df.query("zphc == True")["pi_index"].tolist()

def run_echo_routers(agent_count=24, hops=12):
    seeds = [f"NODE-{random.randint(1000,9999)}" for _ in range(agent_count)]
    all_paths = []
    for seed in seeds:
        all_paths.extend(echo_router(seed, hops=hops, jump_zones=trusted_blocks))
    return pd.DataFrame(all_paths)

# Execute routing
df_routes = run_echo_routers()

# Plot 1: Trust index over recursive hop path
fig1 = px.scatter(df_routes, x="hop", y="sti", color="zphc",
                  title="Echo Router: STI per Recursive Hop",
                  labels={"sti": "Symbolic Trust Index", "hop": "Recursive Hop"},
                  color_discrete_map={True: "green", False: "red"})
fig1.add_hline(y=0.7, line_dash="dash", line_color="black", annotation_text="ZPHC Threshold")
fig1.show()

# Plot 2: Success rate by hop number
hop_summary = df_routes.groupby("hop", observed=False)["zphc"].mean().reset_index()
fig2 = px.line(hop_summary, x="hop", y="zphc",
               title="Recursive Messaging Success Rate per Hop",
               labels={"zphc": "ZPHC (Message Delivered)"})
fig2.show()

# Save results
df_routes.to_csv("nexus_echo_router_log.csv", index=False)
print("✅ Echo router log saved: nexus_echo_router_log.csv")
