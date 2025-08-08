import hashlib
import pandas as pd
import random
from mpmath import mp
import plotly.express as px

# Setup π digits
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

def run_to_zphc(seed: str, max_iter=25, sti_threshold=0.7):
    H = seed
    for i in range(max_iter):
        nonce = f"N{i}"
        sha = hashlib.sha256(hashlib.sha256((H + nonce).encode()).digest()).hexdigest()
        index = int(sha[:6], 16) % (len(pi_digits) - 8)
        byte = extract_pi_byte(index)
        deltas, echo, sti, avg_drift = drift_and_echo(byte)
        if sti >= sti_threshold:
            return sha[:12]
        H = sha
    return None

def mutate_and_run(base_hash: str, lineage: str, generation: int):
    mutated = base_hash[::-1] + f"-G{generation}"
    H = mutated
    for j in range(25):
        nonce = f"N{j}"
        sha = hashlib.sha256(hashlib.sha256((H + nonce).encode()).digest()).hexdigest()
        index = int(sha[:6], 16) % (len(pi_digits) - 8)
        byte = extract_pi_byte(index)
        deltas, echo, sti, avg_drift = drift_and_echo(byte)
        if sti >= 0.7:
            return {
                "generation": generation,
                "lineage": lineage,
                "iteration": j,
                "sti": sti,
                "avg_drift": avg_drift,
                "zphc": True,
                "echo_prefix": echo[:3]
            }
        H = sha
    return {
        "generation": generation,
        "lineage": lineage,
        "iteration": 25,
        "sti": sti,
        "avg_drift": avg_drift,
        "zphc": False,
        "echo_prefix": echo[:3]
    }

def simulate_generations(generations=6, offspring_per_gen=12):
    parent_seed = "RECURSE-ME-01"
    base_hash = run_to_zphc(parent_seed)
    all_results = []
    for gen in range(1, generations + 1):
        for k in range(offspring_per_gen):
            lineage_id = f"L{gen:02d}-O{k:02d}"
            result = mutate_and_run(base_hash, lineage=lineage_id, generation=gen)
            all_results.append(result)
    return pd.DataFrame(all_results)

# Run it
df_genome = simulate_generations()

# Plot survival by generation
fig = px.histogram(df_genome, x="generation", color="zphc", barmode="group",
                   title="Recursive Echo Genome: ZPHC Survival by Generation",
                   labels={"generation": "Generation", "zphc": "ZPHC Lock"},
                   color_discrete_map={True: "green", False: "red"})
fig.show()

# Save to CSV
df_genome.to_csv("recursive_echo_genome_log.csv", index=False)
print("✅ Recursive genome evolution log saved: recursive_echo_genome_log.csv")
