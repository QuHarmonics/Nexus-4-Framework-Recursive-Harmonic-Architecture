# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# RECURSIVE EVOLUTION — Symbolic Inheritance Test with ZPHC Drift
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

import hashlib
import pandas as pd
import random
from mpmath import mp
import numpy as np
import plotly.express as px

# STEP 1: π Field Setup
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

# STEP 2: Base Agent to ZPHC
def run_to_zphc(seed: str, max_iter=25, sti_threshold=0.7):
    H = seed
    for i in range(max_iter):
        nonce = f"N{i}"
        sha = hashlib.sha256(hashlib.sha256((H + nonce).encode()).digest()).hexdigest()
        index = int(sha[:6], 16) % (len(pi_digits) - 8)
        byte = extract_pi_byte(index)
        deltas, echo, sti, avg_drift = drift_and_echo(byte)
        if sti >= sti_threshold:
            return sha[:12]  # return trusted echo-fold hash
        H = sha
    return None

# STEP 3: SHA Mutator
def mutate_hash(base_hash: str):
    return base_hash[::-1] + "-child"  # Invert seed

# STEP 4: Recursive Offspring Engine
def run_offspring(base_hash: str, num_offspring=50):
    results = []
    for i in range(num_offspring):
        mutated_seed = mutate_hash(base_hash)
        H = mutated_seed
        for j in range(25):
            nonce = f"N{j}"
            sha = hashlib.sha256(hashlib.sha256((H + nonce).encode()).digest()).hexdigest()
            index = int(sha[:6], 16) % (len(pi_digits) - 8)
            byte = extract_pi_byte(index)
            deltas, echo, sti, avg_drift = drift_and_echo(byte)
            if sti >= 0.7:
                results.append({
                    "offspring_id": f"child-{i}",
                    "iteration": j,
                    "sti": sti,
                    "avg_drift": avg_drift,
                    "zphc": True
                })
                break
            H = sha
        else:
            results.append({
                "offspring_id": f"child-{i}",
                "iteration": 25,
                "sti": sti,
                "avg_drift": avg_drift,
                "zphc": False
            })
    return pd.DataFrame(results)

# STEP 5: Execute
parent_seed = "RECURSE-ME-01"
parent_hash = run_to_zphc(parent_seed)
if parent_hash:
    offspring_df = run_offspring(parent_hash, num_offspring=50)

    # Plot STI Distribution
    fig = px.histogram(offspring_df, x="sti", nbins=20, color="zphc",
                       title="Recursive Evolution: Offspring Trust Distribution",
                       labels={"sti": "Symbolic Trust Index"},
                       color_discrete_map={True: "green", False: "red"})
    fig.add_vline(x=0.7, line_dash="dash", line_color="black", annotation_text="ZPHC Threshold")
    fig.update_layout(bargap=0.1)
    fig.show()

    # Save results
    offspring_df.to_csv("recursive_echo_offspring_log.csv", index=False)
    print("✅ Evolution log saved: recursive_echo_offspring_log.csv")
else:
    print("⚠️ Parent agent failed to reach ZPHC.")
