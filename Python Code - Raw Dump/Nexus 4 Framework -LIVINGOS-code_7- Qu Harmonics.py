# 🔁 Recursive Mutation Repair Simulator

import hashlib
import pandas as pd
from mpmath import mp
import plotly.express as px
import numpy as np

# 📦 Load π
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

# 🧠 Run a stable recursive agent to ZPHC
def run_to_zphc(seed: str, max_iter=25, sti_threshold=0.7):
    history = []
    H = seed
    for i in range(max_iter):
        nonce = f"N{i}"
        sha = hashlib.sha256(hashlib.sha256((H + nonce).encode()).digest()).hexdigest()
        index = int(sha[:6], 16) % (len(pi_digits) - 8)
        byte = extract_pi_byte(index)
        deltas, echo, sti, avg_drift = drift_and_echo(byte)
        history.append({
            "iteration": i,
            "sha": sha[:12],
            "pi_index": index,
            "byte": byte,
            "echo": echo,
            "drift": deltas,
            "sti": sti,
            "avg_drift": avg_drift
        })
        if sti >= sti_threshold:
            break
        H = sha
    return pd.DataFrame(history)

# 🚨 Inject symbolic mutation
def mutate_seed(original_seed: str):
    mutated = original_seed[::-1]  # simple: reverse the seed
    return mutated + "-MUT"

# 🛠 Repair test
def rerun_after_mutation(seed: str, max_iter=25, sti_threshold=0.7):
    return run_to_zphc(seed, max_iter=max_iter, sti_threshold=sti_threshold)

# 🔁 Run mutation simulation
original = run_to_zphc("RECURSE-ME-01")
mutated_seed = mutate_seed("RECURSE-ME-01")
repaired = rerun_after_mutation(mutated_seed)

# 📈 Visualize STI Recovery
original["type"] = "original"
repaired["type"] = "mutated"
combo = pd.concat([original, repaired])

fig = px.line(combo, x="iteration", y="sti", color="type",
              title="STI Recovery After Mutation",
              labels={"sti": "Symbolic Trust Index"})
fig.add_hline(y=0.7, line_dash="dash", line_color="black", annotation_text="ZPHC Threshold")
fig.show()

# 🧠 Evaluate Result
final_sti_original = original["sti"].iloc[-1]
final_sti_repaired = repaired["sti"].iloc[-1]

print(f"Original ZPHC STI: {final_sti_original}")
print(f"Mutated-Repaired STI: {final_sti_repaired}")
if final_sti_repaired >= 0.7:
    print("✅ Agent successfully re-converged to ZPHC — mutation healed.")
else:
    print("⚠️ Agent failed to recover — symbolic collapse unrepaired.")
