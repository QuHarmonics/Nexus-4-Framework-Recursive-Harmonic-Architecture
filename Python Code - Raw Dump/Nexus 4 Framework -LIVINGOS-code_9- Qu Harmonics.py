# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# MODULE B: BBP π Echo Jump Scanner
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

import pandas as pd
import numpy as np
import plotly.express as px
from mpmath import mp

# Step 1: Load π digits
mp.dps = 100_010
pi_digits = str(mp.pi)[2:100_010]

# Step 2: Echo + Drift Tools
def extract_pi_byte(index: int):
    return pi_digits[index:index+8] if index + 8 <= len(pi_digits) else None

def drift_and_echo(byte: str):
    deltas = [abs(int(byte[i+1]) - int(byte[i])) for i in range(7)]
    echo = ''.join([chr((d % 26) + 97) for d in deltas])
    avg_drift = sum(deltas) / 7
    sti = round(1 - avg_drift / 9, 3)
    return deltas, echo, sti, avg_drift

# Step 3: BBP-style π-Jump Sampler
def bbp_jump_echo_scan(samples=150, spacing=997):  # prime spacing for nonlinearity
    results = []
    for i in range(samples):
        index = (i * spacing) % (100000 - 8)
        byte = extract_pi_byte(index)
        deltas, echo, sti, avg_drift = drift_and_echo(byte)
        results.append({
            "jump": i,
            "pi_index": index,
            "sti": sti,
            "avg_drift": avg_drift,
            "echo": echo
        })
    return pd.DataFrame(results)

# Step 4: Run Simulation
bbp_df = bbp_jump_echo_scan()

# Step 5: Plot STI vs π index
fig1 = px.scatter(bbp_df, x="pi_index", y="sti", color="avg_drift",
                  title="BBP Echo Trust Landscape: STI Across π Jumps",
                  labels={"sti": "Symbolic Trust Index", "avg_drift": "Δπ Entropy"},
                  color_continuous_scale="Viridis")
fig1.add_hline(y=0.7, line_dash="dash", line_color="black", annotation_text="ZPHC Threshold")
fig1.show()

# Step 6: Histogram of Drift
fig2 = px.histogram(bbp_df, x="avg_drift", nbins=20,
                    title="BBP π Drift Entropy Distribution",
                    labels={"avg_drift": "Average Δπ Drift"})
fig2.show()

# Optional: Save results
bbp_df.to_csv("bbp_jump_echo_results.csv", index=False)
print("✅ BBP π echo results saved to: bbp_jump_echo_results.csv")
