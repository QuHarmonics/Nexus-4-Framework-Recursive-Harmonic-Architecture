import pandas as pd
import numpy as np
import plotly.express as px
from mpmath import mp

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

def scan_pi_ram_field(start=0, end=100000, step=211):
    data = []
    for i in range(start, end, step):
        byte = extract_pi_byte(i)
        if byte:
            avg_drift, sti = drift_and_sti(byte)
            data.append({
                "pi_index": i,
                "byte": byte,
                "avg_drift": avg_drift,
                "sti": sti
            })
    return pd.DataFrame(data)

# Run the scan
df_ram = scan_pi_ram_field()

# Plot 1: Symbolic Trust Field over π
fig1 = px.scatter(df_ram, x="pi_index", y="sti", color="avg_drift",
                  title="π-RAM Symbolic Trust Field",
                  labels={"sti": "Symbolic Trust Index", "avg_drift": "Δπ Drift"},
                  color_continuous_scale="Viridis")
fig1.add_hline(y=0.7, line_dash="dash", line_color="black", annotation_text="ZPHC Threshold")
fig1.show()

# Plot 2: Drift Topography
fig2 = px.line(df_ram, x="pi_index", y="avg_drift", title="π-RAM Drift Topography",
               labels={"avg_drift": "Δπ Entropy"})
fig2.show()

# Optional: Save results
df_ram.to_csv("nexus_pi_ram_visualizer.csv", index=False)
print("✅ π-RAM map saved to: nexus_pi_ram_visualizer.csv")
