# ─── Imports ────────────────────────────────────────────────────
import pandas as pd
import plotly.express as px
import numpy as np
import plotly.io as pio

# ─── Force Plotly to render inline in notebook ──────────────────
pio.renderers.default = "notebook"  # Use "notebook" or "inline"

# ─── Load Real Data ─────────────────────────────────────────────
df = pd.read_csv(r"D:\Nexus\JuypterWorkbooks\nexus_recursive_swarm_results.csv")

# ─── Process Fields ─────────────────────────────────────────────
df['drift'] = df['drift'].apply(eval)
df['avg_drift'] = df['drift'].apply(lambda d: sum(d) / len(d))
df['echo_prefix'] = df['echo'].str[:3]

# ────────────────────────────────────────────────────────────────
# 1. STI Histogram
# ────────────────────────────────────────────────────────────────
fig1 = px.histogram(df, x="sti", nbins=20, color="zphc",
                    title="Symbolic Trust Index (STI) Distribution",
                    labels={"sti": "STI", "zphc": "ZPHC"},
                    color_discrete_map={True: "green", False: "red"})
fig1.add_vline(x=0.7, line_dash="dash", line_color="black", annotation_text="ZPHC Threshold")
fig1.update_layout(bargap=0.1)
fig1.show()

# ────────────────────────────────────────────────────────────────
# 2. π Drift Entropy Map
# ────────────────────────────────────────────────────────────────
fig2 = px.scatter(df, x="pi_index", y="avg_drift", color="sti",
                  title="π Drift Entropy Over Recursive Memory",
                  labels={"avg_drift": "Average Δπ Drift"},
                  color_continuous_scale="Viridis")
fig2.show()

# ────────────────────────────────────────────────────────────────
# 3. ZPHC Spiral Plot
# ────────────────────────────────────────────────────────────────
df_zphc = df[df["zphc"] == True].copy()
df_zphc["angle"] = df_zphc.index * 15
df_zphc["radius"] = df_zphc["sti"] * 100

fig3 = px.scatter_polar(df_zphc, r="radius", theta="angle", color="sti",
                        title="ZPHC Phase Spiral: Recursive Trust Convergence",
                        color_continuous_scale="Turbo")
fig3.show()

# ────────────────────────────────────────────────────────────────
# 4. Recursive Echo Identity Tree
# ────────────────────────────────────────────────────────────────
fig4 = px.sunburst(df, path=["echo_prefix", "zphc"], values="sti",
                   title="Recursive Echo Identity Tree",
                   color="sti", color_continuous_scale="Viridis")
fig4.show()
