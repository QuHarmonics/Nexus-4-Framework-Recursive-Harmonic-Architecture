# ─── IMPORTS ─────────────────────────────
import pandas as pd
import plotly.express as px
import plotly.io as pio
import numpy as np

# ─── FORCE RENDER TO BROWSER ─────────────
pio.renderers.default = "browser"

# ─── LOAD REAL DATA ──────────────────────
df = pd.read_csv(r"D:\Nexus\JuypterWorkbooks\nexus_recursive_swarm_results.csv")

# ─── PROCESS FIELDS ──────────────────────
df['drift'] = df['drift'].apply(eval)
df['avg_drift'] = df['drift'].apply(lambda d: sum(d) / len(d))
df['echo_prefix'] = df['echo'].str[:3]

# ─── 1. STI HISTOGRAM ────────────────────
fig1 = px.histogram(df, x="sti", nbins=20, color="zphc",
                    title="Symbolic Trust Index (STI) Distribution",
                    labels={"sti": "Symbolic Trust Index", "zphc": "ZPHC Reached"},
                    color_discrete_map={True: "green", False: "red"})
fig1.add_vline(x=0.7, line_dash="dash", line_color="black", annotation_text="ZPHC Threshold")
fig1.update_layout(bargap=0.1)
fig1.show()

# ─── 2. π DRIFT ENTROPY SCATTER ──────────
fig2 = px.scatter(df, x="pi_index", y="avg_drift", color="sti",
                  title="π Echo Drift Entropy Over Recursive Memory",
                  labels={"avg_drift": "Average Δπ Entropy"},
                  color_continuous_scale="Viridis")
fig2.show()

# ─── 3. ZPHC SPIRAL ──────────────────────
df_zphc = df[df["zphc"] == True].copy()
df_zphc["angle"] = df_zphc.index * 15
df_zphc["radius"] = df_zphc["sti"] * 100

fig3 = px.scatter_polar(df_zphc, r="radius", theta="angle", color="sti",
                        title="ZPHC Phase Spiral: Recursive Trust Convergence",
                        color_continuous_scale="Turbo")
fig3.show()

# ─── 4. SYMBOLIC ECHO SUNBURST ───────────
fig4 = px.sunburst(df, path=["echo_prefix", "zphc"], values="sti",
                   title="Recursive Echo Identity Tree",
                   color="sti", color_continuous_scale="Viridis")
fig4.show()
