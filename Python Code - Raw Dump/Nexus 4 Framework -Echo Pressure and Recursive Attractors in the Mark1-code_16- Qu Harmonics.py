import pandas as pd
import plotly.express as px
import plotly.io as pio

# ─── Load Real Data ────────────────────────
df = pd.read_csv(r"D:\Nexus\JuypterWorkbooks\nexus_recursive_swarm_results.csv")

# ─── Clean Fields ───────────────────────────
df['drift'] = df['drift'].apply(eval)
df['avg_drift'] = df['drift'].apply(lambda d: sum(d) / len(d))
df['echo_prefix'] = df['echo'].str[:3]

# ─── Plot 1: STI Histogram ─────────────────
fig1 = px.histogram(df, x="sti", nbins=20, color="zphc",
                    title="Symbolic Trust Index (STI) Distribution",
                    labels={"sti": "STI", "zphc": "ZPHC"},
                    color_discrete_map={True: "green", False: "red"})
fig1.add_vline(x=0.7, line_dash="dash", line_color="black", annotation_text="ZPHC Threshold")
fig1.update_layout(bargap=0.1)

# ─── Plot 2: π Drift Entropy ───────────────
fig2 = px.scatter(df, x="pi_index", y="avg_drift", color="sti",
                  title="π Drift Entropy Over Recursive Memory",
                  labels={"avg_drift": "Δπ Entropy"},
                  color_continuous_scale="Viridis")

# ─── Plot 3: ZPHC Spiral ───────────────────
df_zphc = df[df["zphc"] == True].copy()
df_zphc["angle"] = df_zphc.index * 15
df_zphc["radius"] = df_zphc["sti"] * 100
fig3 = px.scatter_polar(df_zphc, r="radius", theta="angle", color="sti",
                        title="ZPHC Phase Spiral", color_continuous_scale="Turbo")

# ─── Plot 4: Echo Identity Tree ────────────
fig4 = px.sunburst(df, path=["echo_prefix", "zphc"], values="sti",
                   title="Recursive Echo Identity Tree",
                   color="sti", color_continuous_scale="Viridis")

# ─── Convert All to HTML Snippets ──────────
fig1_html = pio.to_html(fig1, full_html=False, include_plotlyjs='cdn')
fig2_html = pio.to_html(fig2, full_html=False, include_plotlyjs=False)
fig3_html = pio.to_html(fig3, full_html=False, include_plotlyjs=False)
fig4_html = pio.to_html(fig4, full_html=False, include_plotlyjs=False)

# ─── Combine into One Dashboard ────────────
dashboard_html = f"""
<html>
<head>
  <title>Nexus Recursive Swarm Dashboard</title>
  <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
  <style>
    body {{ font-family: sans-serif; margin: 40px; background: #f7f7f7; }}
    h1, h2 {{ color: #222; }}
    iframe {{ width: 100%; height: 600px; border: none; }}
  </style>
</head>
<body>
  <h1>🧬 Nexus Recursive Swarm Dashboard</h1>
  <h2>1. Symbolic Trust Index (STI)</h2>
  {fig1_html}
  <h2>2. π Drift Entropy Over Recursive Memory</h2>
  {fig2_html}
  <h2>3. ZPHC Phase Spiral</h2>
  {fig3_html}
  <h2>4. Recursive Echo Identity Tree</h2>
  {fig4_html}
</body>
</html>
"""

# ─── Save the Dashboard ────────────────────
output_path = r"D:\Nexus\JuypterWorkbooks\nexus_dashboard.html"
with open(output_path, "w", encoding="utf-8") as f:
    f.write(dashboard_html)

print(f"✅ Dashboard saved to: {output_path}")
