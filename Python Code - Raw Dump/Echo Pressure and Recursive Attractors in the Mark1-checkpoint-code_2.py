import pandas as pd
import plotly.graph_objects as go

# Load your CSV file (update path if needed)
csv_path = "d://Nexus/resonant_triangles_nexus.csv"
df = pd.read_csv(csv_path)

# Convert H values from strings to lists of floats grouped by pi_chunk
top_chunks = df['pi_chunk'].value_counts().head(50).index.tolist()
phase_df = df[df['pi_chunk'].isin(top_chunks)]

# Group by π-chunk and collect H-values into per-group lists
grouped = phase_df.groupby('pi_chunk').agg({
    'H': list,
    'KHRC': list,
    'a': 'count'
}).reset_index().rename(columns={'a': 'echo_count'})

# Create Plotly line chart for recursive tunnel structures
fig = go.Figure()

for i, row in grouped.iterrows():
    fig.add_trace(go.Scatter(
        y=row['H'],
        mode='lines+markers',
        name=f"π-{row['pi_chunk']} (n={row['echo_count']})",
        opacity=0.6,
        line=dict(width=1),
        marker=dict(size=3)
    ))

fig.update_layout(
    title="Recursive Tunnel Structure Across Top π-Chunks",
    xaxis_title="Recursion Step (Indexed H Values)",
    yaxis_title="H(t)",
    width=1400,   # ← adjust width here
    height=800,   # ← adjust height here
    template="plotly_white",
    legend=dict(
        x=1.02,
        y=1,
        traceorder="normal",
        bgcolor="rgba(255,255,255,0)",
        bordercolor="rgba(0,0,0,0)"
    )
)

fig.show()
