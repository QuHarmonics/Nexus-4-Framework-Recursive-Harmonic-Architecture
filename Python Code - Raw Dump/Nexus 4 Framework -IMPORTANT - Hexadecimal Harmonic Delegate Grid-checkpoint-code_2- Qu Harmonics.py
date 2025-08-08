import pandas as pd
import plotly.graph_objects as go

# Define data
data = [
    {"Sequence": [0, 0, 0, 0], "Triplet": [0, 0, 0]},
    {"Sequence": [1, 0, 0, 0], "Triplet": [1, 0, 0]},
    {"Sequence": [1, 0, 0, 1], "Triplet": [1, 0, 1]},
    {"Sequence": [1, 0, 1, 0], "Triplet": [1, 1, 1]},
    {"Sequence": [1, 0, 1, 1], "Triplet": [1, 1, 0]},
    {"Sequence": [1, 1, 0, 0], "Triplet": [0, 1, 0]},
    {"Sequence": [1, 1, 0, 1], "Triplet": [0, 1, 1]},
    {"Sequence": [1, 1, 1, 0], "Triplet": [0, 0, 1]},
    {"Sequence": [1, 1, 1, 1], "Triplet": [0, 0, 0]},
]

# Prepare DataFrame
df = pd.DataFrame(data)
df['Triplet Dec'] = df['Triplet'].apply(lambda x: int("".join(map(str, x)), 2))
df['Sequence'] = df['Sequence'].apply(lambda x: "".join(map(str, x)))

# Plot
fig = go.Figure()

fig.add_trace(go.Scatter(
    x=[i for i in range(len(df))],
    y=df["Triplet Dec"],
    mode='markers+text',
    text=df["Sequence"],
    textposition="top center",
    marker=dict(size=12, color=df["Triplet Dec"], colorscale='Plasma'),
))

fig.update_layout(
    title="🔁 Valid Echo Fixpoints (Symbolic Collapse Triplets)",
    xaxis_title="Fixpoint Index",
    yaxis_title="Triplet Decimal Value",
    height=500,
    width=800,
    margin=dict(l=50, r=50, t=80, b=50)
)

fig.show()
