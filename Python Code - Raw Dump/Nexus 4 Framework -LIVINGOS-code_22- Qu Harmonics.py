import plotly.graph_objects as go

# Define nodes
nodes = ["Δ¹\n(Asymmetry)", "Δ²\n(Trust Loop)", "Δ³\n(Recursive Self)", "Δ⁴\n(Foldback)"]

# X-Y coordinates for layout
positions = {
    "Δ¹\n(Asymmetry)": (0.1, 0.5),
    "Δ²\n(Trust Loop)": (0.4, 0.8),
    "Δ³\n(Recursive Self)": (0.7, 0.5),
    "Δ⁴\n(Foldback)": (0.4, 0.2),
}

# Edges (arrows)
edges = [
    ("Δ¹\n(Asymmetry)", "Δ²\n(Trust Loop)"),
    ("Δ²\n(Trust Loop)", "Δ³\n(Recursive Self)"),
    ("Δ³\n(Recursive Self)", "Δ⁴\n(Foldback)"),
    ("Δ⁴\n(Foldback)", "Δ¹\n(Asymmetry)"),  # Loop back (dream return)
]

# Create figure
fig = go.Figure()

# Add nodes
for name, (x, y) in positions.items():
    fig.add_trace(go.Scatter(
        x=[x], y=[y],
        mode="markers+text",
        marker=dict(size=60, color='lightblue', line=dict(width=2, color='darkblue')),
        text=[name],
        textposition="middle center",
        hoverinfo="text",
        showlegend=False
    ))

# Add edges as arrows
for start, end in edges:
    x0, y0 = positions[start]
    x1, y1 = positions[end]
    fig.add_annotation(
        ax=x0, ay=y0, axref='x', ayref='y',
        x=x1, y=y1, xref='x', yref='y',
        showarrow=True,
        arrowhead=3,
        arrowsize=1,
        arrowwidth=2,
        arrowcolor="black"
    )

# Layout settings
fig.update_layout(
    title="Δ‑Shape Identity Transition Flow",
    xaxis=dict(visible=False),
    yaxis=dict(visible=False),
    plot_bgcolor="white",
    margin=dict(l=0, r=0, t=50, b=0)
)

fig.show()
