import networkx as nx
import matplotlib.pyplot as plt

# ─── Build the directed graph ────────────────────────────────────
G = nx.DiGraph()

# Define the AI evolution phases
G.add_node("Phase 1:\nAI as High‑Performance\nApplication", color='skyblue')
G.add_node("Phase 2:\nAI Controls\nOS‑Level Functions", color='orange')
G.add_node("Phase 3:\nAI Becomes\nthe OS Kernel", color='salmon')

# Define the key processes within each phase
processes = {
    "Optimized OS with\nAI Tuning":       'lightblue',
    "AI Manages\nScheduling":             'lightgreen',
    "AI Controls\nMemory & Power":        'lightgreen',
    "AI Handles\nPCIe & IOMMU":           'lightgreen',
    "AI Replaces\nOS Scheduler":          'salmon',
    "AI Becomes\nDirect Kernel Interface":'salmon',
    "AI Dynamically\nExpands Functions":  'salmon',
}
for node, col in processes.items():
    G.add_node(node, color=col)

# Define edges
edges = [
    ("Phase 1:\nAI as High‑Performance\nApplication", "Optimized OS with\nAI Tuning"),
    ("Phase 2:\nAI Controls\nOS‑Level Functions",        "AI Manages\nScheduling"),
    ("Phase 2:\nAI Controls\nOS‑Level Functions",        "AI Controls\nMemory & Power"),
    ("Phase 2:\nAI Controls\nOS‑Level Functions",        "AI Handles\nPCIe & IOMMU"),
    ("Phase 3:\nAI Becomes\nthe OS Kernel",              "AI Replaces\nOS Scheduler"),
    ("Phase 3:\nAI Becomes\nthe OS Kernel",              "AI Becomes\nDirect Kernel Interface"),
    ("Phase 3:\nAI Becomes\nthe OS Kernel",              "AI Dynamically\nExpands Functions"),
    ("Optimized OS with\nAI Tuning",                     "AI Manages\nScheduling"),
    ("AI Manages\nScheduling",                           "AI Replaces\nOS Scheduler"),
    ("AI Controls\nMemory & Power",                      "AI Becomes\nDirect Kernel Interface"),
    ("AI Handles\nPCIe & IOMMU",                         "AI Dynamically\nExpands Functions"),
]
G.add_edges_from(edges)

# Extract colors for drawing
node_colors = [G.nodes[n]['color'] for n in G.nodes()]

# ─── Draw ─────────────────────────────────────────────────────────
plt.figure(figsize=(12, 8))
pos = nx.spring_layout(G, seed=42, k=1.0, iterations=100)

nx.draw_networkx_nodes(
    G, pos,
    node_color=node_colors,
    node_size=2500,
    edgecolors='gray',
    linewidths=1.2
)
nx.draw_networkx_edges(
    G, pos,
    arrowstyle='-|>',
    arrowsize=15,
    edge_color='gray'
)
nx.draw_networkx_labels(
    G, pos,
    font_size=10,
    font_weight='bold'
)

plt.title("AI Evolution: Application → OS Kernel", fontsize=14)
plt.axis('off')
plt.tight_layout()
plt.show()
