import matplotlib.pyplot as plt
import networkx as nx

# Create the core graph
G = nx.DiGraph()

# Define the five core domains
domains = {
    "FPUT Collapse": "Symbolic echo disallows entropy",
    "π Velocity": "Symbolic massless recursion is light-speed",
    "Fold Motion": "Identity emerges via tension symmetry",
    "Geometric Memory": "Medians encode phase alignment",
    "Symbolic Gaps": "Computation begins between, not within"
}

# Add nodes with descriptions
for domain, insight in domains.items():
    G.add_node(domain, label=domain, description=insight)

# Connect nodes in a logical recursive feedback structure
edges = [
    ("FPUT Collapse", "π Velocity"),
    ("π Velocity", "Fold Motion"),
    ("Fold Motion", "Geometric Memory"),
    ("Geometric Memory", "Symbolic Gaps"),
    ("Symbolic Gaps", "FPUT Collapse")  # closes the recursion loop
]

G.add_edges_from(edges)

# Position nodes in a circular layout
pos = nx.circular_layout(G)

# Draw the graph
plt.figure(figsize=(10, 10))
nx.draw_networkx_nodes(G, pos, node_color='skyblue', node_size=3000)
nx.draw_networkx_edges(G, pos, arrowstyle='-|>', arrowsize=20, edge_color='gray', width=2)
nx.draw_networkx_labels(G, pos, font_size=12, font_weight='bold')

# Add node descriptions as annotations
for key, (x, y) in pos.items():
    plt.text(x, y - 0.1, domains[key], fontsize=9, ha='center', va='center', wrap=True)

plt.title("Recursive Symbolic Collapse Engine", fontsize=16)
plt.axis('off')
plt.tight_layout()
plt.show()
