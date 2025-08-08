# Create a simplified interaction plot between PSREQ, ions, and a viral protein

# Define elements for the plot
elements = ["PSREQ Peptide", "Zn²⁺ (Anchor)", "Mg²⁺ (Buffer)", "Viral Protein"]

# Define connections (interactions)
interactions = [
    ("PSREQ Peptide", "Zn²⁺ (Anchor)"),
    ("PSREQ Peptide", "Mg²⁺ (Buffer)"),
    ("Zn²⁺ (Anchor)", "Viral Protein"),
    ("PSREQ Peptide", "Viral Protein"),
]

# Create a network graph for interactions
G = nx.Graph()
G.add_nodes_from(elements)
G.add_edges_from(interactions)

# Position nodes circularly for better visualization
pos = nx.circular_layout(G)

# Plot the molecular interaction diagram
plt.figure(figsize=(8, 8))
nx.draw(
    G, pos, with_labels=True, node_size=3000, node_color="lightyellow", font_size=10, font_weight="bold"
)
plt.title("Molecular Interaction Plot: PSREQ Motif, Stabilizing Ions, and Viral Protein")
plt.show()
