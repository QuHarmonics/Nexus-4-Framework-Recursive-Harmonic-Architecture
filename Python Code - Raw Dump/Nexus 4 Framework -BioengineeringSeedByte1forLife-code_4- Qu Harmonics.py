import matplotlib.pyplot as plt
import networkx as nx

# Create a graph for chemical bonding visualization
G = nx.Graph()

# Adding nodes for elements in the peptide
nodes = [
    "C1", "C2", "C3", "N", "H1", "H2", "H3", "O1", "O2",  # Backbone
    "Zn", "Mg1", "Mg2"  # Stabilizing ions
]

# Adding bonds between nodes
edges = [
    ("C1", "C2"), ("C2", "C3"), ("C3", "N"),  # Backbone
    ("C1", "O1"), ("C1", "H1"), ("C2", "H2"), ("C3", "H3"),  # Hydrogens and oxygens
    ("N", "H1"),  # Amino hydrogen
    ("C3", "O2"),  # Carboxyl oxygen
    ("N", "Zn"), ("C1", "Mg1"), ("O2", "Mg2")  # Coordination with ions
]

# Add nodes and edges to the graph
G.add_nodes_from(nodes)
G.add_edges_from(edges)

# Positioning for better visualization
pos = nx.spring_layout(G)

# Plot the chemical bonding diagram
plt.figure(figsize=(10, 7))
nx.draw(
    G, pos, with_labels=True, node_size=3000, node_color="lightblue", font_size=10, font_weight="bold"
)
plt.title("Simplified Chemical Bonding Diagram for PSREQ with Zn and Mg Stabilizers")
plt.show()
