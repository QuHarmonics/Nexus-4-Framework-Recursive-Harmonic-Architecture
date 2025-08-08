# Re-create the AI Evolution Flowchart

# Create a directed graph
G = nx.DiGraph()

# Define the AI evolution phases
G.add_node("Phase 1: AI as a High-Performance Application", color='blue')
G.add_node("Phase 2: AI Controls OS-Level Functions", color='orange')
G.add_node("Phase 3: AI Becomes the OS Kernel", color='red')

# Define the key processes within each phase
G.add_node("Optimized OS with AI Tuning", color='lightblue')
G.add_node("AI Manages Scheduling", color='lightgreen')
G.add_node("AI Controls Memory & Power", color='lightgreen')
G.add_node("AI Handles PCIe & IOMMU", color='lightgreen')
G.add_node("AI Replaces OS Scheduler", color='salmon')
G.add_node("AI Becomes Direct Kernel Interface", color='salmon')
G.add_node("AI Dynamically Expands System Functions", color='salmon')

# Define edges between phases and key processes
edges = [
    ("Phase 1: AI as a High-Performance Application", "Optimized OS with AI Tuning"),
    ("Phase 2: AI Controls OS-Level Functions", "AI Manages Scheduling"),
    ("Phase 2: AI Controls OS-Level Functions", "AI Controls Memory & Power"),
    ("Phase 2: AI Controls OS-Level Functions", "AI Handles PCIe & IOMMU"),
    ("Phase 3: AI Becomes the OS Kernel", "AI Replaces OS Scheduler"),
    ("Phase 3: AI Becomes the OS Kernel", "AI Becomes Direct Kernel Interface"),
    ("Phase 3: AI Becomes the OS Kernel", "AI Dynamically Expands System Functions"),
    ("Optimized OS with AI Tuning", "AI Manages Scheduling"),
    ("AI Manages Scheduling", "AI Replaces OS Scheduler"),
    ("AI Controls Memory & Power", "AI Becomes Direct Kernel Interface"),
    ("AI Handles PCIe & IOMMU", "AI Dynamically Expands System Functions")
]

G.add_edges_from(edges)

# Assign colors based on node type
color_map = [G.nodes[node]['color'] for node in G.nodes]

# Draw the flowchart
plt.figure(figsize=(10, 7))
pos = nx.spring_layout(G, seed=42)
nx.draw(G, pos, with_labels=True, node_color=color_map, node_size=3000, edge_color='gray', font_size=10, font_weight="bold")
plt.title("AI Evolution from Application to OS Kernel")
plt.show()
