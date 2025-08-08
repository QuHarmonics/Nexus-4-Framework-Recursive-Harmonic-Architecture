import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow

# Set up the figure and axis
fig, ax = plt.subplots(figsize=(10, 6))
ax.axis("off")

# Define lifecycle stages and PSREQ disruptions
stages = [
    ("Viral Entry", "PSREQ peptide binds to glycoproteins, blocking receptor interactions."),
    ("Replication", "Inhibits HSV DNA polymerase, halting genome replication."),
    ("Latency", "Targets conserved mechanisms, preventing latency establishment."),
    ("Reactivation", "Destabilizes epigenetic changes, minimizing recurrence."),
    ("Virion Assembly", "Disrupts capsid assembly and viral particle production."),
]

# Coordinates for flowchart
positions = [
    (0.1, 0.8),  # Viral Entry
    (0.4, 0.8),  # Replication
    (0.7, 0.8),  # Latency
    (0.4, 0.6),  # Reactivation
    (0.4, 0.4),  # Virion Assembly
]

# Draw boxes for stages
for pos, (stage, disruption) in zip(positions, stages):
    ax.text(pos[0], pos[1] + 0.05, stage, fontsize=10, fontweight="bold",
            ha="center", bbox=dict(boxstyle="round,pad=0.5", facecolor="lightblue"))
    ax.text(pos[0], pos[1] - 0.05, disruption, fontsize=8, ha="center", wrap=True)

# Add arrows to connect the stages
arrows = [
    ((0.1, 0.8), (0.4, 0.8)),
    ((0.4, 0.8), (0.7, 0.8)),
    ((0.7, 0.8), (0.4, 0.6)),
    ((0.4, 0.6), (0.4, 0.4)),
]

for start, end in arrows:
    ax.add_patch(FancyArrow(start[0], start[1], end[0] - start[0], end[1] - start[1],
                             width=0.01, color="black", head_width=0.03, head_length=0.02))

# Add a title
ax.set_title("Lifecycle Disruptions Caused by the PSREQ Pathway", fontsize=12, fontweight="bold")

# Display the flowchart
plt.show()
