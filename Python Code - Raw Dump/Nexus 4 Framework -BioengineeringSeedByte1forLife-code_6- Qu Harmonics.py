import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# Define the sequence and regions of interest
sequence = "PSREQGAPAPAGDAPAGHGE"
conserved_regions = [slice(0, 5), slice(8, 12)]  # Conserved motifs in sequence
functional_sites = [15, 16, 17]  # Key functional residues

# Plot settings
fig, ax = plt.subplots(figsize=(12, 2))
ax.set_xlim(0, len(sequence))
ax.set_ylim(0, 1)
ax.axis('off')

# Plot sequence as rectangles with functional and conserved highlights
for i, char in enumerate(sequence):
    # Highlight conserved regions
    if any(i in range(region.start, region.stop) for region in conserved_regions):
        ax.add_patch(Rectangle((i, 0.25), 1, 0.5, color="lightgreen"))
    # Highlight functional sites
    if i in functional_sites:
        ax.add_patch(Rectangle((i, 0.25), 1, 0.5, color="lightblue"))
    # Add residue text
    ax.text(i + 0.5, 0.5, char, ha="center", va="center", fontsize=12, color="black")

# Add annotations for clarity
ax.text(0, 0.8, "Conserved Regions", color="green", fontsize=10)
ax.text(len(sequence) - 5, 0.8, "Functional Sites", color="blue", fontsize=10)

plt.title("Sequence Plot: PSREQ Motif with Conserved and Functional Highlights")
plt.show()
