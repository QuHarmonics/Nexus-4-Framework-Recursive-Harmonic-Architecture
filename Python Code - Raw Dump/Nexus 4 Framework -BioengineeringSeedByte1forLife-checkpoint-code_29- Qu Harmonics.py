import matplotlib.pyplot as plt
from Bio.Seq import Seq

# Step 1: Genetic Alignment (Hex to nucleotide sequence approximation)
hex_sequence = """
33 33 20 43 30 0A 46 37 20 45 32 0A 43 37 20 34 35 20 46 43 20 
30 30 20 30 30 20 30 30 20 30 30 0A 45 45 0A 45 42 20 46 45 0A 
33 31 20 44 42 0A 38 42 20 43 33 0A 32 30 32 30 0A 37 35 20 46 36 0A
"""
# Remove whitespaces and convert to list of hex pairs
hex_pairs = hex_sequence.replace("\n", "").split()
nucleotides = []

# Map Hex to Nucleotides A, T, G, C based on hex value mod 4
hex_to_nucleotide = ["A", "T", "G", "C"]
for pair in hex_pairs:
    value = int(pair, 16)
    nucleotides.append(hex_to_nucleotide[value % 4])

# Convert to string for sequence representation
nucleotide_sequence = "".join(nucleotides)

# Step 2: Protein Folding Representation
# Use angles derived from hex values to simulate folding
angles = [int(pair, 16) % 360 for pair in hex_pairs]

# Generate folding diagram (X, Y coordinates)
x, y = [0], [0]
for angle in angles:
    x.append(x[-1] + 1 * (angle % 2))
    y.append(y[-1] + 1 * ((angle + 1) % 2))

# Step 3: Molecular Visualization (Mock-up of molecular arrangement)
# Simplified depiction of molecular nodes
nodes = list(range(len(hex_pairs)))
positions = [(int(pair, 16) % 10, int(pair, 16) // 10) for pair in hex_pairs]

# Plotting the results
fig, axs = plt.subplots(3, 1, figsize=(10, 15))

# Genetic Alignment
axs[0].text(0.5, 0.5, nucleotide_sequence, fontsize=12, wrap=True, ha='center', va='center')
axs[0].set_title("Genetic Alignment (Hex to Nucleotides)")
axs[0].axis("off")

# Protein Folding Representation
axs[1].plot(x, y, marker='o')
axs[1].set_title("Protein Folding Representation")
axs[1].axis("equal")

# Molecular Visualization
for pos, label in zip(positions, nodes):
    axs[2].scatter(pos[0], pos[1], s=100)
    axs[2].text(pos[0], pos[1], str(label), fontsize=8, ha='center')
axs[2].set_title("Molecular Visualization")
axs[2].axis("equal")

plt.tight_layout()
plt.show()

# Displaying the nucleotide sequence separately as well
nucleotide_sequence
