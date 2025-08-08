import matplotlib.pyplot as plt

# Step 1: Define the nucleotide sequence
nucleotide_sequence = "ATGCATGC"

# Step 2: Convert nucleotides to hexadecimal (based on ASCII encoding)
hex_sequence = [hex(ord(nuc))[2:].upper() for nuc in nucleotide_sequence]

# Step 3: Convert hexadecimal back to nucleotide sequence (for validation)
reconstructed_nucleotides = "".join([chr(int(hx, 16)) for hx in hex_sequence])

# Step 4: Generate protein folding angles (mock-up using ASCII value modulus for simplicity)
angles = [int(hx, 16) % 360 for hx in hex_sequence]

# Generate folding diagram (X, Y coordinates)
x, y = [0], [0]
for angle in angles:
    x.append(x[-1] + 1 * (angle % 2))
    y.append(y[-1] + 1 * ((angle + 1) % 2))

# Step 5: Molecular visualization (mock-up of molecular arrangement)
nodes = list(range(len(hex_sequence)))
positions = [(i, int(hx, 16) % 10) for i, hx in enumerate(hex_sequence)]

# Plotting the results
fig, axs = plt.subplots(3, 1, figsize=(10, 15))

# Genetic Alignment
axs[0].text(0.5, 0.5, f"Original Nucleotide Sequence: {nucleotide_sequence}\n"
                      f"Hexadecimal Sequence: {' '.join(hex_sequence)}\n"
                      f"Reconstructed Nucleotide Sequence: {reconstructed_nucleotides}",
             fontsize=12, wrap=True, ha='center', va='center')
axs[0].set_title("Genetic Alignment (Nucleotide to Hexadecimal)")
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
