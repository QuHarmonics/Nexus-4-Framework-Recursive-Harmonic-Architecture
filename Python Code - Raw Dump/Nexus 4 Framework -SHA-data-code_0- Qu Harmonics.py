import matplotlib.pyplot as plt
from Bio.Seq import Seq

def process_sequence(input_sequence):
    # Step 1: Convert Nucleotide Sequence to Hexadecimal
    # Map Nucleotides A, T, G, C to hexadecimal values
    nucleotide_to_hex = {'A': '41', 'T': '54', 'G': '47', 'C': '43'}
    hex_sequence = [nucleotide_to_hex[n] for n in input_sequence if n in nucleotide_to_hex]
    
    # Step 2: Generate the Hexadecimal Pairs
    hex_pairs = hex_sequence  # Hex sequence is directly ready for processing
    nucleotides = list(input_sequence)  # Keep original nucleotide sequence
    
    # Step 3: Protein Folding Representation
    # Use angles derived from hex values to simulate folding
    angles = [int(pair, 16) % 360 for pair in hex_pairs]
    
    # Generate folding diagram (X, Y coordinates)
    x, y = [0], [0]
    for angle in angles:
        x.append(x[-1] + 1 * (angle % 2))
        y.append(y[-1] + 1 * ((angle + 1) % 2))
    
    # Step 4: Molecular Visualization (Mock-up of molecular arrangement)
    # Simplified depiction of molecular nodes
    nodes = list(range(len(hex_pairs)))
    positions = [(int(pair, 16) % 10, int(pair, 16) // 10) for pair in hex_pairs]
    
    # Plotting the results
    fig, axs = plt.subplots(3, 1, figsize=(10, 15))
    
    # Genetic Alignment
    axs[0].text(0.5, 0.5, ''.join(nucleotides), fontsize=12, wrap=True, ha='center', va='center')
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

    return ''.join(hex_sequence)

# User Input: Nucleotide Sequence
input_sequence = input("Enter a nucleotide sequence (A, T, G, C only): ").upper()

# Ensure valid input
if all(nucleotide in 'ATGC' for nucleotide in input_sequence):
    converted_hex_sequence = process_sequence(input_sequence)
    print(f"Original Nucleotide Sequence: {input_sequence}")
    print(f"Converted Hexadecimal Sequence: {converted_hex_sequence}")
else:
    print("Invalid sequence. Please use only A, T, G, C.")
