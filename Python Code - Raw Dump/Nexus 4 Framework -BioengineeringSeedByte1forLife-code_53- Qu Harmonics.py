# Define the hex to nucleotide mapping
hex_to_nucleotide = {
    "41": "A",
    "43": "C",
    "47": "G",
    "54": "T"
}

# Example provided sequence for E. coli (partial for demonstration)
asm_sequence = [
    "41", "43", "47", "43", "43", "54", "47", "47", "43",
    "54", "41", "47", "43", "54", "54", "47", "41", "43"
]

# Convert hex sequence into nucleotide sequence by triplets
nucleotide_triplets = [
    hex_to_nucleotide.get(byte, "?") for byte in asm_sequence
]

# Group nucleotides into codons (triplets)
codons = ["".join(nucleotide_triplets[i:i+3]) for i in range(0, len(nucleotide_triplets), 3)]

codons
