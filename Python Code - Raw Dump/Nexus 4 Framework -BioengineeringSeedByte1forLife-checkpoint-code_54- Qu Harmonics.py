# Hex to Nucleotide Mapping
hex_to_nucleotide = {
    "41": "A",
    "43": "C",
    "47": "G",
    "54": "T"
}

# Genetic Code Mapping (Codon to Amino Acid)
codon_to_amino_acid = {
    "ATA": "I", "ATC": "I", "ATT": "I", "ATG": "M",
    "ACA": "T", "ACC": "T", "ACG": "T", "ACT": "T",
    "AAC": "N", "AAT": "N", "AAA": "K", "AAG": "K",
    "AGC": "S", "AGT": "S", "AGA": "R", "AGG": "R",
    "CTA": "L", "CTC": "L", "CTG": "L", "CTT": "L",
    "CCA": "P", "CCC": "P", "CCG": "P", "CCT": "P",
    "CAC": "H", "CAT": "H", "CAA": "Q", "CAG": "Q",
    "CGA": "R", "CGC": "R", "CGG": "R", "CGT": "R",
    "GTA": "V", "GTC": "V", "GTG": "V", "GTT": "V",
    "GCA": "A", "GCC": "A", "GCG": "A", "GCT": "A",
    "GAC": "D", "GAT": "D", "GAA": "E", "GAG": "E",
    "GGA": "G", "GGC": "G", "GGG": "G", "GGT": "G",
    "TCA": "S", "TCC": "S", "TCG": "S", "TCT": "S",
    "TTC": "F", "TTT": "F", "TTA": "L", "TTG": "L",
    "TAC": "Y", "TAT": "Y", "TAA": "*", "TAG": "*",
    "TGC": "C", "TGT": "C", "TGA": "*", "TGG": "W"
}

# Function to process ASM and convert to amino acids
def asm_to_amino_acids(asm_sequence):
    # Convert ASM to nucleotides
    nucleotide_sequence = [hex_to_nucleotide.get(byte, "?") for byte in asm_sequence]

    # Group nucleotides into codons (triplets)
    codons = ["".join(nucleotide_sequence[i:i+3]) for i in range(0, len(nucleotide_sequence), 3)]

    # Translate codons into amino acids
    amino_acid_sequence = [codon_to_amino_acid.get(codon, "?") for codon in codons]

    return codons, "".join(amino_acid_sequence)

# Example ASM sequence (from provided disassembly)
asm_sequence = [
    "41", "54", "47", "41", "43", "43", "41", "54", "47", "41", "54", "54", 
    "41", "43", "47", "43", "43", "41", "41", "47", "43", "54", "54", "47", 
    "41", "54", "54", "43", "54", "54", "54", "54", "41", "43", "41", "47", 
    "43", "54", "54", "47", "54", "43", "43", "41", "47", "47", "47", "47", 
    "47", "41", "54", "43", "43", "41", "54", "47", "54", "41", "41"
]

# Generate codons and amino acid sequence
codons, amino_acids = asm_to_amino_acids(asm_sequence)

# Output results
print("Codons:", codons)
print("Amino Acid Sequence:", amino_acids)
