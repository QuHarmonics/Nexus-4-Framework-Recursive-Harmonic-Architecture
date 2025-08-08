from Bio.SeqUtils.ProtParam import ProteinAnalysis

# Define the peptide sequence
peptide_sequence = "CCACAGGCATGGCCAATAGCA"

# Analyze the sequence
analysis = ProteinAnalysis(peptide_sequence)

# Calculate molecular formula components
amino_acids_count = analysis.count_amino_acids()
total_c = sum([amino_acids_count[aa] * {"A": 3, "C": 3, "G": 2, "T": 4, "S": 3, "N": 4, "E": 5, "Q": 5}.get(aa, 0) for aa in amino_acids_count])
total_h = sum([amino_acids_count[aa] * {"A": 5, "C": 7, "G": 4, "T": 5, "S": 5, "N": 6, "E": 7, "Q": 8}.get(aa, 0) for aa in amino_acids_count])
total_o = sum([amino_acids_count[aa] * {"A": 1, "C": 2, "G": 1, "T": 2, "S": 2, "N": 2, "E": 4, "Q": 3}.get(aa, 0) for aa in amino_acids_count])
total_n = sum([amino_acids_count[aa] * {"A": 1, "C": 1, "G": 1, "T": 1, "S": 1, "N": 2, "E": 1, "Q": 1}.get(aa, 0) for aa in amino_acids_count])
total_s = sum([amino_acids_count[aa] * {"C": 1}.get(aa, 0) for aa in amino_acids_count])

# Adjust for peptide bonds (H2O loss per bond)
num_peptide_bonds = len(peptide_sequence) - 1
total_h -= num_peptide_bonds * 2
total_o -= num_peptide_bonds

# Generate the chemical formula
chemical_formula = f"C{total_c}H{total_h}O{total_o}N{total_n}S{total_s}"
chemical_formula
