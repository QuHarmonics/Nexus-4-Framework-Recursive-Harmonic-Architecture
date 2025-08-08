from rdkit import Chem
from rdkit.Chem import AllChem

# Define the peptide sequence
peptide_sequence = "CCACAGGCATGGCCAATAGCA"

# Define a map of amino acids to their SMILES representation
amino_acid_smiles = {
    'A': 'N[C@@H](C)C(=O)',    # Alanine
    'C': 'N[C@@H](CS)C(=O)',  # Cysteine
    'G': 'NCC(=O)',           # Glycine
    'T': 'N[C@@H](C(O)C)C(=O)', # Threonine
}

# Build the SMILES representation of the peptide
smiles_peptide = ""
for i, amino_acid in enumerate(peptide_sequence):
    smiles_peptide += amino_acid_smiles.get(amino_acid, '')
    # Remove the terminal =O of the current amino acid if it's not the last
    if i < len(peptide_sequence) - 1:
        smiles_peptide = smiles_peptide[:-4]  # Remove the last =O to connect with the next amino acid

# Add N-terminal and C-terminal groups
smiles_peptide = "N" + smiles_peptide + "C(=O)O"

# Display the result
print(smiles_peptide)
