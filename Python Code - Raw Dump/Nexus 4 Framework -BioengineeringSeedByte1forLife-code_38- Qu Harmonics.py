from rdkit import Chem
from rdkit.Chem import AllChem

# Define the peptide sequence
peptide_sequence = "CCACAGGCATGGCCAATAGCA"

# Use a peptide to SMILES translator function
# Define a map of amino acids to their SMILES representation
amino_acid_smiles = {
    'A': 'N[C@@H](C)C(=O)',    # Alanine
    'C': 'N[C@@H](CS)C(=O)',  # Cysteine
    'G': 'NCC(=O)',           # Glycine
    'T': 'N[C@@H](C(O)C)C(=O)', # Threonine
}

# Build the SMILES representation of the peptide
smiles_peptide = ""
for amino_acid in peptide_sequence:
    smiles_peptide += amino_acid_smiles.get(amino_acid, '')

# Add proper connectivity
smiles_peptide = smiles_peptide[:-6]  # Remove the last peptide bond fragment

# Display the result
smiles_peptide
