from rdkit import Chem
from rdkit.Chem import Draw

# Create the PSREQ peptide structure (simplified for visualization purposes)
psreq_smiles = "C(C(=O)NCC(C(=O)NCC(C(=O)NCC(C(=O)NCC(C(=O)O)C)C)C)C)C"  # Simplified peptide backbone
psreq_mol = Chem.MolFromSmiles(psreq_smiles)

# Render the molecular diagram
img = Draw.MolToImage(psreq_mol, size=(500, 500))
img.show()
