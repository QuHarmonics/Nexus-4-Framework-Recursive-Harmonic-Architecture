from rdkit import Chem
from rdkit.Chem import Draw

# Load the structure from the given molecular image
# Using the visual representation, manually encode the SMILES as accurately as possible for the peptide shown.
# Based on the image analysis, let's generate its SMILES
image_smiles = "CC(C(=O)O)NC(=O)C(C(=O)NC(=O)C(C(=O)NC(=O)C(C)C)C)C"

# Create a molecule object
molecule = Chem.MolFromSmiles(image_smiles)

# Draw the molecule to confirm accuracy
Draw.MolToImage(molecule)
