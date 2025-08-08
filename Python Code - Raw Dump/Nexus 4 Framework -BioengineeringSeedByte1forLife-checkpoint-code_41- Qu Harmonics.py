from rdkit import Chem
from rdkit.Chem import AllChem

# Define the SMILES for the smaller and larger compounds
smaller_smiles = "NC(=O)C1=CC=C(C=C1)O"  # Example SMILES for the smaller compound
larger_smiles = "N[C@@H](CS)C(=O)N[C@@H](CS)C(=O)N[C@@H](C)C(=O)N[C@@H](CS)C(=O)N[C@@H](C)C(=O)NCC(=O)NCC(=O)N[C@@H](CS)C(=O)N[C@@H](C)C(=O)N[C@@H](C(O)C)C(=O)NCC(=O)NCC(=O)N[C@@H](CS)C(=O)N[C@@H](CS)C(=O)N[C@@H](C)C(=O)N[C@@H](C)C(=O)N[C@@H](C(O)C)C(=O)N[C@@H](C)C(=O)NCC(=O)N[C@@H](CS)C(=O)N[C@@H](C)"

# Convert the SMILES strings to RDKit molecule objects
smaller_mol = Chem.MolFromSmiles(smaller_smiles)
larger_mol = Chem.MolFromSmiles(larger_smiles)

# Ensure molecules have hydrogens added
smaller_mol = Chem.AddHs(smaller_mol)
larger_mol = Chem.AddHs(larger_mol)

# Define the functional groups for reaction
amine_smarts = "[NH2]"  # Amine group in smaller compound
carboxylic_acid_smarts = "[C(=O)O]"  # Carboxylic acid group in larger compound

# Find the amine in the smaller compound
amine_match = smaller_mol.GetSubstructMatch(Chem.MolFromSmarts(amine_smarts))

# Find the carboxylic acid in the larger compound
carboxylic_match = larger_mol.GetSubstructMatch(Chem.MolFromSmarts(carboxylic_acid_smarts))

# If matches are found, proceed with the reaction
if amine_match and carboxylic_match:
    rxn = AllChem.ReactionFromSmarts("[NH2:1].[C(=O)O:2]>>[N:1][C:2](=O)")
    reactants = (smaller_mol, larger_mol)
    products = rxn.RunReactants(reactants)

    # Check if the reaction was successful
    if products:
        # Convert the first product to a SMILES string
        combined_mol = products[0][0]
        combined_smiles = Chem.MolToSmiles(combined_mol)
        print("Combined Compound SMILES:", combined_smiles)
    else:
        print("No products were formed from the reaction.")
else:
    print("No functional groups found for bonding.")
