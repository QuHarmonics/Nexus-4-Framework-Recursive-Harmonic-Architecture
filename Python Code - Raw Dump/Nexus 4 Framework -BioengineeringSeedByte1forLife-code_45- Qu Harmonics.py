from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

# Input: SMILES representation of the chemical structure
smiles = "[Mg+2].[Mg+2].[Zn+2].[Zn+2].C1=CC(C(=O)O)=CC=C1NCC(C(=O)O)NC(C(=O)O)CCN"

# Create a molecule object
mol = Chem.MolFromSmiles(smiles)

# Check if molecule is valid
if mol is None:
    print("Invalid SMILES string. Please check the input.")
else:
    # Generate molecule formula
    formula = CalcMolFormula(mol)
    print(f"Molecular Formula: {formula}")
    
    # Visualize the molecule
    print("Drawing molecule...")
    Draw.MolToFile(mol, "molecule.png")
    
    # Define SMARTS patterns for functional groups
    functional_groups = {
        "Amine": "[NX3;H2,H1;!$(NC=O)]",  # Primary or secondary amines
        "Carboxylic Acid": "[CX3](=O)[OX1H1]",  # -COOH
        "Sulfhydryl": "[SX2H]",  # -SH
        "Aromatic Ring": "a",  # Any aromatic atoms
    }
    
    # Identify functional groups
    group_matches = {}
    for group, smarts in functional_groups.items():
        pattern = Chem.MolFromSmarts(smarts)
        if pattern:
            matches = mol.GetSubstructMatches(pattern)
            group_matches[group] = matches
            print(f"{group} Matches: {len(matches)}")

    # Generate isomeric pairs
    print("\nIsomeric Pairs:")
    for group1, matches1 in group_matches.items():
        for group2, matches2 in group_matches.items():
            if group1 != group2:  # Pair different functional groups
                print(f"Pairing {group1} with {group2}:")
                for m1 in matches1:
                    for m2 in matches2:
                        print(f"  {group1} at {m1}, {group2} at {m2}")

    # Save visualization of functional group locations
    for group, matches in group_matches.items():
        highlight_atoms = [atom[0] for atom in matches]
        img = Draw.MolToImage(mol, highlightAtoms=highlight_atoms)
        img.save(f"{group}_highlight.png")

    print("Functional group highlights saved.")
