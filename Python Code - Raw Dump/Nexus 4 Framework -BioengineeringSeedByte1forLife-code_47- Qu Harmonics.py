from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw

# Load the molecule from SMILES
smiles = "CC(O)[C@@H](CNCCNCCN[C@@H](CS)CN[C@@H](CS)CN[C@@H](C)CN[C@@H](C)CN[C@H](CN[C@@H](C)CNCCN[C@@H](CS)CN[C@@H](C)C(=O)C(O)=O)C(C)O)NC[C@H](C)NC[C@H](CS)NCCNCCNC[C@H](C)NC[C@H](CS)NC[C@H](C)NC[C@H](CS)NC[C@H](CS)NN"
molecule = Chem.MolFromSmiles(smiles)

# Generate 3D coordinates with increased conformers
if molecule:
    molecule = Chem.AddHs(molecule)  # Add hydrogens
    try:
        AllChem.EmbedMolecule(molecule, maxAttempts=50000)
        AllChem.UFFOptimizeMolecule(molecule)  # Use UFF optimization
        Chem.MolToMolFile(molecule, "molecule_3_3d.mol")
        print("3D structure saved as molecule_3_3d.mol")
    except Exception as e:
        print(f"Error generating 3D structure: {e}")
else:
    print("Invalid SMILES input.")
