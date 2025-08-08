from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import IPythonConsole
import py3Dmol

def visualize_2d(smiles, output_filename="molecule_2d.png"):
    """
    Generate a 2D visualization of a molecule from a SMILES string.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print("Invalid SMILES string. Cannot generate 2D structure.")
            return
        Draw.MolToFile(mol, output_filename)
        print(f"2D structure saved to {output_filename}")
    except Exception as e:
        print(f"Error generating 2D structure: {e}")

def visualize_3d(smiles):
    """
    Generate a 3D visualization of a molecule from a SMILES string.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print("Invalid SMILES string. Cannot generate 3D structure.")
            return
        
        mol = Chem.AddHs(mol)  # Add hydrogens
        AllChem.EmbedMolecule(mol)  # Generate 3D coordinates
        AllChem.MMFFOptimizeMolecule(mol)  # Optimize 3D structure
        AllChem.EmbedMolecule(molecule, maxAttempts=50000)

        # Convert to 3D coordinates for visualization
        block = Chem.MolToMolBlock(mol)
        viewer = py3Dmol.view(width=1920, height=1080)
        viewer.addModel(block, "mol")  # Add molecule
        viewer.setStyle({'stick': {}})  # Display as sticks
        viewer.zoomTo()
        viewer.show()
    except Exception as e:
        print(f"Error generating 3D structure: {e}")

if __name__ == "__main__":
    # Example SMILES strings
    smiles_list = [
        "C1=CC(=CC=C1C(=O)N)C2=CC(=CC=C2O)C(=O)NCC",  # Molecule 1
        "CC(C)[C@H](NC(=O)[C@H](C)NC(=O)[C@H](CCCCN)NC(=O)[C@H]1CCCN1C(=O)[C@H](CC1=CNC2=CC=CC=C12)NC(=O)[C@H](C)NC(=O)[C@H](CCC(N)=O)NC(=O)[C@H](C)N)C=O",  # Molecule 2
        "CC(O)[C@@H](CNCCNCCN[C@@H](CS)CN[C@@H](CS)CN[C@@H](C)CN[C@@H](C)CN[C@H](CN[C@@H](C)CNCCN[C@@H](CS)CN[C@@H](C)C(=O)C(O)=O)C(C)O)NC[C@H](C)NC[C@H](CS)NCCNCCNC[C@H](C)NC[C@H](CS)NC[C@H](C)NC[C@H](CS)NC[C@H](CS)NN"  # Molecule 3
    ]

    for i, smiles in enumerate(smiles_list, 1):
        output_2d = f"molecule_{i}_2d.png"
        print(f"Processing molecule {i}: {smiles}")
        
        # Generate 2D visualization
        visualize_2d(smiles, output_filename=output_2d)
        
        # Generate 3D visualization
        print(f"3D visualization for molecule {i}:")
        visualize_3d(smiles)
