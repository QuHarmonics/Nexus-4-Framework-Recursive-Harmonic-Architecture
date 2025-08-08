from rdkit import Chem
from rdkit.Chem import AllChem

# Load the molecule from SMILES
smiles = "CC(O)[C@@H](CNCCNCCN[C@@H](CS)CN[C@@H](CS)CN[C@@H](C)CN[C@@H](C)CN[C@H](CN[C@@H](C)CNCCN[C@@H](CS)CN[C@@H](C)C(=O)C(O)=O)C(C)O)NC[C@H](C)NC[C@H](CS)NCCNCCNC[C@H](C)NC[C@H](CS)NC[C@H](C)NC[C@H](CS)NC[C@H](CS)NN"
molecule = Chem.MolFromSmiles(smiles)

if molecule:
    molecule = Chem.AddHs(molecule)  # Add explicit hydrogens

    try:
        # Define embedding parameters with additional options
        params = AllChem.ETKDGv3()
        params.maxAttempts = 1000000
        params.randomSeed = -1
        params.boxSizeMult = 5.0  # Increase embedding box size
        params.useRandomCoords = True  # Start with random coordinates
        params.useMacrocycleTorsions = True  # Improve sampling for large molecules
        params.enforceChirality = True

        # Embed the molecule
        status = AllChem.EmbedMolecule(molecule, params)

        if status == 0:  # Success
            # Perform energy minimization
            optimize_status = AllChem.UFFOptimizeMolecule(molecule, maxIters=5000)

            if optimize_status == 0:  # Optimization succeeded
                Chem.MolToMolFile(molecule, "molecule_3_3d.mol")
                print("3D structure saved as molecule_3_3d.mol")
            else:
                print("Optimization failed or did not converge.")
        else:
            print("Embedding failed after maximum attempts. Consider refining the molecule or parameters.")

    except Exception as e:
        print(f"Error during 3D generation: {e}")
else:
    print("Invalid SMILES input.")
