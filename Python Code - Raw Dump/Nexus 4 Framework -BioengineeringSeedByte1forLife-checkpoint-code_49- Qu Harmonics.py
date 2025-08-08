from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D

# Reload the molecule from the saved file
mol_3d = Chem.MolFromMolFile("molecule_3_3d.mol")

if mol_3d:
    # Generate 2D depiction for visualization
    drawer = rdMolDraw2D.MolDraw2DSVG(500, 500)
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol_3d)
    drawer.FinishDrawing()

    # Save the visualization as an SVG file
    with open("molecule_3_2d.svg", "w") as f:
        f.write(drawer.GetDrawingText())

    print("2D visualization saved as molecule_3_2d.svg")
else:
    print("Failed to load the 3D molecule for visualization.")
