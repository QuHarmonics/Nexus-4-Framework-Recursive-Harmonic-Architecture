from rdkit import Chem
from rdkit.Chem import rdChemReactions, Draw

# Define SMILES for the smaller and larger compounds
smaller_smiles = "Smaller compound SMILES here"
larger_smiles = "Larger compound SMILES here"

# Create molecule objects
smaller_mol = Chem.MolFromSmiles(smaller_smiles)
larger_mol = Chem.MolFromSmiles(larger_smiles)

# Define a potential reaction mechanism
reaction_smarts = "[C:1](=[O:2])-[O:3].[N:4]>>[C:1](=[O:2])[N:4]"
reaction = rdChemReactions.ReactionFromSmarts(reaction_smarts)

# Perform the reaction
products = reaction.RunReactants((smaller_mol, larger_mol))

# Visualize the results
for product_set in products:
    for product in product_set:
        Chem.SanitizeMol(product)
        Draw.MolToImage(product).show()

# Save results to a file
Draw.MolToFile(product, "combined_product.png")

# Display the products' SMILES
for product_set in products:
    for product in product_set:
        print(Chem.MolToSmiles(product))
