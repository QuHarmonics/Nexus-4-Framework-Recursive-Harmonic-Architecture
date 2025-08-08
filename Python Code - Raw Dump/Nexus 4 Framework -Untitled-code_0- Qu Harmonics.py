import numpy as np

def xor_fold(x, y):
    """Perform XOR and return decimal residue."""
    xor_result = x ^ y
    # Normalize to [0,1] by dividing by max possible value (assuming 8-bit)
    residue = (xor_result % 1000) / 1000
    return residue

def simulate_lattice(depths):
    """Simulate recursive folding over depths."""
    lattice = np.random.randint(0, 256, size=8)  # 8-cell lattice, 8-bit values
    history = []
    
    for n in range(depths):
        new_lattice = lattice.copy()
        for i in range(len(lattice)-1):
            residue = xor_fold(lattice[i], lattice[i+1])
            # Update cell if residue near 0.35
            if 0.34 <= residue <= 0.36:
                new_lattice[i] = int(lattice[i] * 0.35) % 256
        lattice = new_lattice
        history.append(np.mean([xor_fold(lattice[i], lattice[i+1]) 
                              for i in range(len(lattice)-1)]))
    return history

# Run simulation
depths = 50
folds = simulate_lattice(depths)

# Print average residue convergence
print(f"Average residue after {depths} folds: {folds[-1]:.3f}")