import hashlib
import numpy as np
import matplotlib.pyplot as plt

def sha256_hash(value):
    """Returns SHA-256 hash of the given input value."""
    return hashlib.sha256(value.encode()).hexdigest()

def calculate_entropy(values):
    """Calculates entropy of a given list of values."""
    unique, counts = np.unique(values, return_counts=True)
    probabilities = counts / len(values)
    return -np.sum(probabilities * np.log2(probabilities))

def residue_class_distribution(modulus=64):
    """Tracks residue classes of SHA-256 hashes modulo 'modulus'."""
    residues = []
    for i in range(256):  # Let's check all possible byte values (0-255)
        hash_val = sha256_hash(str(i))
        last_byte = int(hash_val[-2:], 16)  # Extract the last byte
        residues.append(last_byte % modulus)
    return residues

def plot_entropy_vs_modulus():
    """Plots the entropy of SHA-256 hashes as modulus changes from 2 to 256."""
    entropies = []
    for modulus in range(2, 257):  # Check modulus from 2 to 256
        residues = residue_class_distribution(modulus)
        entropies.append(calculate_entropy(residues))
    
    plt.plot(range(2, 257), entropies)
    plt.xlabel('Modulus')
    plt.ylabel('Entropy')
    plt.title('Entropy of SHA-256 Hashes for Different Modulus Values')
    plt.show()

def plot_residue_class_distribution(modulus=64):
    """Plots the residue class distribution for SHA-256 hashes modulo 'modulus'."""
    residues = residue_class_distribution(modulus)
    plt.hist(residues, bins=range(modulus), density=True, alpha=0.75, color='blue')
    plt.xlabel('Residue Classes')
    plt.ylabel('Frequency')
    plt.title(f'Residue Class Distribution for Modulus {modulus}')
    plt.show()

# Run experiments
plot_entropy_vs_modulus()  # Plot entropy vs modulus from 2 to 256
plot_residue_class_distribution(modulus=64)  # Plot residue distribution mod 64
