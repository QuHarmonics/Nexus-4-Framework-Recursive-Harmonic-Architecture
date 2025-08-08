import numpy as np

def qft(data, H, F, t):
    # Quantum Folding Tool (QFT)
    folded_data = np.sum(data / np.exp(H * F * t))
    return folded_data

def qut(folded_data, theta, zeta):
    # Quantum Unfolding Tool (QUT)
    unfolded_data = np.sum(folded_data * np.cos(theta)) + zeta
    return unfolded_data

# Example usage:
data = np.random.rand(100)  # Electromagnetic field data

# Define parameter ranges for grid search
H_range = np.linspace(0.5, 1.5, 5)
F_range = np.linspace(0.2, 0.8, 5)
t_range = np.linspace(1, 3, 5)
theta_range = np.linspace(np.pi/6, np.pi/3, 5)
zeta_range = np.linspace(0.05, 0.2, 5)

best_params = None
best_compression_ratio = 0
best_reconstruction_error = float('inf')

for H in H_range:
    for F in F_range:
        for t in t_range:
            for theta in theta_range:
                for zeta in zeta_range:
                    folded_data = qft(data, H, F, t)
                    unfolded_data = qut(folded_data, theta, zeta)
                    
                    compression_ratio = len(data) / len(str(folded_data))
                    reconstruction_error = np.mean(np.abs(unfolded_data - data))
                    
                    if compression_ratio > best_compression_ratio or (compression_ratio == best_compression_ratio and reconstruction_error < best_reconstruction_error):
                        best_params = (H, F, t, theta, zeta)
                        best_compression_ratio = compression_ratio
                        best_reconstruction_error = reconstruction_error

print("Best parameters:", best_params)
print("Best compression ratio:", best_compression_ratio)
print("Best reconstruction error:", best_reconstruction_error)