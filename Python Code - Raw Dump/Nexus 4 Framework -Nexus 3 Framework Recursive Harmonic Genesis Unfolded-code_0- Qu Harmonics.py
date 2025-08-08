import numpy as np
import matplotlib.pyplot as plt

# Define the complex function
def f(x):
    return x**4 - 4*x**2 + 2*x + np.sin(10*x)

# Harmonic optimization function
def harmonic_optimize(x_range, steps=100, freq_scale=10):
    x = np.linspace(x_range[0], x_range[1], 1000)
    y = f(x)
    
    # Initialize wave amplitudes (representing "resonance" at each x)
    amplitudes = np.ones_like(x)
    
    # Simulate harmonic collapse
    for _ in range(steps):
        # Compute "frequencies" based on position and function value
        frequencies = freq_scale * (x + y / np.max(np.abs(y)))
        
        # Update amplitudes: lower energy (y) increases amplitude (constructive),
        # higher energy decreases it (destructive interference)
        amplitudes *= np.exp(-y / np.max(np.abs(y)))
        
        # Normalize to prevent divergence
        amplitudes /= np.max(amplitudes) + 1e-10
    
    # Find the "resonant" point (max amplitude = global minimum)
    min_idx = np.argmax(amplitudes)
    return x[min_idx], y[min_idx], x, y, amplitudes

# Run the optimization
x_min, y_min, x_vals, y_vals, amps = harmonic_optimize([-2, 2])

# Plot results
plt.figure(figsize=(10, 6))
plt.plot(x_vals, y_vals, label='f(x)')
plt.plot(x_vals, amps * np.max(y_vals), label='Amplitude (Resonance)')
plt.scatter(x_min, y_min, color='red', label=f'Global Min: ({x_min:.2f}, {y_min:.2f})')
plt.legend()
plt.title('Harmonic Optimization of f(x)')
plt.xlabel('x')
plt.ylabel('f(x) / Amplitude')
plt.grid(True)
plt.show()

print(f"Global minimum found at x = {x_min:.2f}, f(x) = {y_min:.2f}")