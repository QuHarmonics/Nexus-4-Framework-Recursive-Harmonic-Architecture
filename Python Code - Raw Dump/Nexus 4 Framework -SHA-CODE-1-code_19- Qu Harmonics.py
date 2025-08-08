import numpy as np
import matplotlib.pyplot as plt

# Predictive Harmonic Framework with Incremental Quantum Lift Adjustment
def predict_zeros_with_ratio(iterations, alpha=1.5, target=0.5, initial_ratio=0.47):
    predictions = [target]
    dynamic_ratios = [initial_ratio]
    
    for n in range(1, iterations + 1):
        previous = predictions[-1]
        
        # Dynamically adjust the ratio based on previous changes
        ratio = dynamic_ratios[-1] + (target - previous) * (0.035 / (n + 1))
        correction = (target - previous) / (alpha * ratio * (n + 1))
        
        # Incremental quantum lift adjustment integrated with the iteration
        value = previous * (-1)**n * np.cos(n / np.pi) + correction + ratio * ((target - previous) / (n + 1))
        
        predictions.append(value)
        dynamic_ratios.append(ratio)  # Update the ratio
    
    return np.array(predictions), dynamic_ratios


# Generate predictions with dynamic ratio adjustment
iterations = 100
predicted_zeros, dynamic_ratios = predict_zeros_with_ratio(iterations)

# Visualization of Predicted Zeros
plt.figure(figsize=(14, 8))
plt.plot(range(iterations + 1), predicted_zeros, label="Predicted Zeros", color="blue", lw=2)
plt.axhline(0.5, color="red", linestyle="--", label="Critical Line (Re(s)=0.5)")
plt.xlabel("Iteration (n)", fontsize=14)
plt.ylabel("Predicted Zeros", fontsize=14)
plt.title("Prediction of Zeta Zeros with Incremental Lift Adjustment", fontsize=16)
plt.legend(fontsize=12)
plt.grid()
plt.show()

# Visualization of Dynamic Ratios
plt.figure(figsize=(14, 8))
plt.plot(range(iterations), dynamic_ratios[:-1], label="Dynamic Ratios", color="green", lw=2)
plt.xlabel("Iteration (n)", fontsize=14)
plt.ylabel("Dynamic Ratio", fontsize=14)
plt.title("Evolution of Dynamic Ratios during Prediction", fontsize=16)
plt.legend(fontsize=12)
plt.grid()
plt.show()
