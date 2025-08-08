import numpy as np
import matplotlib.pyplot as plt

# Predictive Harmonic Framework with Internal Ratio Calculation
def predict_zeros_with_ratio(iterations, alpha=1.5, target=0.5):
    predictions = [target]
    
    for n in range(1, iterations + 1):
        previous = predictions[-1]
        correction = (target - previous) / (alpha * (n + 1))
        value = previous * (-1)**n * np.cos(n / np.pi) + correction
        
        # Calculate ratio internally
        ratio = 0  # Default ratio for the first step
        if n > 1:  # Only calculate ratio after the first iteration
            ratio = abs(value - predictions[-1]) / abs(predictions[-1]) if abs(predictions[-1]) > 0 else 0
        
        # Use the ratio for dynamic adjustment (if needed)
        # Modify `value` based on internal ratio if this affects your process
        #value=value * ratio
        predictions.append(value)
    
    return np.array(predictions)

# Generate predictions
iterations = 1000  # Increase iterations
predicted_zeros = predict_zeros_with_ratio(iterations)

# Visualization of Predicted Zeros
plt.figure(figsize=(14, 8))
plt.plot(range(iterations + 1), predicted_zeros, label="Predicted Zeros", color="blue", lw=2)
plt.axhline(0.5, color="red", linestyle="--", label="Critical Line (Re(s)=0.5)")
plt.xlabel("Iteration (n)", fontsize=14)
plt.ylabel("Predicted Zeros", fontsize=14)
plt.title("Prediction of Zeta Zeros using Harmonic Framework", fontsize=16)
plt.legend(fontsize=12)
plt.grid()
plt.show()
