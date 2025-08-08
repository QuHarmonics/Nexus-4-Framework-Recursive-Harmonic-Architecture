import numpy as np
import plotly.graph_objects as go

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
fig1 = go.Figure()
fig1.add_trace(go.Scatter(
    x=list(range(iterations + 1)),
    y=predicted_zeros,
    mode='lines',
    name='Predicted Zeros',
    line=dict(color='blue', width=2)
))
fig1.add_trace(go.Scatter(
    x=[0, iterations],
    y=[0.5, 0.5],
    mode='lines',
    name='Critical Line (Re(s)=0.5)',
    line=dict(color='red', dash='dash')
))
fig1.update_layout(
    title="Prediction of Zeta Zeros with Incremental Lift Adjustment",
    xaxis_title="Iteration (n)",
    yaxis_title="Predicted Zeros",
    legend=dict(font=dict(size=12)),
    template="plotly_white"
)
fig1.show()

# Visualization of Dynamic Ratios
fig2 = go.Figure()
fig2.add_trace(go.Scatter(
    x=list(range(iterations)),
    y=dynamic_ratios[:-1],
    mode='lines',
    name='Dynamic Ratios',
    line=dict(color='green', width=2)
))
fig2.update_layout(
    title="Evolution of Dynamic Ratios during Prediction",
    xaxis_title="Iteration (n)",
    yaxis_title="Dynamic Ratio",
    legend=dict(font=dict(size=12)),
    template="plotly_white"
)
fig2.show()
