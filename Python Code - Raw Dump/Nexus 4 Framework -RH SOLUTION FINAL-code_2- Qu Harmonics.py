import numpy as np
import plotly.graph_objects as go

# Predictive Harmonic Framework with Dynamic Lift Seeding
def predict_zeros_with_dynamic_seed(iterations, alpha=1.5, target=0.5, initial_ratio=0.47):
    predictions = [target]
    dynamic_ratios = [initial_ratio]
    calculated_lifts = []  # To store calculated lift values for observation

    for n in range(1, iterations + 1):
        previous = predictions[-1]
        
        # Calculate the .35 ratio adjustment
        ratio = dynamic_ratios[-1] + (target - previous) * (0.035 / (n + 1))
        correction = (target - previous) / (alpha * ratio * (n + 1))
        
        # Calculate the oscillatory component (harmonic framework without lift)
        value = previous * (-1)**n * np.cos(n / np.pi) + correction
        
        # Dynamic lift for stabilization
        if n == 1:
            # Seed the first lift dynamically
            lift = (target * 0.035) / (n + 1)
        else:
            # Normal lift calculation for subsequent iterations
            lift = ratio * ((target - previous) / (n + 1))
        calculated_lifts.append(lift)  # Store lift for analysis
        
        # Combine oscillation with lift for the final prediction
        value += lift
        predictions.append(value)
        dynamic_ratios.append(ratio)  # Update the ratio
    
    return np.array(predictions), dynamic_ratios, calculated_lifts

# Generate predictions
iterations = 100
predicted_zeros, dynamic_ratios, calculated_lifts = predict_zeros_with_dynamic_seed(iterations)

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
    title="Prediction of Zeta Zeros with Dynamic Lift Seeding",
    xaxis_title="Iteration (n)",
    yaxis_title="Predicted Zeros",
    legend=dict(font=dict(size=12)),
    template="plotly_dark",
    xaxis=dict(showgrid=True),
    yaxis=dict(showgrid=True)
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
    template="plotly_dark",
    xaxis=dict(showgrid=True),
    yaxis=dict(showgrid=True)
)
fig2.show()

# Visualization of Lift Component
fig3 = go.Figure()
fig3.add_trace(go.Scatter(
    x=list(range(iterations)),
    y=calculated_lifts,
    mode='lines',
    name='Calculated Lift',
    line=dict(color='orange', width=2)
))
fig3.update_layout(
    title="Isolated Lift Component during Prediction",
    xaxis_title="Iteration (n)",
    yaxis_title="Calculated Lift",
    legend=dict(font=dict(size=12)),
    template="plotly_dark",
    xaxis=dict(showgrid=True),
    yaxis=dict(showgrid=True)
)
fig3.show()
