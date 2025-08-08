import numpy as np
import plotly.graph_objects as go

# Assume M_series is your list of M(t) values
M = np.array(M_series)
t = np.arange(len(M))

# 1. Estimate asymptote
A = M[-1]

# 2. Compute residuals and select positive ones above a tiny epsilon
y = M - A
epsilon_fit = 1e-6
mask = y > epsilon_fit
t_fit = t[mask]
y_fit = y[mask]

# 3. Only fit if we have enough points
if len(t_fit) >= 2:
    # Linearize: log(y) = log(B) - t/τ
    slope, intercept = np.polyfit(t_fit, np.log(y_fit), 1)
    tau = -1 / slope
    B = np.exp(intercept)
    M_fit = A + B * np.exp(-t / tau)
    fit_label = f'Fit τ={tau:.2f}'
else:
    tau = None
    M_fit = None
    fit_label = 'Insufficient data for fit'

# 4. Plot
fig = go.Figure()
fig.add_trace(go.Scatter(x=t, y=M, mode='lines+markers', name='M(t)'))
if M_fit is not None:
    fig.add_trace(go.Scatter(x=t, y=M_fit, mode='lines', name=fit_label))
fig.update_layout(
    title='Global Misalignment with Optional Exponential Fit',
    xaxis_title='Tick',
    yaxis_title='M(t)'
)
fig.show()

# 5. Report
if tau is not None:
    print(f"Estimated relaxation time τ = {tau:.2f}")
else:
    print("Not enough positive residuals to fit an exponential.")
