import numpy as np
import plotly.graph_objects as go

# Physical constants
G = 6.67430e-11         # m³·kg⁻¹·s⁻²
M_sun = 1.98847e30      # kg
M_earth = 5.9722e24     # kg
m_p = 1.6726219e-27     # kg

# Recursive-gravity hyperparameters
gamma = 0.35            # harmonic attractor
beta = 1e-10            # decay constant (1/m)
N = 100                 # recursion depth

# Recursive sum S(r) = ∑ₙ₌₀ᴺ (γ e^{-βr})ⁿ
def S(r, gamma, beta, N):
    x = gamma * np.exp(-beta * r)
    return (1 - x**(N+1)) / (1 - x)

# Separation range (log-spaced)
r_vals = np.logspace(7, 16, 200)  # 10⁷ m … 10¹⁶ m

# Compute ratio F_RG / F_N = S(r)
ratio = S(r_vals, gamma, beta, N)

# Benchmark points
r_ES = 1.496e11   # Earth–Sun distance (m)
r_PP = 1e-15      # Proton–proton separation (m)
ratio_ES = S(r_ES, gamma, beta, N)
ratio_PP = S(r_PP, gamma, beta, N)

# Plot
fig = go.Figure([
    go.Scatter(x=r_vals, y=ratio, mode='lines', name='F_RG/F_N'),
    go.Scatter(x=[r_ES], y=[ratio_ES], mode='markers+text', name='Earth–Sun',
               text=[f"{ratio_ES:.3f}"], textposition='top right'),
    go.Scatter(x=[r_PP], y=[ratio_PP], mode='markers+text', name='Proton–Proton',
               text=[f"{ratio_PP:.3f}"], textposition='bottom left'),
])
fig.update_xaxes(type='log', title='Separation r (m)')
fig.update_yaxes(title='F_RG / F_N (dimensionless)')
fig.update_layout(title='Recursive vs. Newtonian Gravity Ratio')
fig.show()
