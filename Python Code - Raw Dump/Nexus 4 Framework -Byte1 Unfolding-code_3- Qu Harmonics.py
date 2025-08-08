import numpy as np
import plotly.graph_objects as go

# Physical constants
G = 6.67430e-11  # gravitational constant, m^3 kg^-1 s^-2
M_sun = 1.98847e30  # kg
M_earth = 5.9722e24  # kg
m_p = 1.6726219e-27  # proton mass, kg

# Recursive gravity parameters (tunable)
gamma = 0.35
beta = 1e-10  # in 1/m
N = 100       # recursion depth

# Sum operator S(r) = sum_{n=0}^N (gamma * exp(-beta * r))^n
def S(r, gamma, beta, N):
    x = gamma * np.exp(-beta * r)
    return (1 - x**(N+1)) / (1 - x)

# Recursive gravity force F_RG(r)
def F_recursive(r, m1, m2, gamma, beta, N):
    # Base Newtonian term factor
    return G * m1 * m2 / r**2 * S(r, gamma, beta, N)

# Newtonian force F_N(r)
def F_newton(r, m1, m2):
    return G * m1 * m2 / r**2

# Range of r values (log scale)
r_vals = np.logspace(7, 16, 200)  # from 10^7 m to 10^16 m

# Compute ratios for Earth-Sun system
ratio_ES = S(r_vals, gamma, beta, N)

# Markers at specific benchmarks
r_ES = 1.496e11  # 1 AU
r_PP = 1e-15     # ~proton-proton separation

ratio_ES_point = S(r_ES, gamma, beta, N)
ratio_PP_point = S(r_PP, gamma, beta, N)

# Plot
fig = go.Figure()

# Main ratio curve
fig.add_trace(go.Scatter(
    x=r_vals, y=ratio_ES, mode='lines', name='Ratio F_RG/F_N'
))

# Benchmark markers
fig.add_trace(go.Scatter(
    x=[r_ES], y=[ratio_ES_point], mode='markers+text', name='Earth–Sun',
    text=[f"{ratio_ES_point:.3f}"], textposition='top right'
))
fig.add_trace(go.Scatter(
    x=[r_PP], y=[ratio_PP_point], mode='markers+text', name='Proton–Proton',
    text=[f"{ratio_PP_point:.3f}"], textposition='bottom left'
))

# Plot layout
fig.update_xaxes(type='log', title='Separation r (m)')
fig.update_yaxes(title='F_RG / F_N (dimensionless)')
fig.update_layout(
    title='Recursive Gravity vs. Newtonian Gravity Ratio',
    legend=dict(x=0.7, y=0.9)
)

fig.show()
