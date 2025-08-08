import numpy as np
import plotly.graph_objects as go

# Re-run simulation to collect M_series
N, T_max = 30, 50
alpha, beta = 0.1, 1.0
eps, eps_crit = 0.01, 0.1
neighbors = [(-1,0),(1,0),(0,-1),(0,1)]

np.random.seed(42)
H = np.random.rand(N, N)
R = np.random.rand(N, N)

def tick(H, R):
    new_R = R.copy()
    for i in range(N):
        for j in range(N):
            for di, dj in neighbors:
                ni, nj = i+di, j+dj
                if 0 <= ni < N and 0 <= nj < N and H[i,j] > H[ni,nj]:
                    d = abs(H[i,j] - H[ni,nj])
                    new_R[ni,nj] -= alpha * np.tanh(beta * d)
    mis_total = 0.0
    for i in range(N):
        for j in range(N):
            E = H[i,j] * new_R[i,j]
            mis, count = 0.0, 0
            for di, dj in neighbors:
                ni, nj = i+di, j+dj
                if 0 <= ni < N and 0 <= nj < N:
                    E_n = H[ni,nj] * new_R[ni,nj]
                    mis += abs(E - E_n)
                    count += 1
            mis_total += mis / max(count,1)
    M = mis_total / (N*N)
    return new_R, M

M_series = []
for t in range(T_max):
    R, M = tick(H, R)
    M_series.append(M)

# Exponential fit setup
A = M_series[-1]
y = np.array(M_series) - A
t = np.arange(len(y))

# Use all but last point for fitting
t_fit = t[:-1]
y_fit = y[:-1]
mask = y_fit > 1e-6
t_fit = t_fit[mask]
y_fit = y_fit[mask]

if len(t_fit) > 1:
    coeffs = np.polyfit(t_fit, np.log(y_fit), 1)
    slope, intercept = coeffs
    tau = -1 / slope
    B = np.exp(intercept)
    M_fit = A + B * np.exp(-t / tau)
    fit_label = f'Fit τ={tau:.2f}'
else:
    tau = np.nan
    M_fit = None
    fit_label = 'Fit unavailable'

# Plot original and fitted if available
fig = go.Figure()
fig.add_trace(go.Scatter(x=t, y=M_series, mode='lines+markers', name='M(t)'))
if M_fit is not None:
    fig.add_trace(go.Scatter(x=t, y=M_fit, mode='lines', name=fit_label))
fig.update_layout(
    title='Global Misalignment with Exponential Fit',
    xaxis_title='Tick',
    yaxis_title='M(t)',
    legend=dict(x=0.7, y=0.95)
)
fig.show()

print(f"Estimated tau: {tau:.2f}")
