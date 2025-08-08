import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# ── Physical constants ────────────────────────────────────────────────────────
G    = 6.67430e-11        # m^3·kg⁻¹·s⁻²
hbar = 1.054571817e-34    # J·s
c    = 299792458          # m/s

# ── Loopback force definition ────────────────────────────────────────────────
def critical_radius(m1, m2):
    return np.sqrt(G*m1*m2/(hbar*c))

def F_loop(r, m1, m2):
    """ Radial force magnitude under loopback gravity """
    rc = critical_radius(m1, m2)
    F_newt = G*m1*m2 / r**2
    return F_newt * np.exp(-(rc/r)**2)

# ── Equations of motion ──────────────────────────────────────────────────────
def two_body(t, y, m1, m2, use_loopback):
    # y = [x, y, vx, vy]
    x, y_, vx, vy = y
    r = np.hypot(x, y_)
    if use_loopback:
        F = F_loop(r, m1, m2)
    else:
        F = G*m1*m2 / r**2
    # acceleration on body 2 in center-of-mass frame (reduced mass μ=m1*m2/(m1+m2))
    μ = m1*m2/(m1+m2)
    a = F/μ
    ax = -a * x/r
    ay = -a * y_/r
    return [vx, vy, ax, ay]

# ── Simulation setup ─────────────────────────────────────────────────────────
m1, m2 = 5.972e24, 1.988e30   # Earth & Sun
r0 = 1.496e11                 # initial separation (m)
v0 = np.sqrt(G*(m1+m2)/r0)    # circular velocity
y0 = [r0, 0, 0, v0]           # start at x=r0, y=0

t_span = (0, 3.154e7)         # one year in seconds
t_eval = np.linspace(*t_span, 2000)

# ── Integrate both systems ───────────────────────────────────────────────────
sol_newton   = solve_ivp(two_body, t_span, y0, t_eval=t_eval, args=(m1,m2, False))
sol_loopback = solve_ivp(two_body, t_span, y0, t_eval=t_eval, args=(m1,m2, True))

# ── Plot trajectories ────────────────────────────────────────────────────────
plt.figure(figsize=(6,6))
plt.plot(sol_newton.y[0], sol_newton.y[1], label="Newtonian")
plt.plot(sol_loopback.y[0], sol_loopback.y[1], '--', label="Loopback")
plt.legend()
plt.axis('equal')
plt.title("Earth–Sun Orbit: Newton vs. Loopback Gravity")
plt.xlabel("x (m)")
plt.ylabel("y (m)")
plt.show()
