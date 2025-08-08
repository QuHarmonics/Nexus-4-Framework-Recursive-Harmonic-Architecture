# deterministic lattice update
def step(U_t, psi, G, Phi):
    inside = G.constrain(U_t + psi)        # spatial clamp
    return Phi.filter(inside)              # phase/angle filter

# dream loop
def dream(U_seq, psi0, G, Phi, T):
    psi = psi0
    for t in range(T):
        psi = step(U_seq[t], psi, G, Phi)
    return psi

# life genesis: iterate until ||psi_{t+1}-psi_t|| <= eps
def genesis(U_seq, G, Phi, eps=1e-12, max_iter=1_000_000):
    psi = np.zeros(G.shape)
    for _ in range(max_iter):
        new_psi = step(U_seq[0], psi, G, Phi)   # time‑slice ≈ stationary exterior
        if np.linalg.norm(new_psi-psi) <= eps:
            return new_psi
        psi = new_psi
    raise ConvergenceError
