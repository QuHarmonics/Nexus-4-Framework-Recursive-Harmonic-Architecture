import numpy as np

L = 256                          # lattice length
tau = 1.8                        # curvature threshold
lam = 0.35                       # KRRB gain
steps = 5000

a = np.zeros((L, 3))
b = np.zeros((L, 3))
g = np.random.randn(L, 3) * 0.1  # tiny random vacuum

def lap(arr):                    # 1-D laplacian helper
    return np.roll(arr,1,0) - 2*arr + np.roll(arr,-1,0)

hits = []

for t in range(steps):
    # γ-layer update (wave mixing + small diffusion)
    g += lam * lap(g) + 0.05 * np.random.randn(*g.shape)

    # curvature magnitude per site
    curv = np.linalg.norm(lap(g), axis=1)

    # detect compression couplers (twin-prime analogues)
    triggers = (curv > tau) & (np.concatenate(([False], curv[:-1] <= tau)))

    # collapse: write into β-layer as integer pair (p, p+2)
    idx = np.where(triggers)[0]
    for k in idx:
        p = 3 + 2*len(hits)      # fake prime sequence for toy demo
        b[k] = (p, p+2, 0)
        hits.append((t, k, p))

    # α-layer PID: simple proportional push toward target H=0.35
    H = np.mean(np.linalg.norm(a, axis=1))
    a += 0.05 * (0.35 - H)

print(f"compression events logged: {len(hits)}")
