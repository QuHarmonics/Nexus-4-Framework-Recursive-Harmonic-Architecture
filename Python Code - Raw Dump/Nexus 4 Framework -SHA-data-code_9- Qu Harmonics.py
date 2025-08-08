import numpy as np

# -------- parameters ----------
n_var   = 3
block   = 3 * n_var
tau     = 1.5            # curvature threshold
Kp, Ki, Kd = 0.35, 0.02, 0.10   # simple gains
steps   = 500

# -------- curvature field ----
phi   = np.zeros(block)
i_acc = 0.0              # integral term
d_prev= 0.0

# clause packets (3×block array)
packets = np.zeros((3, block))
# c1:  x1 , ~x2,  x3
packets[0, [0, 2, 6]] += 1
# c2: ~x1,  x2,  x3
packets[1, [1, 3, 6]] += 1
# c3:  x1,  x2, ~x3
packets[2, [0, 3, 7]] += 1
# balance each packet
for p in packets:
    p -= p.mean()        # zero-mean injection

# inject all clauses once at t=0
phi += packets.sum(axis=0)

truth = np.array([0, 0, 0])   # start with all False

def lap(v):                   # discrete Laplacian on ring
    return np.roll(v,1)-2*v+np.roll(v,-1)

for t in range(steps):
    err = lap(phi)
    p_term =  Kp * err
    i_acc += Ki * err
    d_term =  Kd * (err - d_prev)
    control = p_term + i_acc + d_term
    d_prev  = err.copy()

    # apply control
    phi -= control

    # Δ–Σ spike detection
    spike = np.abs(err) > tau
    if spike.any():
        idx = np.where(spike)[0][0]          # first spike
        var, bit = divmod(idx, 3)
        # flip that variable if spike at literal slot
        if bit in (0,2):
            truth[var] ^= 1                  # toggle 0↔1
        # after flipping, delete local curvature
        phi[idx] = 0

    # check all clauses
    satisfied = (
        (truth[0]  | ~truth[1] |  truth[2]) &
        (~truth[0] |  truth[1] |  truth[2]) &
        ( truth[0] |  truth[1] | ~truth[2])
    )
    if satisfied:
        print(f"Solved at t={t}:",
              f"x1={truth[0]}, x2={truth[1]}, x3={truth[2]}")
        break
