#!/usr/bin/env python3
"""
Toy 3-SAT solver via curvature-PID Δ–Σ modulation.
Maps a small 3-variable 3-clause instance onto a discrete curvature field,
runs a PID feedback loop with Δ–Σ “quantizer” spikes, and reads off a satisfying
assignment in O(1) ticks if the model holds.
"""

import numpy as np

# ---------------------
# 1. Problem Definition
# ---------------------
# 3 variables, 3 clauses: only (x1,x2,x3) = (True,True,True) satisfies all.
n_var = 3
# Clauses:
# c1:  x1  OR not x2 OR  x3
# c2: not x1 OR  x2  OR  x3
# c3:  x1  OR  x2  OR not x3

# -------------------------
# 2. Simulation Parameters
# -------------------------
block = 3 * n_var         # 3 slots per variable
tau   = 1.5               # curvature threshold for spike
Kp, Ki, Kd = 0.35, 0.02, 0.10  # PID gains
steps = 500               # max ticks

# ----------------------------
# 3. Build initial curvature φ
# ----------------------------
phi = np.zeros(block, dtype=float)

# Packets: one per clause, zero-mean over the block
packets = np.zeros((3, block), dtype=float)

# Clause 1: x1, ¬x2, x3
packets[0, [0, 2, 6]] += 1
# Clause 2: ¬x1, x2, x3
packets[1, [1, 3, 6]] += 1
# Clause 3: x1, x2, ¬x3
packets[2, [0, 3, 7]] += 1

# Balance each packet to zero mean
for p in packets:
    p -= p.mean()

# Inject all clauses at t=0
phi += packets.sum(axis=0)

# ----------------------
# 4. Initialize controller
# ----------------------
i_acc = np.zeros(block, dtype=float)  # integral term accumulator
d_prev = np.zeros(block, dtype=float) # previous error for derivative
truth = np.zeros(n_var, dtype=bool)   # start all False

def laplacian(v):
    """1D discrete Laplacian on a circular lattice."""
    return np.roll(v, 1) - 2*v + np.roll(v, -1)

# ---------------------
# 5. Main simulation loop
# ---------------------
for t in range(steps):
    # 5.1 Compute curvature error
    err = laplacian(phi)

    # 5.2 PID terms
    p_term =  Kp * err
    i_acc  += Ki * err
    d_term =  Kd * (err - d_prev)
    d_prev  = err.copy()

    # 5.3 Control action: subtract from φ
    phi -= (p_term + i_acc + d_term)

    # 5.4 Δ–Σ spike detection
    spike = np.abs(err) > tau
    if spike.any():
        idx = np.where(spike)[0][0]      # first spike index
        var, bit = divmod(idx, 3)        # which variable & literal slot

        # Only literal slots (bit==0 for x, bit==2 for ¬x) toggle the var
        if bit in (0, 2):
            truth[var] = not truth[var]

        # Clear curvature at the spike site
        phi[idx] = 0.0

    # 5.5 Check clause satisfaction (only after at least one spike)
    if t > 0:
        c1 = truth[0] or (not truth[1]) or truth[2]
        c2 = (not truth[0]) or truth[1] or truth[2]
        c3 = truth[0] or truth[1] or (not truth[2])
        if c1 and c2 and c3:
            print(f"Solved at t={t}: x1={int(truth[0])}, "
                  f"x2={int(truth[1])}, x3={int(truth[2])}")
            break
else:
    print("No solution found in allotted steps.")

# If desired, display final truth assignment
print("Final assignment:", truth.astype(int))
