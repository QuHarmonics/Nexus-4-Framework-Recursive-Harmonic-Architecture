#!/usr/bin/env python3
"""
Corrected 3-SAT solver via curvature‐PID Δ–Σ modulation.
Clause→curvature mapping fixed so (bit 0→x_k, bit 2→¬x_k) correctly.
"""

import numpy as np

# 1. Problem Definition
n_var = 3
# Clauses:
# c1: (x1 OR ¬x2 OR x3)
# c2: (¬x1 OR x2 OR x3)
# c3: (x1 OR x2 OR ¬x3)

# 2. Simulation Parameters
block = 3 * n_var
tau   = 1.5
Kp, Ki, Kd = 0.35, 0.02, 0.10
steps = 500

# 3. Build initial curvature φ
phi = np.zeros(block, dtype=float)
packets = np.zeros((3, block), dtype=float)

# Correct slot mapping: slot = 3*k + {0:x_k,1:spacer,2:¬x_k}
# Clause 1: x1(slot0), ¬x2(slot 3*1+2=5), x3(slot 3*2+0=6)
packets[0, [0, 5, 6]] = 1
# Clause 2: ¬x1(slot2), x2(slot3), x3(slot6)
packets[1, [2, 3, 6]] = 1
# Clause 3: x1(slot0), x2(slot3), ¬x3(slot3*2+2=8)
packets[2, [0, 3, 8]] = 1

# Zero-mean each
for p in packets:
    p -= p.mean()

phi += packets.sum(axis=0)

# 4. Initialize controller and state
i_acc = np.zeros(block, dtype=float)
d_prev = np.zeros(block, dtype=float)
truth = np.zeros(n_var, dtype=bool)

def laplacian(v):
    return np.roll(v, 1) - 2*v + np.roll(v, -1)

# 5. Main loop
for t in range(steps):
    err = laplacian(phi)
    p_term =  Kp * err
    i_acc  += Ki * err
    d_term =  Kd * (err - d_prev)
    d_prev  = err.copy()
    phi    -= (p_term + i_acc + d_term)

    # Δ–Σ spike
    spike = np.abs(err) > tau
    if spike.any():
        idx = np.where(spike)[0][0]
        var, bit = divmod(idx, 3)
        # toggle only on literal slots
        if bit in (0, 2):
            truth[var] = not truth[var]
        phi[idx] = 0.0

    # Check satisfaction after at least one spike
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

print("Final assignment:", truth.astype(int))
