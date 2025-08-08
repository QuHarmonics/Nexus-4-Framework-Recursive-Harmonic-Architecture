#!/usr/bin/env python3
import numpy as np

# 1. Problem Definition
n_var = 3
# 2. Simulation Parameters
block = 3 * n_var
tau   = 1.5
Kp, Ki, Kd = 0.35, 0.02, 0.10
steps = 500

# 3. Build initial curvature φ
phi = np.zeros(block, dtype=float)
packets = np.zeros((3, block), dtype=float)

# Clause → slot mapping
# slot = 3*k + 0 for x_k, +2 for ¬x_k
packets[0, [0, 5, 6]] = 1   # c1: x1, ¬x2, x3
packets[1, [2, 3, 6]] = 1   # c2: ¬x1, x2, x3
packets[2, [0, 3, 8]] = 1   # c3: x1, x2, ¬x3

for p in packets:
    p -= p.mean()
phi += packets.sum(axis=0)

# 4. Initialize controller and state
i_acc = np.zeros(block, dtype=float)
d_prev = np.zeros(block, dtype=float)
truth = np.zeros(n_var, dtype=bool)

def laplacian(v):
    return np.roll(v, 1) - 2*v + np.roll(v, -1)

# 5. Main simulation loop
for t in range(steps):
    # 5.1 Compute error
    err = laplacian(phi)

    # 5.2 PID components
    p_term =  Kp * err
    i_acc  += Ki * err
    d_term =  Kd * (err - d_prev)
    d_prev  = err.copy()

    # 5.3 Apply control
    phi -= (p_term + i_acc + d_term)

    # 5.4 Δ–Σ spike detection & variable toggle
    spike = np.abs(err) > tau
    if spike.any():
        idx = np.where(spike)[0][0]
        var, bit = divmod(idx, 3)
        print(f"\nTick {t}: Spike at slot {idx} (var {var+1}, bit {bit})")
        print("  Before:", truth.astype(int))

        if bit in (0, 2):
            truth[var] = not truth[var]
            print(f"  Toggled x{var+1} →", int(truth[var]))
        else:
            print("  Spacer spike—no toggle")

        phi[idx] = 0.0

    # 5.5 Check for a satisfying assignment *inside* the loop
    #      (only after at least one spike)
    if t > 0 and spike.any():
        c1 = truth[0] or (not truth[1]) or truth[2]
        c2 = (not truth[0]) or truth[1] or truth[2]
        c3 = truth[0] or truth[1] or (not truth[2])
        print(f"Tick {t}: Clause status – c1={c1}, c2={c2}, c3={c3}")
        if c1 and c2 and c3:
            print(f"\nSolved at t={t}: x1={int(truth[0])}, "
                  f"x2={int(truth[1])}, x3={int(truth[2])}")
            break
else:
    print("\nNo solution found in allotted steps.")

print("\nFinal assignment:", truth.astype(int))
