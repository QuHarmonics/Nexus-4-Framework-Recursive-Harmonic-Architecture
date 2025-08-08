import numpy as np

# --- 1  encode sequence ---------------------------------
seq = "AGCTGACT"
P = {'A':2.0,'T':2.0,'G':3.0,'C':3.0}
A = {'A':1.9,'T':1.8,'G':2.7,'C':2.6}
P_vec = np.array([P[b] for b in seq])
A_vec = np.array([A[b] for b in seq])

# --- 2  global parameters -------------------------------
H, F_factor, t = 0.35, 0.5, 1
exp_term = np.exp(H*F_factor*t)

# --- 3  fold --------------------------------------------
fold_contrib = (P_vec / A_vec) * exp_term
F_Q = fold_contrib.sum()
print("Folded scalar:", F_Q)

# --- 4  unfold ------------------------------------------
theta = np.arange(len(seq)) * np.pi/2
U_Q = (fold_contrib * np.cos(theta)).sum()
print("Unfolded checksum:", U_Q)
