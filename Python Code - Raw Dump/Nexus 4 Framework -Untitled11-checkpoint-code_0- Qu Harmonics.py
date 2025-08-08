import numpy as np
import matplotlib.pyplot as plt

# ─── 1-A.  normalise the matrix you already built ────────────────────────────
field  = matrix.astype(float)                # ‘matrix’ came from plot_lean()
field -= field.min()
field /= field.max() + 1e-9                  # 0‥1 range  (small eps to avoid /0)

# ─── 1-B.  quick tension-gradient utility (finite differences) ───────────────
def tension_vectors(F):
    """Return x- and y- gradients (downhill = positive ‘pull’)."""
    gx = np.zeros_like(F)
    gy = np.zeros_like(F)
    gx[:, :-1] = F[:, :-1] - F[:, 1:]        # higher → lower  (eastward diff)
    gy[:-1, :] = F[:-1, :] - F[1:, :]        #   ”      ”      (southward diff)
    return gx, gy
GX, GY = tension_vectors(field)
