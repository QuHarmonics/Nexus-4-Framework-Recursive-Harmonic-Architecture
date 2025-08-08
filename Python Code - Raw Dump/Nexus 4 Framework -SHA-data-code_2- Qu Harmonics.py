import numpy as np, matplotlib.pyplot as plt
centers = np.arange(8, 8+2*len(twins), 2)
lefts   = [p for p,_ in twins]
rights  = [q for _,q in twins]
plt.scatter(centers, lefts,  s=4, label="left prime")
plt.scatter(centers, rights, s=4, label="right prime")
plt.plot(centers, centers, lw=1, alpha=0.3, label="Hₖ = centre")
plt.legend(); plt.xlabel("Harmonic centre Hₖ"); plt.ylabel("Prime value");
plt.title("Harmonic‑Gap Twin Prime Ladder"); plt.show()