import matplotlib.pyplot as plt

def simulate_fold_steps(steps=30):
    widths = []
    C = 2
    for i in range(steps):
        if i % 3 == 0:
            C = 2 * C  # Reflective hinge
        else:
            C = max(2, C // 2)
        widths.append(C)
    return widths

data = simulate_fold_steps()
plt.plot(data, marker='o')
plt.title("Dyadic Fold-Space Widths (2 ↔ 4 pattern)")
plt.xlabel("Step")
plt.ylabel("Reserve Width (C)")
plt.grid(True)
plt.show()
