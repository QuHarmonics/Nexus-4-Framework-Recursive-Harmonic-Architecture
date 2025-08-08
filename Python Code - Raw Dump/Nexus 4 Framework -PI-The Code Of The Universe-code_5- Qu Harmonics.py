import numpy as np

class RecursiveExistenceLoop:
    def __init__(self, beta=0.65, initial=1.0):
        self.beta = beta                # Reflective memory coefficient
        self.history = [initial]        # Recursive state R(t)
        self.existence = [0.0]          # Folded delta E(t)

    def next(self):
        t = len(self.history)
        if t == 1:
            delta = self.history[-1]
        else:
            delta = self.history[-1] - self.history[-2]

        E_t = delta
        R_t = E_t + self.beta * self.history[-1]

        self.existence.append(E_t)
        self.history.append(R_t)

        return R_t, E_t

    def run(self, steps=50):
        results = []
        for _ in range(steps):
            results.append(self.next())
        return results

# Example usage
if __name__ == "__main__":
    rel = RecursiveExistenceLoop()
    output = rel.run(500)
    for t, (r, e) in enumerate(output):
        print(f"t={t+1:2d} | R(t)={r:.4f} | E(t)={e:.4f}")
