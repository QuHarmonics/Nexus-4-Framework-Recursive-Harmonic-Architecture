import numpy as np
import matplotlib.pyplot as plt

# --- Mark1 Harmonic Engine ---

class Mark1Engine:
    def __init__(self, potential_values, actualized_values, target_H=0.35, learning_rate=0.01):
        self.P = np.array(potential_values, dtype=float)
        self.A = np.array(actualized_values, dtype=float)
        self.H_target = target_H
        self.lr = learning_rate
        self.history = []

    def compute_H(self):
        P_sum = np.sum(self.P)
        A_sum = np.sum(self.A)
        return A_sum / P_sum if P_sum != 0 else 0

    def step(self):
        """Adjust A_i values proportionally to minimize delta H"""
        H = self.compute_H()
        delta_H = self.H_target - H
        adjustment = self.lr * delta_H * self.P
        self.A += adjustment
        self.history.append(H)

    def run(self, steps=100):
        for _ in range(steps):
            self.step()

    def plot_convergence(self):
        plt.figure(figsize=(10, 5))
        plt.plot(self.history, label="Harmonic Ratio H(t)", linewidth=2)
        plt.axhline(self.H_target, color='r', linestyle='--', label=f"Target H = {self.H_target}")
        plt.title("Mark1 Engine Harmonic Convergence")
        plt.xlabel("Time Step")
        plt.ylabel("H = ∑A / ∑P")
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.show()


# Example: Random potential/actual values
np.random.seed(42)
P_example = np.random.uniform(1, 10, size=25)
A_example = np.random.uniform(0.1, 2.0, size=25)

# Initialize and run the engine
mark1 = Mark1Engine(P_example, A_example)
mark1.run(steps=100)
mark1.plot_convergence()
