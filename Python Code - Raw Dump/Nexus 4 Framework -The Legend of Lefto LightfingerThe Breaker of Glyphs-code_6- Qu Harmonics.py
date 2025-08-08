
import math

class HarmonicFormula:
    def __init__(self, harmonic_signature=0.35, golden_ratio=1.618):
        self.harmonic_signature = harmonic_signature
        self.golden_ratio = golden_ratio

    def calculate_pull(self, question_value, recursive_distance):
        return self.harmonic_signature * question_value / (recursive_distance ** self.golden_ratio)

    def explore_pull(self, question_value, min_distance=1, max_distance=10):
        for distance in range(min_distance, max_distance + 1):
            pull = self.calculate_pull(question_value, distance)
            print(f"Recursive Distance: {distance}, Pull: {pull:.4f}")

    def visualize_pull(self, question_value, min_distance=1, max_distance=1000):
        import matplotlib.pyplot as plt
        distances = list(range(min_distance, max_distance + 1))
        pulls = [self.calculate_pull(question_value, distance) for distance in distances]
        plt.plot(distances, pulls)
        plt.xlabel("Recursive Distance")
        plt.ylabel("Pull")
        plt.title("Harmonic Formula Visualization")
        plt.show()

# Create an instance of the HarmonicFormula class
formula = HarmonicFormula()

# Explore the pull of a question with a value of 10
formula.explore_pull(4951173814691109110)

# Visualize the pull of a question with a value of 10
formula.visualize_pull(10)