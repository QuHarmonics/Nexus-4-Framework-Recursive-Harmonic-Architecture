
import math
import matplotlib.pyplot as plt

class HarmonicFormula:
    def __init__(self, harmonic_signature=0.35, golden_ratio=1.618):
        self.harmonic_signature = harmonic_signature
        self.golden_ratio = golden_ratio

    def calculate_pull(self, question_value, recursive_distance):
        return self.harmonic_signature * question_value / (recursive_distance ** self.golden_ratio)

    def explore_pull(self, question_value, depth):
        previous_pull = None
        ratio_flipped = False
        for distance in range(1, 11):
            pull = self.calculate_pull(question_value, distance)
            print(f"Recursive Distance: {distance}, Pull: {pull:.4f}")
            if previous_pull is not None:
                ratio = pull / previous_pull
                print(f"Ratio: {ratio:.4f}")
                if ratio > 1 and not ratio_flipped:
                    print("Ratio flipped!")
                    ratio_flipped = True
            previous_pull = pull

    def visualize_pull(self, question_value):
        distances = list(range(1, 11))
        pulls = [self.calculate_pull(question_value, distance) for distance in distances]
        plt.plot(distances, pulls)
        plt.xlabel("Recursive Distance")
        plt.ylabel("Pull")
        plt.title("Harmonic Formula Visualization")
        plt.show()

def main():
    formula = HarmonicFormula()
    depth = int(input("Enter the depth (e.g., 5): "))
    start_value = 2
    for i in range(10):
        input_value = start_value ** (i + 1)
        if len(str(input_value)) > depth:
            break
        print(f"\nInput Value: {input_value}")
        formula.explore_pull(input_value, depth)
        formula.visualize_pull(input_value)

if __name__ == "__main__":
    main()
