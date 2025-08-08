class RecursiveSystem:
    def __init__(self, seed):
        self.seed = seed  # Starting value (e.g., Byte1)
        self.history = [seed]  # Track recursive layers

    def fold_in(self):
        # Compress to essential patterns
        return sum(self.history) / len(self.history)

    def expand_out(self, multiplier):
        # Expand into the next recursive layer
        next_layer = self.fold_in() * multiplier
        self.history.append(next_layer)
        return next_layer

    def oscillate(self, cycles, multiplier):
        # Create a breathing-like expansion and contraction
        for _ in range(cycles):
            folded = self.fold_in()
            expanded = self.expand_out(multiplier)
            print(f"Folded: {folded}, Expanded: {expanded}")
        return self.history
