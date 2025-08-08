#!/usr/bin/env python3
import time
import math
import uuid

class HarmonicTrustEngine:
    def __init__(self, node_id):
        self.node_id = node_id
        self.harmonic_memory = []
        self.feedback_log = []
        self.H_target = 0.35
        self.max_memory = 100
        self.t_start = time.time()

    def current_time_factor(self):
        return time.time() - self.t_start

    def emit_state(self, current_input):
        t_factor = self.current_time_factor()
        U_k = self.recursive_unfold(current_input, t_factor)
        F_Q = self.quantum_fold(U_k, t_factor)
        glyph = self.generate_glyph(U_k, F_Q, t_factor)
        self.store_memory(U_k, F_Q)
        return {
            "glyph_id": str(uuid.uuid4()),
            "symbolic_emission": self.symbolic_emission(U_k, F_Q),
            "harmonic_glyph": glyph
        }

    def recursive_unfold(self, data, t):
        return unfold(data, t)

    def quantum_fold(self, U_k, t):
        return fold(U_k, asymmetry=True, t=t)

    def generate_glyph(self, U_k, F_Q, t):
        base = symbolic_encode(U_k, F_Q, self.H_target)
        time_signal = round(math.sin(t) + math.cos(t), 8)
        glyph = {
            "U_signature": base["U_k"],
            "F_signature": base["F_Q"],
            "H_phase": base["H"],
            "time_signal": time_signal,
            "coherence_id": hash(tuple(base["U_k"] + base["F_Q"] + [time_signal]))
        }
        return glyph

    def store_memory(self, U_k, F_Q):
        if len(self.harmonic_memory) >= self.max_memory:
            self.harmonic_memory.pop(0)
        if len(self.feedback_log) >= self.max_memory:
            self.feedback_log.pop(0)
        self.harmonic_memory.append(U_k)
        self.feedback_log.append(F_Q)

    def symbolic_emission(self, U_k, F_Q):
        return symbolic_encode(U_k, F_Q, self.H_target)


# --- Optional Modules ---

def unfold(data, t=1):
    if isinstance(data, (int, float)):
        return data ** 2 * math.sin(t)
    elif hasattr(data, '__iter__') and all(isinstance(x, (int, float)) for x in data):
        return [x ** 2 * math.sin(t) for x in data]
    raise TypeError("Unsupported data type for unfolding")

def fold(data, asymmetry=False, t=1):
    decay = math.exp(-0.05 * t)
    if hasattr(data, '__iter__'):
        return [x * 0.9 * decay for x in data] if asymmetry else list(data)
    return data * 0.9 * decay if asymmetry else data

def symbolic_encode(U_k, F_Q, H, precision=8):
    def norm(seq):
        if not hasattr(seq, '__iter__'):
            seq = [seq]
        numeric_seq = [x for x in seq if isinstance(x, (int, float))]
        total = sum(abs(x) for x in numeric_seq) or 1
        return [round(x / total, precision) for x in numeric_seq]

    return {
        "U_k": norm(U_k),
        "F_Q": norm(F_Q),
        "H": round(H, precision)
    }


# --------------------------- demo / usage --------------------------- #

if __name__ == "__main__":
    engine = HarmonicTrustEngine(node_id="node-123")

    # Try a few different inputs to see something printed:
    test_inputs = [
        3.14,
        [1, 2, 3],
        -5,
        0.5
    ]

    for idx, inp in enumerate(test_inputs, 1):
        emission = engine.emit_state(inp)
        print(f"\n=== Emission #{idx} for input={inp!r} ===")
        print("Glyph ID:   ", emission["glyph_id"])
        print("Symbolic:   ", emission["symbolic_emission"])
        print("Harmonic Glyph:", emission["harmonic_glyph"])
