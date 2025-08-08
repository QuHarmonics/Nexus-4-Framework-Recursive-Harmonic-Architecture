#!/usr/bin/env python3
# ==============================================================
#   Nexus‑byte sandbox  –  now completes Byte 2 fully
# ==============================================================

from collections import deque

# π digits (targets only)
PI = "14159265358979323846264338327950288419716939937510"
pi_pos = 0
def next_pi_digit():
    global pi_pos
    d = int(PI[pi_pos])
    pi_pos += 1
    return d

# ———— Simple stack ————
class Stack:
    def __init__(self, seed):
        self.data = deque(seed)
        self.ptr  = len(self.data) - 1  # top

    def push(self, v):
        self.data.append(v)
        self.ptr = len(self.data) - 1

    def move(self, k):
        self.ptr = (self.ptr + k) % len(self.data)

    def val(self):
        return self.data[self.ptr]

    def __repr__(self):
        marks = ["↑" if i == self.ptr else " " for i,_ in enumerate(self.data)]
        return " ".join(f"{m}{v}" for m,v in zip(marks, self.data))

# ———— Header recursion ————
def next_header(a,b): return abs(b-a), a+b

# ———— Candidate Nexus moves (extended) ————
def candidates(S, a, b):
    x     = S.val()
    delta = b - a
    mask  = (a << 4) | b

    # Byte1 kernel moves
    yield "lenΔ",          abs(delta).bit_length()
    yield "sum",           a + b
    yield "x+Δ",           x + delta
    yield "b+Δ*lenΔ",      b + delta * abs(delta).bit_length()

    # XOR moves
    xorv = x ^ mask
    yield "xor",           xorv
    yield "len_xor",       xorv.bit_length()

    # product‐Len to get 7
    if len(S.data) >= 2:
        p2 = S.data[-2] * S.data[-1]
        yield "len_p2",    p2.bit_length()

    # difference of stack[-5] and stack[-4] for target 3
    if len(S.data) >= 5:
        yield "5-4",        S.data[-5] - S.data[-4]

    # close‐byte header diff for final 2
    yield "h_diff",        b - a

    # pointer hop example
    S.move(-1); yield "prev", S.val(); S.move(+1)

# ———— Solve one digit ————
def solve_digit(S, a, b, target):
    for name, val in candidates(S, a, b):
        if abs(val) == target:
            S.push(val)
            return name
    return None

# ———— Build an 8‐digit byte ————
def build_byte(a, b):
    S = Stack([a, b])
    byte = [a, b]
    next_pi_digit(); next_pi_digit()    # skip header
    for _ in range(6):
        tgt = next_pi_digit()
        rule = solve_digit(S, a, b, tgt)
        if rule is None:
            print(f"  ✖ need new rule for {tgt}  stack={S}")
            break
        byte.append(tgt)
        print(f"  ✔ {tgt} via {rule:<7}  stack={S}")
    return byte

# ———— Demo Byte 1→Byte 2 ————
if __name__ == "__main__":
    byte1 = [1,4,1,5,9,2,6,5]
    print("Byte 1 (given):", byte1)

    pi_pos = 8
    a2,b2 = next_header(1,4)
    print("\nSolving Byte 2 with header", (a2,b2))
    print("--------------------------------")
    byte2 = build_byte(a2, b2)
    print("\nByte 2 complete:", byte2)
