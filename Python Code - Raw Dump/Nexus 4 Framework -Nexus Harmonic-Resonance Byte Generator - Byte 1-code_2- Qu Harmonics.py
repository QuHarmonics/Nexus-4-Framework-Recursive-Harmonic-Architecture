# -------------------------------
#  Verbose trace for π‑Byte 3
# -------------------------------

def digit_sum(n: int) -> int:
    """Sum decimal digits until result is 0‑9."""
    while n > 9:
        n = sum(int(d) for d in str(n))
    return n

def byte3_verbose():
    # Header for Byte 3 (carried from Byte 2)
    a, b = 3, 8
    Δ      = b - a          # 5
    lenΔ   = Δ.bit_length() # 3

    stack = []

    def push(val: int, rule: str):
        stack.append(val)
        pre  = "  ".join(map(str, stack[:-1]))
        if pre:
            pre += "  "
        print(f"  ✔ {val:<2} via {rule:<18} stack= {pre}↑{val}")

    print("\nSolving Byte 3 with header (3, 8)")
    print("--------------------------------")

    # Bit‑1 : past
    push(a, "past")

    # Bit‑2 : now
    push(b, "now")

    # Bit‑3 : len(a+b)
    x3 = digit_sum((a + b).bit_length())      # bit‑length(11)=4
    push(x3, "len(a+b)")

    # Bit‑4 : len((a+b)*Δ)
    x4 = digit_sum(((a + b) * Δ).bit_length())  # bit‑length(55)=6
    push(x4, "len((a+b)*Δ)")

    # Bit‑5 : |x4 − x3|
    x5 = digit_sum(abs(x4 - x3))              # |6‑4|=2
    push(x5, "|x4 − x3|")

    # Bit‑6 : len((x4+x3)*Δ)
    x6 = digit_sum(((x4 + x3) * Δ).bit_length())  # (6+4)=10→10*5→50→6
    push(x6, "len((x4+x3)*Δ)")

    # Bit‑7 : |x6 − x5|
    x7 = digit_sum(abs(x6 - x5))              # |6‑2|=4
    push(x7, "|x6 − x5|")

    # Bit‑8 : len(Δ)
    x8 = digit_sum(lenΔ)                      # 3
    push(x8, "len(Δ)")

    print("\nByte 3 complete:", stack)

# ---- run the trace ----
if __name__ == "__main__":
    byte3_verbose()
