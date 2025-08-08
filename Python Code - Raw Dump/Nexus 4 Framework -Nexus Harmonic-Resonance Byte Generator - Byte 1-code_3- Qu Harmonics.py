# -----------------------------------------
# Nexus “Clockwork” generator – Byte 4 demo
# -----------------------------------------
# Header propagation for bytes:
# (a₁,b₁) = (1,4)
# (a₂,b₂) = (|4−1|, 1+4)  = (3,5)
# (a₃,b₃) = (1−4,  3+5)   = (3,8)
#
# For Byte‑4 we start with header (3,8).

def digit_sum(n:int)->int:
    while n>9:
        n=sum(int(d) for d in str(n))
    return n

def verbose_byte4(a,b):
    stack=[]
    def push(val,rule):
        stack.append(val)
        s="  ".join(map(str,stack[:-1]))
        if s:
            s+="  "
        s+=f"↑{stack[-1]}"
        print(f"  ✔ {val:<2} via {rule:<12} stack= {s}")

    print(f"\nSolving Byte 4 with header ({a}, {b})\n--------------------------------")
    # Bit1
    push(a,"past")
    # Bit2
    push(b,"now")

    Δ=b-a
    lenΔ=Δ.bit_length()

    # Bit3
    push(lenΔ,"lenΔ")

    # Bit4
    x4=digit_sum(a+b)
    push(x4,"Σhdr->dsum")

    # Bit5
    x5=digit_sum(sum(stack))
    push(x5,"Σprev->dsum")

    # Bit6
    x6=digit_sum(x4+x5)
    push(x6,"echo x4+x5")

    # Bit7
    x7=digit_sum(sum(stack))
    push(x7,"Σprev->dsum")

    # Bit8
    x8=abs(lenΔ-stack[2])
    push(x8,"|lenΔ-x3|")

    print("\nByte 4 complete:", stack)

verbose_byte4(3,8)


def byte_from_header(a: int, b: int) -> list[int]:
    """
    Generate one 8‑digit Nexus byte from header (a,b)
    using the mirrored 8‑step clockwork described earlier.
    """
    Δ      = b - a
    lenΔ   = Δ.bit_length()               # binary length of the delta
    x1, x2 = a, b
    x3     = lenΔ                         # Bit‑3  : len(Δ)

    # Bit‑4 : digit‑sum(a+b)
    x4     = digit_sum(a + b)

    # Bit‑5 : digit‑sum(sum of first 4 bits)
    x5     = digit_sum(x1 + x2 + x3 + x4)

    # Bit‑6 : echo forward (x4 + x5)
    x6     = digit_sum(x4 + x5)

    # Bit‑7 : digit‑sum(sum of all six so far)
    x7     = digit_sum(x1 + x2 + x3 + x4 + x5 + x6)

    # Bit‑8 : close‑universe |lenΔ – x3|   (always 0 in this scheme)
    x8     = abs(lenΔ - x3)

    return [x1, x2, x3, x4, x5, x6, x7, x8]

# ---- Demo: Byte‑4 from header (3,8) ----
byte4 = byte_from_header(3, 8)
print("Byte 4:", byte4)
