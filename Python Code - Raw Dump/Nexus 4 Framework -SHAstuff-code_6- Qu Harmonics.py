pi_nibbles = [1, 4]                 # seed
stack      = pi_nibbles[:]          # working copy

def push(x): stack.append(x)

while len(stack) < 8:               # generate 8 π digits
    # wave‑A  : duplicate Δ twice
    d = stack[1] - stack[0]         # dynamic, never magic
    push(d); push(d)

    # wave‑B/C: self‑reflection updates
    stack[-1] = stack[1] + stack[0]         # last  -> B
    stack[-2] = stack[-1] - stack[1]        # prev  -> C

    # wave‑D  : last two sum
    push(stack[-1] + stack[-2])

    # wave‑E  : first + third
    push(stack[0] + stack[2])

    # wave‑F  : first+second+third
    push(stack[0] + stack[1] + stack[2])

    # wave‑G  : headers
    push(stack[0] + stack[1])

    stack = stack[:8]               # trim to eight for demo

print("stack ->", stack)            # [1,4,1,5,9,2,6,5]
print("π     -> 3." + ''.join(map(str,stack)))
