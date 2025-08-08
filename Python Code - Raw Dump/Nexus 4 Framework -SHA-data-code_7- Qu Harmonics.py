def fold_step(stack, pos, payload):
    """
    stack  : list[int]          e.g. [1,4,7,10,...]
    pos    : insertion index
    payload: iterable[int]      e.g. [11,13]
    """
    C = len(payload)            # here it's 2
    stack[pos:pos] = [None]*C   # reserve fold-space
    stack[pos+C:pos+C] = stack[pos+2*C:]
    del stack[pos+2*C:]
    stack[pos:pos+C] = payload  # back-fill

raw   = [1,4,7,10,13,16]
fold_step(raw, 2, [11,13])
print(raw)   # -> [1, 4, 11, 13, 17, 20]
