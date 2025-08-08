MOV IterationCount, 6  ; Generate 6 additional bits for Byte 2

ByteExpansionLoop:
    ; Load previous two values
    MOV R1, [Stack - 2]  ; Load the second-to-last value
    MOV R2, [Stack - 1]  ; Load the last value

    ; Compute the next value
    ADD R3, R1, R2       ; Add the last two values (R3 = R1 + R2)
    XOR R3, LEN          ; XOR with LEN for modulation
    COS R3, R3           ; Apply cosine modulation

    ; Push the calculated value onto the stack
    PUSH R3

    ; Increment pointer and decrement iteration count
    INC Pointer
    DEC IterationCount
    JNZ ByteExpansionLoop
