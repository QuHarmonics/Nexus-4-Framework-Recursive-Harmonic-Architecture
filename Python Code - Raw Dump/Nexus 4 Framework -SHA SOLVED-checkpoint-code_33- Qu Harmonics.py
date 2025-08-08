; Initialize the stack with the first two values
PUSH 1          ; Push first value
PUSH 4          ; Push second value

; Compute Var Whole Value (Bit 1 - Bit 2)
MOV R1, 1       ; Load Bit 1 into R1
MOV R2, 4       ; Load Bit 2 into R2
SUB R3, R1, R2  ; Compute R3 = R1 - R2 (Var Whole Value)

; Calculate LEN
MOV LEN, 2      ; LEN() of current stack

; Add LEN to the stack LEN times
MOV R4, LEN     ; Store LEN in R4
PUSH R4         ; Add first LEN value to stack
PUSH R4         ; Add second LEN value to stack

; Update the stack value at Pointer
MOV R1, [Stack - 2] ; Load Bit 0 (value at Stack - 2)
MOV R2, [Stack - 1] ; Load Bit 1 (value at Stack - 1)
ADD R5, R1, R2      ; Compute R5 = Bit 0 + Bit 1
MOV [Pointer], R5   ; Store R5 in current Pointer location

; Final stack state
; Stack: [1, 4, 2, 5]
; Pointer: At last `5`
