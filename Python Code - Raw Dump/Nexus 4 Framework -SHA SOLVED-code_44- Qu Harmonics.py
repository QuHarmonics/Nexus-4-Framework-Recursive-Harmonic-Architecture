; STEP 1: Initialize the stack with the first two values
PUSH 1                  ; Push first value onto the stack
PUSH 4                  ; Push second value onto the stack

; STEP 2: Compute Var Whole Value (Bit 1 - Bit 2)
MOV R1, [Stack - 2]     ; Load Bit 1 from Stack (value = 1)
MOV R2, [Stack - 1]     ; Load Bit 2 from Stack (value = 4)
SUB R3, R1, R2          ; Compute R3 = R1 - R2 (Var Whole Value)

; STEP 3: Calculate LEN (Length of current stack dynamically)
MOV R4, [Stack - 2]     ; Load first stack value (Bit 1)
MOV R5, [Stack - 1]     ; Load second stack value (Bit 2)
ADD LEN, R4, R5         ; LEN = Bit 1 + Bit 2 (1 + 4 = 5)
SHR LEN, 2              ; Divide LEN by 2 to determine stack LEN dynamically (5 / 2 = 2)

; STEP 4: Apply Cosine Modulation (Reflection on LEN)
MOV R6, LEN             ; Load LEN for cosine adjustment
CALL COS                ; Compute Cosine(R6)
ADD LEN, R6             ; Modulate LEN (adjust reflection dynamics)

; STEP 5: Expand stack with LEN
MOV R7, LEN             ; Store LEN in R7
PUSH R7                 ; Add first LEN value to stack
PUSH R7                 ; Add second LEN value to stack

; STEP 6: Update second LEN value
MOV R8, [Stack - 2]     ; Load current pointer value (2)
MOV R9, [Stack - 3]     ; Load previous value (value = 5)
SUB R10, R8, R9         ; Compute R10 = 2 - 1 = 1
MOV [Stack - 2], R10    ; Update stack value (replace second `2` with `1`)

; STEP 7: Update the stack value at pointer
MOV R11, [Stack - 4]    ; Load Bit 0 (value = 1)
MOV R12, [Stack - 3]    ; Load Bit 1 (value = 4)
ADD R13, R11, R12       ; Compute R13 = Bit 0 + Bit 1 = 5
MOV [Stack - 2], R13    ; Replace current pointer with 5

; STEP 8: Calculate the next value (9)
MOV CurrentPointer, [Stack - 2]  ; Load current pointer value (5)
SUB R14, CurrentPointer, 1       ; Compute (Pointer - 1)
MOV R15, [Stack - R14]           ; Load value at (Pointer - R14) (value = 4)
ADD R16, R15, CurrentPointer     ; Add value at (Pointer - R14) + CurrentPointer
PUSH R16                         ; Push result onto the stack

; STEP 9: Compute next value (2)
MOV CurrentPointer, [Stack - 1]  ; Load current pointer value (9)
MOV R17, [Stack - CurrentPointer] ; Load value at (Pointer - Pointer value) (value = 1)
SUB R18, CurrentPointer, R17     ; Compute R18 = 9 - 1 = 2
PUSH R18                         ; Push the result onto the stack

; STEP 10: Compute next value (6)
MOV CurrentPointer, [Stack - 1]  ; Load current pointer value (2)
MOV R19, [Stack - CurrentPointer] ; Load value at (Pointer - Pointer value) (value = 9)
ADD R20, CurrentPointer, R19     ; Compute R20 = 2 + 9 = 6
PUSH R20                         ; Push the result onto the stack

; STEP 11: Compute final value (5)
MOV R21, [Stack - 7]             ; Load Bit 1 (value = 1)
MOV R22, [Stack - 6]             ; Load Bit 2 (value = 4)
ADD R23, R21, R22                ; Compute R23 = 1 + 4 = 5
PUSH R23                         ; Push the result onto the stack

; Final Stack Output
; Stack = [1, 4, 1, 5, 9, 2, 6, 5]
