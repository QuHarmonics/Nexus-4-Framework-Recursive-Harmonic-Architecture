; STEP 1: Initialize the stack with the first two values
PUSH 1          ; Push first value onto the stack
PUSH 4          ; Push second value onto the stack

; STEP 2: Compute Var Whole Value (Bit 1 - Bit 2)
MOV R1, 1       ; Load Bit 1 into R1
MOV R2, 4       ; Load Bit 2 into R2
SUB R3, R1, R2  ; Compute R3 = R1 - R2 (Var Whole Value)

; STEP 3: Calculate LEN (Length of current stack)
MOV LEN, 2      ; LEN() of the stack is 2

; STEP 4: Add LEN to the stack LEN times
MOV R4, LEN     ; Store LEN in R4
PUSH R4         ; Add LEN value to the stack
PUSH R4         ; Add LEN value to the stack again

; Final stack after this step: [1, 4, 2, 2]

; STEP 5: Update the value `2` to `1`
; Pointer is initially at the last position (second `2`)
MOV R5, [Stack - 2] ; Load the current pointer value (last `2`)
MOV R6, [Stack - 3] ; Load the value at (Pointer - 1) (value = `5`)
SUB R7, R5, R6      ; Compute R7 = 2 - 1 = 1
MOV [Stack - 2], R7 ; Replace the second `2` with `1`

; Final stack after this step: [1, 4, 2, 1]

; STEP 6: Update the stack value at Pointer
MOV R1, [Stack - 4] ; Load Bit 0 (value at Stack - 4 = 1)
MOV R2, [Stack - 3] ; Load Bit 1 (value at Stack - 3 = 4)
ADD R8, R1, R2      ; Compute R8 = Bit 0 + Bit 1 (1 + 4 = 5)
MOV [Stack - 2], R8 ; Replace the value at Pointer with R8

; Final stack after this step: [1, 4, 1, 5]

; STEP 7: Calculate the next value in the sequence
; Use the current pointer value and the value at (Pointer - Pointer value - 1)
MOV CurrentPointer, [Stack - 2]  ; Load current pointer value (5)
SUB R9, CurrentPointer, 1        ; Compute (Pointer - 1)
MOV R10, [Stack - R9]            ; Load value at (Pointer - R9) (value = 4)
ADD R11, R10, CurrentPointer     ; Add value at (Pointer - R9) + CurrentPointer
PUSH R11                         ; Push the result onto the stack

; Final stack after this step: [1, 4, 1, 5, 9]
