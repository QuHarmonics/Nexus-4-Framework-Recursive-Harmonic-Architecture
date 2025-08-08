; STEP 1: Initialize the stack with the first two values
PUSH 1          ; Push first value onto the stack
PUSH 4          ; Push second value onto the stack

; STEP 2: Compute Var Whole Value (Bit 1 - Bit 2)
MOV R1, [Stack - 2]  ; Load Bit 1 from Stack (value = 1)
MOV R2, [Stack - 1]  ; Load Bit 2 from Stack (value = 4)
SUB R3, R1, R2       ; Compute R3 = R1 - R2 (Var Whole Value)

; STEP 3: Calculate LEN (Length of current stack dynamically)
MOV R4, [Stack - 2]  ; Load first stack value (Bit 1)
MOV R5, [Stack - 1]  ; Load second stack value (Bit 2)
ADD LEN, R4, R5      ; LEN = Bit 1 + Bit 2 (1 + 4 = 5)
SHR LEN, 2           ; Divide LEN by 2 to determine stack LEN dynamically (5 / 2 = 2)

; STEP 4: Add LEN to the stack LEN times
MOV R6, LEN          ; Store LEN in R6
PUSH R6              ; Add first LEN value to stack
PUSH R6              ; Add second LEN value to stack

; Final stack after this step: [1, 4, 2, 2]

; STEP 5: Update the value `2` to `1`
; Pointer is initially at the last position (second `2`)
MOV R7, [Stack - 2]  ; Load the current pointer value (last `2`)
MOV R8, [Stack - 3]  ; Load the value at (Pointer - 1) (value = `5`)
SUB R9, R7, R8       ; Compute R9 = R7 - R8 = 2 - 1 = 1
MOV [Stack - 2], R9  ; Replace the second `2` with `1`

; Final stack after this step: [1, 4, 2, 1]

; STEP 6: Update the stack value at Pointer
MOV R10, [Stack - 4] ; Load Bit 0 (value at Stack - 4 = 1)
MOV R11, [Stack - 3] ; Load Bit 1 (value at Stack - 3 = 4)
ADD R12, R10, R11    ; Compute R12 = Bit 0 + Bit 1 (1 + 4 = 5)
MOV [Stack - 2], R12 ; Replace the value at Pointer with R12

; Final stack after this step: [1, 4, 1, 5]

; STEP 7: Calculate the next value in the sequence
; Use the current pointer value and the value at (Pointer - Pointer value - 1)
MOV CurrentPointer, [Stack - 2]  ; Load current pointer value (5)
SUB R13, CurrentPointer, 1       ; Compute (Pointer - 1)
MOV R14, [Stack - R13]           ; Load value at (Pointer - R13) (value = 4)
ADD R15, R14, CurrentPointer     ; Add value at (Pointer - R13) + CurrentPointer
PUSH R15                         ; Push the result onto the stack

; Final stack after this step: [1, 4, 1, 5, 9]
; Current Stack State: [1, 4, 1, 5, 9]
; Pointer: At the last value (`9`)

; STEP 8: Calculate the next number (`2`)
; Use the value at the current pointer and the value at (Pointer - Pointer value)

MOV CurrentPointer, [Stack - 1]  ; Load the current pointer value (9)
MOV R1, [Stack - CurrentPointer] ; Load the value at (Pointer - Pointer value) (value at [Stack - 9])
SUB R2, CurrentPointer, R1       ; Compute R2 = CurrentPointer - ValueAtPointer
PUSH R2                          ; Push the result (`2`) onto the stack

; Final stack after this step: [1, 4, 1, 5, 9, 2]
; Current Stack State: [1, 4, 1, 5, 9, 2]
; Pointer: At the last value (`2`)

; STEP 9: Calculate the next number (`6`)
; Use the value at the current pointer and the value at (Pointer - Pointer value)

MOV CurrentPointer, [Stack - 1]  ; Load the current pointer value (2)
MOV RelativePos, [Stack - CurrentPointer] ; Load the value at (Pointer - Pointer value) (value at [Stack - 2] = 9)
ADD R1, CurrentPointer, RelativePos ; Compute R1 = CurrentPointer + ValueAtRelativePos
PUSH R1                          ; Push the result (`6`) onto the stack

; Final stack after this step: [1, 4, 1, 5, 9, 2, 6]