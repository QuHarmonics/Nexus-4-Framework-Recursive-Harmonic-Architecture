; Initialize Stack
MOV R0, 1             ; Load Bit 1 (First Value)
MOV R1, 4             ; Load Bit 2 (Second Value)
PUSH R0               ; Push Bit 1 to the stack
PUSH R1               ; Push Bit 2 to the stack

; Step 1: Compute C dynamically
POP R1                ; Pop Bit 2 (4)
POP R0                ; Pop Bit 1 (1)
SUB R2, R1, R0        ; R2 = Bit 2 - Bit 1 (4 - 1 = 3)
PUSH R0               ; Restore Bit 1
PUSH R1               ; Restore Bit 2
PUSH R2               ; Push result (C = 3)
PUSH R2               ; Push C again (Same value, repeated)

; Step 2: Add Last Two Values and Replace Last Position
POP R3                ; Last value (C = 3)
POP R2                ; Second-to-Last value (C = 3)
ADD R4, R3, R2        ; Add last two values (3 + 3 = 6)
PUSH R2               ; Restore Second-to-Last value
PUSH R4               ; Push result (6)

; Step 3: Add Headers Dynamically
POP R3                ; Pop last value (6)
POP R2                ; Pop Second-to-Last value (Bit 2 = 4)
ADD R5, R3, R2        ; Add headers (6 + 4 = 10)
PUSH R3               ; Restore Last Value
PUSH R5               ; Push Headers Result (10)

; Continue Further Operations Dynamically...
