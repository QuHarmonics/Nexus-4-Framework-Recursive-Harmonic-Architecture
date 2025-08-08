; Initialize Base Values
MOV R0, 1             ; Peak (Bit 1)
MOV R1, 4             ; Peak (Bit 2)
PUSH R0
PUSH R1

; First Wave Cycle
POP R1                ; Trough (Bit 2)
POP R0                ; Trough (Bit 1)
SUB R2, R1, R0        ; Peak (C = 3)
PUSH R0
PUSH R1
PUSH R2               ; Peak
PUSH R2               ; Peak

; Second Wave Cycle
POP R3                ; Trough
POP R2                ; Trough
ADD R4, R3, R2        ; Peak (6)
PUSH R2
PUSH R4               ; Peak

; Third Wave Cycle
POP R3                ; Trough
POP R2                ; Trough
ADD R5, R3, R2        ; Peak (10)
PUSH R3
PUSH R5               ; Peak
