Data point 1:
Difficulty Level to ASM

0:  1d                      .byte 0x1d
1:  00 ff                   add    bh,bh
3:  ff                      .byte 0xff

The second instruction (00 ff) is an ADD operation that adds the value of register BH to itself. This effectively doubles the value stored in BH.
