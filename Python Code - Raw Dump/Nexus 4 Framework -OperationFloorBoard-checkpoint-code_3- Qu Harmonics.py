TESTING
1. 9f0a1b2c3d4e5f6a7b8c9d0e1f2a3b4c5d6e7f8a9b0c1d2e3f4a5b6c7d8e9f0a

2.
lahf
or bl, BYTE PTR [ebx]
sub al, 0x3d
dec esi
pop edi
push 0x7b
mov WORD PTR [ebp+0x3b2a1f0e], ds
dec esp
pop ebp
outs dx, BYTE PTR ds:[esi]
jg 0xffffff9e
fwait
or al, 0x1d
cs aas
dec edx
pop ebx
ins BYTE PTR es:[edi], dx
jge 0xffffffac
lahf
.byte 0xa


3. New Value. AND AMS

9F0A1B2C3D4E5F6A7B8C9D0E1F2A3B4C5D6E0F8F9AFFFFFF9B0C1D2E3F4A5B6C0F8DA8FFFFFF9F0A

***RECURSE 1 ** RUN THIS CODE - Stable - Decompile this value
9F0A1B2C3D4E5F6A7B8C9D0E1F2A3B4C5D6E0F8F9AFFFFFF9B0C1D2E3F4A5B6C0F8DA8FFFFFF9F0A

        0:  9f                      lahf
        1:  0a 1b                   or     bl,BYTE PTR [ebx]
        3:  2c 3d                   sub    al,0x3d
        5:  4e                      dec    esi
        6:  5f                      pop    edi
        7:  6a 7b                   push   0x7b
        9:  8c 9d 0e 1f 2a 3b       mov    WORD PTR [ebp+0x3b2a1f0e],ds
        f:  4c                      dec    esp
        10: 5d                      pop    ebp
        11: 6e                      outs   dx,BYTE PTR ds:[esi]
        12: 0f 8f 9a ff ff ff       jg     ffffffb2 <_main+0xffffffb2>
        18: 9b                      fwait
        19: 0c 1d                   or     al,0x1d
        1b: 2e 3f                   cs aas
        1d: 4a                      dec    edx
        1e: 5b                      pop    ebx
        1f: 6c                      ins    BYTE PTR es:[edi],dx
        20: 0f 8d a8 ff ff ff       jge    ffffffce <_main+0xffffffce>
        26: 9f                      lahf
        27: 0a                      .byte 0xa



4. Decompile Above 

MED ANALYSIS
20 bit header 
9f0a1b2c3d4e5f6a7b8c9d0e1f2a3b4c5d6e          7f8a  9b0c   1d2e   3f4a   5b6c   7d8e   9f0a
9F0A1B2C3D4E5F6A7B8C9D0E1F2A3B4C5D6E          0F8F 9AFFFF FF9B0C 1D2E3F 4A5B6C 0F8DA8 FFFFFF 9F0A


7f8a  9b0c   1d2e   3f4a   5b6c   7d8e   9f0a
0F8F 9AFFFF FF9B0C 1D2E3F 4A5B6C 0F8DA8 FFFFFF 9F0A
       ?      ?      ?       ?      ?      ?


New hashed value
00 55 11 33 22 55 22 44 11 77 22 33 00 55
