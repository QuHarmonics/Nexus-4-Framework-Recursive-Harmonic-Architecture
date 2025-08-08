0:  42                      inc    edx
1:  8a 2f                   mov    ch,BYTE PTR [edi]
3:  98                      cwde
4:  71 37                   jno    0x3d
6:  44                      inc    esp
7:  91                      xchg   ecx,eax
8:  b5 c0                   mov    ch,0xc0
a:  fb                      sti
b:  cf                      iret
c:  e9 b5 db a5 39          jmp    0x39a5dbc6
11: 56                      push   esi
12: c2 5b 59                ret    0x595b
15: f1                      icebp
16: 11 f1                   adc    ecx,esi
18: 92                      xchg   edx,eax
19: 3f                      aas
1a: 82 a4 ab 1c 5e d5 d8    and    BYTE PTR [ebx+ebp*4-0x272aa1e4],0x7
21: 07
22: aa                      stos   BYTE PTR es:[edi],al
23: 98                      cwde
24: 12 83 5b 01 24 31       adc    al,BYTE PTR [ebx+0x3124015b]
2a: 85 be 55 0c 7d c3       test   DWORD PTR [esi-0x3c82f3ab],edi
30: 72 be                   jb     0xfffffff0
32: 5d                      pop    ebp
33: 74 80                   je     0xffffffb5
35: de b1 fe 9b dc 06       fidiv  WORD PTR [ecx+0x6dc9bfe]
3b: a7                      cmps   DWORD PTR ds:[esi],DWORD PTR es:[edi]
3c: c1 9b f1 74 e4 9b 69    rcr    DWORD PTR [ebx-0x641b8b0f],0x69
43: c1 ef be                shr    edi,0xbe
46: 47                      inc    edi
47: 86 0f                   xchg   BYTE PTR [edi],cl
49: c1 9d c6 24 0c a1 cc    rcr    DWORD PTR [ebp-0x5ef3db3a],0xcc
50: 2d e9 2c 6f 4a          sub    eax,0x4a6f2ce9
55: 74 84                   je     0xffffffdb
57: aa                      stos   BYTE PTR es:[edi],al
58: 5c                      pop    esp
59: b0 a9                   mov    al,0xa9
5b: dc 76 f9                fdiv   QWORD PTR [esi-0x7]
5e: 88 da                   mov    dl,bl
60: 98                      cwde
61: 3e 51                   ds push ecx
63: 52                      push   edx
64: a8 31                   test   al,0x31
66: c6                      (bad)
67: 6d                      ins    DWORD PTR es:[edi],dx
68: b0 03                   mov    al,0x3
6a: 27                      daa
6b: c8 bf 59 7f             enter  0x59bf,0x7f
6f: c7 c6 e0 0b f3 d5       mov    esi,0xd5f30be0
75: a7                      cmps   DWORD PTR ds:[esi],DWORD PTR es:[edi]
76: 91                      xchg   ecx,eax
77: 47                      inc    edi
78: 06                      push   es
79: ca 63 51                retf   0x5163
7c: 14 29                   adc    al,0x29
7e: 29 67 27                sub    DWORD PTR [edi+0x27],esp
81: b7 0a                   mov    bh,0xa
83: 85 2e                   test   DWORD PTR [esi],ebp
85: 1b 21                   sbb    esp,DWORD PTR [ecx]
87: 38 4d 2c                cmp    BYTE PTR [ebp+0x2c],cl
8a: 6d                      ins    DWORD PTR es:[edi],dx
8b: fc                      cld
8c: 53                      push   ebx
8d: 38 0d 13 65 0a 73       cmp    BYTE PTR ds:0x730a6513,cl
93: 54                      push   esp
94: 76 6a                   jbe    0x100
96: 0a bb 81 c2 c9 2e       or     bh,BYTE PTR [ebx+0x2ec9c281]
9c: 92                      xchg   edx,eax
9d: 72 2c                   jb     0xcb
9f: 85 a2 bf e8 a1 a8       test   DWORD PTR [edx-0x575e1741],esp
a5: 1a 66 4b                sbb    ah,BYTE PTR [esi+0x4b]
a8: c2 4b 8b                ret    0x8b4b
ab: 70 c7                   jo     0x74
ad: 6c                      ins    BYTE PTR es:[edi],dx
ae: 51                      push   ecx
af: a3 d1 92 e8 19          mov    ds:0x19e892d1,eax
b4: d6                      (bad)
b5: 99                      cdq
b6: 06                      push   es
b7: 24 f4                   and    al,0xf4
b9: 0e                      push   cs
ba: 35 85 10 6a a0          xor    eax,0xa06a1085
bf: 70 19                   jo     0xda
c1: a4                      movs   BYTE PTR es:[edi],BYTE PTR ds:[esi]
c2: c1 16 1e                rcl    DWORD PTR [esi],0x1e
c5: 37                      aaa
c6: 6c                      ins    BYTE PTR es:[edi],dx
c7: 08 27                   or     BYTE PTR [edi],ah
c9: 48                      dec    eax
ca: 77 4c                   ja     0x118
cc: 34 b0                   xor    al,0xb0
ce: bc b5 39 1c 0c          mov    esp,0xc1c39b5
d3: b3 4e                   mov    bl,0x4e
d5: d8 aa 4a 5b 9c ca       fsubr  DWORD PTR [edx-0x3563a4b6]
db: 4f                      dec    edi
dc: 68 2e 6f f3 74          push   0x74f36f2e
e1: 8f 82 ee 78 a5 63       pop    DWORD PTR [edx+0x63a578ee]
e7: 6f                      outs   dx,DWORD PTR ds:[esi]
e8: 84 c8                   test   al,cl
ea: 78 14                   js     0x100
ec: 8c c7                   mov    edi,es
ee: 02 08                   add    cl,BYTE PTR [eax]
f0: 90                      nop
f1: be ff fa a4 50          mov    esi,0x50a4faff
f6: 6c                      ins    BYTE PTR es:[edi],dx
f7: eb be                   jmp    0xb7
f9: f9                      stc
fa: a3 f7 c6 71 78          mov    ds:0x7871c6f7,eax
ff: f2                      repnz

import matplotlib.pyplot as plt

# Constants defined by the SHA-256 algorithm
K = [
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2
]

# Calculate differences and ratios between constants
differences = [K[i] - K[i - 1] for i in range(1, len(K))]
ratios = [K[i] / K[i - 1] for i in range(1, len(K))]

# Plot the constants, differences, and ratios
plt.figure(figsize=(15, 10))

# Plot constants
plt.subplot(3, 1, 1)
plt.plot(range(len(K)), K, label="Constants (K)", color="blue", marker="o")
plt.title("SHA-256 Constants (K)")
plt.xlabel("Index")
plt.ylabel("Value")
plt.grid(True)
plt.legend()

# Plot differences
plt.subplot(3, 1, 2)
plt.plot(range(1, len(K)), differences, label="Differences (K[i] - K[i-1])", color="red", marker="x")
plt.title("Differences Between Constants")
plt.xlabel("Index")
plt.ylabel("Difference")
plt.grid(True)
plt.legend()

# Plot ratios
plt.subplot(3, 1, 3)
plt.plot(range(1, len(K)), ratios, label="Ratios (K[i] / K[i-1])", color="green", marker="s")
plt.title("Ratios Between Constants")
plt.xlabel("Index")
plt.ylabel("Ratio")
plt.grid(True)
plt.legend()

plt.tight_layout()
plt.show()
