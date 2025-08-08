HARMONIE HIV URE
0:  33 33                   xor    esi,DWORD PTR [ebx]
2:  20 43 30                and    BYTE PTR [ebx+0x30],al
5:  0a 46 37                or     al,BYTE PTR [esi+0x37]
8:  20 45 32                and    BYTE PTR [ebp+0x32],al
b:  0a 43 37                or     al,BYTE PTR [ebx+0x37]
e:  20 34 35 20 46 43 20    and    BYTE PTR [esi*1+0x20434620],dh
15: 30 30                   xor    BYTE PTR [eax],dh
17: 20 30                   and    BYTE PTR [eax],dh
19: 30 20                   xor    BYTE PTR [eax],ah
1b: 30 30                   xor    BYTE PTR [eax],dh
1d: 20 30                   and    BYTE PTR [eax],dh
1f: 30 0a                   xor    BYTE PTR [edx],cl
21: 45                      inc    ebp
22: 45                      inc    ebp
23: 0a 45 42                or     al,BYTE PTR [ebp+0x42]
26: 20 46 45                and    BYTE PTR [esi+0x45],al
29: 0a 33                   or     dh,BYTE PTR [ebx]
2b: 31 20                   xor    DWORD PTR [eax],esp
2d: 44                      inc    esp
2e: 42                      inc    edx
2f: 0a 38                   or     bh,BYTE PTR [eax]
31: 42                      inc    edx
32: 20 43 33                and    BYTE PTR [ebx+0x33],al
35: 0a 32                   or     dh,BYTE PTR [edx]
37: 30 20                   xor    BYTE PTR [eax],ah
39: 32 30                   xor    dh,BYTE PTR [eax]
3b: 0a 37                   or     dh,BYTE PTR [edi]
3d: 35 20 46 36 0a          xor    eax,0xa364620