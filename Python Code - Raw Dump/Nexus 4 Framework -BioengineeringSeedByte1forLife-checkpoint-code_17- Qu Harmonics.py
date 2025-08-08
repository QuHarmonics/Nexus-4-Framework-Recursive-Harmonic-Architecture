0:  30 3a                   xor    BYTE PTR [edx],bh
2:  20 20                   and    BYTE PTR [eax],ah
4:  33 33                   xor    esi,DWORD PTR [ebx]
6:  20 63 30                and    BYTE PTR [ebx+0x30],ah
9:  20 20                   and    BYTE PTR [eax],ah
b:  20 20                   and    BYTE PTR [eax],ah
d:  20 20                   and    BYTE PTR [eax],ah
f:  20 20                   and    BYTE PTR [eax],ah
11: 20 20                   and    BYTE PTR [eax],ah
13: 20 20                   and    BYTE PTR [eax],ah
15: 20 20                   and    BYTE PTR [eax],ah
17: 20 20                   and    BYTE PTR [eax],ah
19: 20 20                   and    BYTE PTR [eax],ah
1b: 20 78 6f                and    BYTE PTR [eax+0x6f],bh
1e: 72 20                   jb     0x40
20: 20 20                   and    BYTE PTR [eax],ah
22: 20 65 61                and    BYTE PTR [ebp+0x61],ah
25: 78 2c                   js     0x53
27: 65 61                   gs popa
29: 78 0a                   js     0x35
2b: 32 3a                   xor    bh,BYTE PTR [edx]
2d: 20 20                   and    BYTE PTR [eax],ah
2f: 66 37                   data16 aaa
31: 20 65 32                and    BYTE PTR [ebp+0x32],ah
34: 20 20                   and    BYTE PTR [eax],ah
36: 20 20                   and    BYTE PTR [eax],ah
38: 20 20                   and    BYTE PTR [eax],ah
3a: 20 20                   and    BYTE PTR [eax],ah
3c: 20 20                   and    BYTE PTR [eax],ah
3e: 20 20                   and    BYTE PTR [eax],ah
40: 20 20                   and    BYTE PTR [eax],ah
42: 20 20                   and    BYTE PTR [eax],ah
44: 20 20                   and    BYTE PTR [eax],ah
46: 20 6d 75                and    BYTE PTR [ebp+0x75],ch
49: 6c                      ins    BYTE PTR es:[edi],dx
4a: 20 20                   and    BYTE PTR [eax],ah
4c: 20 20                   and    BYTE PTR [eax],ah
4e: 65 64 78 0a             gs fs js 0x5c
52: 34 3a                   xor    al,0x3a
54: 20 20                   and    BYTE PTR [eax],ah
56: 63 37                   arpl   WORD PTR [edi],si
58: 20 34 35 20 66 63 20    and    BYTE PTR [esi*1+0x20636620],dh
5f: 30 30                   xor    BYTE PTR [eax],dh
61: 20 30                   and    BYTE PTR [eax],dh
63: 30 20                   xor    BYTE PTR [eax],ah
65: 30 30                   xor    BYTE PTR [eax],dh
67: 20 30                   and    BYTE PTR [eax],dh
69: 30 20                   xor    BYTE PTR [eax],ah
6b: 20 20                   and    BYTE PTR [eax],ah
6d: 20 6d 6f                and    BYTE PTR [ebp+0x6f],ch
70: 76 20                   jbe    0x92
72: 20 20                   and    BYTE PTR [eax],ah
74: 20 44 57 4f             and    BYTE PTR [edi+edx*2+0x4f],al
78: 52                      push   edx
79: 44                      inc    esp
7a: 20 50 54                and    BYTE PTR [eax+0x54],dl
7d: 52                      push   edx
7e: 20 5b 65                and    BYTE PTR [ebx+0x65],bl
81: 62 70 2d                bound  esi,QWORD PTR [eax+0x2d]
84: 30 78 34                xor    BYTE PTR [eax+0x34],bh
87: 5d                      pop    ebp
88: 2c 30                   sub    al,0x30
8a: 78 30                   js     0xbc
8c: 0a 62 3a                or     ah,BYTE PTR [edx+0x3a]
8f: 20 20                   and    BYTE PTR [eax],ah
91: 65 65 20 20             gs and BYTE PTR gs:[eax],ah
95: 20 20                   and    BYTE PTR [eax],ah
97: 20 20                   and    BYTE PTR [eax],ah
99: 20 20                   and    BYTE PTR [eax],ah
9b: 20 20                   and    BYTE PTR [eax],ah
9d: 20 20                   and    BYTE PTR [eax],ah
9f: 20 20                   and    BYTE PTR [eax],ah
a1: 20 20                   and    BYTE PTR [eax],ah
a3: 20 20                   and    BYTE PTR [eax],ah
a5: 20 20                   and    BYTE PTR [eax],ah
a7: 20 20                   and    BYTE PTR [eax],ah
a9: 6f                      outs   dx,DWORD PTR ds:[esi]
aa: 75 74                   jne    0x120
ac: 20 20                   and    BYTE PTR [eax],ah
ae: 20 20                   and    BYTE PTR [eax],ah
b0: 64 78 2c                fs js  0xdf
b3: 61                      popa
b4: 6c                      ins    BYTE PTR es:[edi],dx
b5: 0a 63 3a                or     ah,BYTE PTR [ebx+0x3a]
b8: 20 20                   and    BYTE PTR [eax],ah
ba: 65 62 20                bound  esp,QWORD PTR gs:[eax]
bd: 66 65 20 20             data16 and BYTE PTR gs:[eax],ah
c1: 20 20                   and    BYTE PTR [eax],ah
c3: 20 20                   and    BYTE PTR [eax],ah
c5: 20 20                   and    BYTE PTR [eax],ah
c7: 20 20                   and    BYTE PTR [eax],ah
c9: 20 20                   and    BYTE PTR [eax],ah
cb: 20 20                   and    BYTE PTR [eax],ah
cd: 20 20                   and    BYTE PTR [eax],ah
cf: 20 20                   and    BYTE PTR [eax],ah
d1: 20 6a 6d                and    BYTE PTR [edx+0x6d],ch
d4: 70 20                   jo     0xf6
d6: 20 20                   and    BYTE PTR [eax],ah
d8: 20 30                   and    BYTE PTR [eax],dh
da: 78 63                   js     0x13f
dc: 0a 65 3a                or     ah,BYTE PTR [ebp+0x3a]
df: 20 20                   and    BYTE PTR [eax],ah
e1: 33 31                   xor    esi,DWORD PTR [ecx]
e3: 20 64 62 20             and    BYTE PTR [edx+eiz*2+0x20],ah
e7: 20 20                   and    BYTE PTR [eax],ah
e9: 20 20                   and    BYTE PTR [eax],ah
eb: 20 20                   and    BYTE PTR [eax],ah
ed: 20 20                   and    BYTE PTR [eax],ah
ef: 20 20                   and    BYTE PTR [eax],ah
f1: 20 20                   and    BYTE PTR [eax],ah
f3: 20 20                   and    BYTE PTR [eax],ah
f5: 20 20                   and    BYTE PTR [eax],ah
f7: 20 20                   and    BYTE PTR [eax],ah
f9: 78 6f                   js     0x16a
fb: 72 20                   jb     0x11d
fd: 20 20                   and    BYTE PTR [eax],ah
ff: 20 65 62                and    BYTE PTR [ebp+0x62],ah
102:    78 2c                   js     0x130
104:    65 62 78 0a             bound  edi,QWORD PTR gs:[eax+0xa]
108:    31 30                   xor    DWORD PTR [eax],esi
10a:    3a 20                   cmp    ah,BYTE PTR [eax]
10c:    38 62 20                cmp    BYTE PTR [edx+0x20],ah
10f:    63 33                   arpl   WORD PTR [ebx],si
111:    20 20                   and    BYTE PTR [eax],ah
113:    20 20                   and    BYTE PTR [eax],ah
115:    20 20                   and    BYTE PTR [eax],ah
117:    20 20                   and    BYTE PTR [eax],ah
119:    20 20                   and    BYTE PTR [eax],ah
11b:    20 20                   and    BYTE PTR [eax],ah
11d:    20 20                   and    BYTE PTR [eax],ah
11f:    20 20                   and    BYTE PTR [eax],ah
121:    20 20                   and    BYTE PTR [eax],ah
123:    20 6d 6f                and    BYTE PTR [ebp+0x6f],ch
126:    76 20                   jbe    0x148
128:    20 20                   and    BYTE PTR [eax],ah
12a:    20 65 61                and    BYTE PTR [ebp+0x61],ah
12d:    78 2c                   js     0x15b
12f:    65 62 78 0a             bound  edi,QWORD PTR gs:[eax+0xa]
133:    31 32                   xor    DWORD PTR [edx],esi
135:    3a 20                   cmp    ah,BYTE PTR [eax]
137:    32 30                   xor    dh,BYTE PTR [eax]
139:    20 32                   and    BYTE PTR [edx],dh
13b:    30 20                   xor    BYTE PTR [eax],ah
13d:    20 20                   and    BYTE PTR [eax],ah
13f:    20 20                   and    BYTE PTR [eax],ah
141:    20 20                   and    BYTE PTR [eax],ah
143:    20 20                   and    BYTE PTR [eax],ah
145:    20 20                   and    BYTE PTR [eax],ah
147:    20 20                   and    BYTE PTR [eax],ah
149:    20 20                   and    BYTE PTR [eax],ah
14b:    20 20                   and    BYTE PTR [eax],ah
14d:    20 20                   and    BYTE PTR [eax],ah
14f:    61                      popa
150:    6e                      outs   dx,BYTE PTR ds:[esi]
151:    64 20 20                and    BYTE PTR fs:[eax],ah
154:    20 20                   and    BYTE PTR [eax],ah
156:    42                      inc    edx
157:    59                      pop    ecx
158:    54                      push   esp
159:    45                      inc    ebp
15a:    20 50 54                and    BYTE PTR [eax+0x54],dl
15d:    52                      push   edx
15e:    20 5b 65                and    BYTE PTR [ebx+0x65],bl
161:    61                      popa
162:    78 5d                   js     0x1c1
164:    2c 61                   sub    al,0x61
166:    68 0a 31 34 3a          push   0x3a34310a
16b:    20 37                   and    BYTE PTR [edi],dh
16d:    35 20 66 36 20          xor    eax,0x20366620
172:    20 20                   and    BYTE PTR [eax],ah
174:    20 20                   and    BYTE PTR [eax],ah
176:    20 20                   and    BYTE PTR [eax],ah
178:    20 20                   and    BYTE PTR [eax],ah
17a:    20 20                   and    BYTE PTR [eax],ah
17c:    20 20                   and    BYTE PTR [eax],ah
17e:    20 20                   and    BYTE PTR [eax],ah
180:    20 20                   and    BYTE PTR [eax],ah
182:    20 20                   and    BYTE PTR [eax],ah
184:    6a 6e                   push   0x6e
186:    65 20 20                and    BYTE PTR gs:[eax],ah
189:    20 20                   and    BYTE PTR [eax],ah
18b:    30 78 63                xor    BYTE PTR [eax+0x63],bh