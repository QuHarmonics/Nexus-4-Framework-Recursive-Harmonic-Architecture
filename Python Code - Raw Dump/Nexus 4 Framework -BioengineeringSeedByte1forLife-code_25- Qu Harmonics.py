0:  33 33                   xor    esi,DWORD PTR [ebx]
2:  20 33                   and    BYTE PTR [ebx],dh
4:  33 32                   xor    esi,DWORD PTR [edx]
6:  30 20                   xor    BYTE PTR [eax],ah
8:  33 33                   xor    esi,DWORD PTR [ebx]
a:  33 33                   xor    esi,DWORD PTR [ebx]
c:  20 30                   and    BYTE PTR [eax],dh
e:  41                      inc    ecx
f:  33 32                   xor    esi,DWORD PTR [edx]
11: 20 33                   and    BYTE PTR [ebx],dh
13: 30 32                   xor    BYTE PTR [edx],dh
15: 30 20                   xor    BYTE PTR [eax],ah
17: 33 34 20                xor    esi,DWORD PTR [eax+eiz*1]
1a: 33 33                   xor    esi,DWORD PTR [ebx]
1c: 32 30                   xor    dh,BYTE PTR [eax]
1e: 20 33                   and    BYTE PTR [ebx],dh
20: 33 33                   xor    esi,DWORD PTR [ebx]
22: 30 20                   xor    BYTE PTR [eax],ah
24: 30 41 33                xor    BYTE PTR [ecx+0x33],al
27: 30 20                   xor    BYTE PTR [eax],ah
29: 34 31                   xor    al,0x31
2b: 20 32                   and    BYTE PTR [edx],dh
2d: 30 33                   xor    BYTE PTR [ebx],dh
2f: 34 20                   xor    al,0x20
31: 33 36                   xor    esi,DWORD PTR [esi]
33: 32 30                   xor    dh,BYTE PTR [eax]
35: 20 33                   and    BYTE PTR [ebx],dh
37: 33 33                   xor    esi,DWORD PTR [ebx]
39: 37                      aaa
3a: 30 41 20                xor    BYTE PTR [ecx+0x20],al
3d: 33 32                   xor    esi,DWORD PTR [edx]
3f: 33 30                   xor    esi,DWORD PTR [eax]
41: 20 32                   and    BYTE PTR [edx],dh
43: 30 33                   xor    BYTE PTR [ebx],dh
45: 34 20                   xor    al,0x20
47: 33 35 32 30 20 33       xor    esi,DWORD PTR ds:0x33203032
4d: 33 33                   xor    esi,DWORD PTR [ebx]
4f: 32 20                   xor    ah,BYTE PTR [eax]
51: 30 41 33                xor    BYTE PTR [ecx+0x33],al
54: 30 20                   xor    BYTE PTR [eax],ah
56: 34 31                   xor    al,0x31
58: 20 32                   and    BYTE PTR [edx],dh
5a: 30 33                   xor    BYTE PTR [ebx],dh
5c: 34 20                   xor    al,0x20
5e: 33 33                   xor    esi,DWORD PTR [ebx]
60: 32 30                   xor    dh,BYTE PTR [eax]
62: 20 33                   and    BYTE PTR [ebx],dh
64: 33 33                   xor    esi,DWORD PTR [ebx]
66: 37                      aaa
67: 30 41 20                xor    BYTE PTR [ecx+0x20],al
6a: 33 32                   xor    esi,DWORD PTR [edx]
6c: 33 30                   xor    esi,DWORD PTR [eax]
6e: 20 32                   and    BYTE PTR [edx],dh
70: 30 33                   xor    BYTE PTR [ebx],dh
72: 33 20                   xor    esp,DWORD PTR [eax]
74: 33 34 20                xor    esi,DWORD PTR [eax+eiz*1]
77: 32 30                   xor    dh,BYTE PTR [eax]
79: 33 33                   xor    esi,DWORD PTR [ebx]
7b: 20 33                   and    BYTE PTR [ebx],dh
7d: 35 20 32 30 20          xor    eax,0x20303220
82: 33 32                   xor    esi,DWORD PTR [edx]
84: 20 33                   and    BYTE PTR [ebx],dh
86: 30 20                   xor    BYTE PTR [eax],ah
88: 32 30                   xor    dh,BYTE PTR [eax]
8a: 33 34 20                xor    esi,DWORD PTR [eax+eiz*1]
8d: 33 36                   xor    esi,DWORD PTR [esi]
8f: 32 30                   xor    dh,BYTE PTR [eax]
91: 20 33                   and    BYTE PTR [ebx],dh
93: 34 20                   xor    al,0x20
95: 33 33                   xor    esi,DWORD PTR [ebx]
97: 32 30                   xor    dh,BYTE PTR [eax]
99: 20 33                   and    BYTE PTR [ebx],dh
9b: 32 33                   xor    dh,BYTE PTR [ebx]
9d: 30 20                   xor    BYTE PTR [eax],ah
9f: 30 41 33                xor    BYTE PTR [ecx+0x33],al
a2: 33 20                   xor    esp,DWORD PTR [eax]
a4: 33 30                   xor    esi,DWORD PTR [eax]
a6: 32 30                   xor    dh,BYTE PTR [eax]
a8: 20 33                   and    BYTE PTR [ebx],dh
aa: 33 33                   xor    esi,DWORD PTR [ebx]
ac: 30 20                   xor    BYTE PTR [eax],ah
ae: 30 41 33                xor    BYTE PTR [ecx+0x33],al
b1: 32 20                   xor    ah,BYTE PTR [eax]
b3: 33 30                   xor    esi,DWORD PTR [eax]
b5: 32 30                   xor    dh,BYTE PTR [eax]
b7: 20 33                   and    BYTE PTR [ebx],dh
b9: 33 33                   xor    esi,DWORD PTR [ebx]
bb: 30 20                   xor    BYTE PTR [eax],ah
bd: 30 41 33                xor    BYTE PTR [ecx+0x33],al
c0: 33 20                   xor    esp,DWORD PTR [eax]
c2: 33 30                   xor    esi,DWORD PTR [eax]
c4: 32 30                   xor    dh,BYTE PTR [eax]
c6: 20 33                   and    BYTE PTR [ebx],dh
c8: 32 33                   xor    dh,BYTE PTR [ebx]
ca: 30 20                   xor    BYTE PTR [eax],ah
cc: 30 41 33                xor    BYTE PTR [ecx+0x33],al
cf: 33 20                   xor    esp,DWORD PTR [eax]
d1: 33 30                   xor    esi,DWORD PTR [eax]
d3: 32 30                   xor    dh,BYTE PTR [eax]
d5: 20 33                   and    BYTE PTR [ebx],dh
d7: 33 33                   xor    esi,DWORD PTR [ebx]
d9: 30 20                   xor    BYTE PTR [eax],ah
db: 30 41 33                xor    BYTE PTR [ecx+0x33],al
de: 32 20                   xor    ah,BYTE PTR [eax]
e0: 33 30                   xor    esi,DWORD PTR [eax]
e2: 32 30                   xor    dh,BYTE PTR [eax]
e4: 20 33                   and    BYTE PTR [ebx],dh
e6: 33 33                   xor    esi,DWORD PTR [ebx]
e8: 30 20                   xor    BYTE PTR [eax],ah
ea: 30 41 33                xor    BYTE PTR [ecx+0x33],al
ed: 33 20                   xor    esp,DWORD PTR [eax]
ef: 33 30                   xor    esi,DWORD PTR [eax]
f1: 32 30                   xor    dh,BYTE PTR [eax]
f3: 20 33                   and    BYTE PTR [ebx],dh
f5: 30 34 31                xor    BYTE PTR [ecx+esi*1],dh
f8: 30 41 20                xor    BYTE PTR [ecx+0x20],al
fb: 33 34 20                xor    esi,DWORD PTR [eax+eiz*1]
fe: 33 35 20 30 41 20       xor    esi,DWORD PTR ds:0x20413020
104:    33 34 20                xor    esi,DWORD PTR [eax+eiz*1]
107:    33 35 20 30 41 33       xor    esi,DWORD PTR ds:0x33413020
10d:    30 20                   xor    BYTE PTR [eax],ah
10f:    34 31                   xor    al,0x31
111:    20 32                   and    BYTE PTR [edx],dh
113:    30 33                   xor    BYTE PTR [ebx],dh
115:    34 20                   xor    al,0x20
117:    33 35 32 30 20 33       xor    esi,DWORD PTR ds:0x33203032
11d:    34 20                   xor    al,0x20
11f:    33 32                   xor    esi,DWORD PTR [edx]
121:    30 41 20                xor    BYTE PTR [ecx+0x20],al
124:    33 32                   xor    esi,DWORD PTR [edx]
126:    33 30                   xor    esi,DWORD PTR [eax]
128:    20 32                   and    BYTE PTR [edx],dh
12a:    30 33                   xor    BYTE PTR [ebx],dh
12c:    34 20                   xor    al,0x20
12e:    33 36                   xor    esi,DWORD PTR [esi]
130:    32 30                   xor    dh,BYTE PTR [eax]
132:    20 33                   and    BYTE PTR [ebx],dh
134:    34 20                   xor    al,0x20
136:    33 35 20 30 41 20       xor    esi,DWORD PTR ds:0x20413020
13c:    33 30                   xor    esi,DWORD PTR [eax]
13e:    20 34 31                and    BYTE PTR [ecx+esi*1],dh
141:    20 32                   and    BYTE PTR [edx],dh
143:    30 33                   xor    BYTE PTR [ebx],dh
145:    33 20                   xor    esp,DWORD PTR [eax]
147:    33 33                   xor    esi,DWORD PTR [ebx]
149:    30 41 20                xor    BYTE PTR [ecx+0x20],al
14c:    33 33                   xor    esi,DWORD PTR [ebx]
14e:    33 31                   xor    esi,DWORD PTR [ecx]
150:    20 32                   and    BYTE PTR [edx],dh
152:    30 33                   xor    BYTE PTR [ebx],dh
154:    32 20                   xor    ah,BYTE PTR [eax]
156:    33 30                   xor    esi,DWORD PTR [eax]
158:    30 41 20                xor    BYTE PTR [ecx+0x20],al
15b:    33 34 20                xor    esi,DWORD PTR [eax+eiz*1]
15e:    33 34 30                xor    esi,DWORD PTR [eax+esi*1]
161:    41                      inc    ecx
162:    20 33                   and    BYTE PTR [ebx],dh
164:    34 20                   xor    al,0x20
166:    33 32                   xor    esi,DWORD PTR [edx]
168:    30 41 20                xor    BYTE PTR [ecx+0x20],al
16b:    33 30                   xor    esi,DWORD PTR [eax]
16d:    34 31                   xor    al,0x31
16f:    32 30                   xor    dh,BYTE PTR [eax]
171:    20 33                   and    BYTE PTR [ebx],dh
173:    33 33                   xor    esi,DWORD PTR [ebx]
175:    38 20                   cmp    BYTE PTR [eax],ah
177:    30 41 33                xor    BYTE PTR [ecx+0x33],al
17a:    34 20                   xor    al,0x20
17c:    33 32                   xor    esi,DWORD PTR [edx]
17e:    30 41 20                xor    BYTE PTR [ecx+0x20],al
181:    33 32                   xor    esi,DWORD PTR [edx]
183:    33 30                   xor    esi,DWORD PTR [eax]
185:    20 32                   and    BYTE PTR [edx],dh
187:    30 33                   xor    BYTE PTR [ebx],dh
189:    34 20                   xor    al,0x20
18b:    33 33                   xor    esi,DWORD PTR [ebx]
18d:    32 30                   xor    dh,BYTE PTR [eax]
18f:    20 33                   and    BYTE PTR [ebx],dh
191:    33 33                   xor    esi,DWORD PTR [ebx]
193:    33 20                   xor    esp,DWORD PTR [eax]
195:    30 41 33                xor    BYTE PTR [ecx+0x33],al
198:    30 20                   xor    BYTE PTR [eax],ah
19a:    34 31                   xor    al,0x31
19c:    20 32                   and    BYTE PTR [edx],dh
19e:    30 33                   xor    BYTE PTR [ebx],dh
1a0:    33 20                   xor    esp,DWORD PTR [eax]
1a2:    33 32                   xor    esi,DWORD PTR [edx]
1a4:    30 41 20                xor    BYTE PTR [ecx+0x20],al
1a7:    33 33                   xor    esi,DWORD PTR [ebx]
1a9:    33 30                   xor    esi,DWORD PTR [eax]
1ab:    20 32                   and    BYTE PTR [edx],dh
1ad:    30 33                   xor    BYTE PTR [ebx],dh
1af:    32 20                   xor    ah,BYTE PTR [eax]
1b1:    33 30                   xor    esi,DWORD PTR [eax]
1b3:    30 41 20                xor    BYTE PTR [ecx+0x20],al
1b6:    33 33                   xor    esi,DWORD PTR [ebx]
1b8:    33 32                   xor    esi,DWORD PTR [edx]
1ba:    20 32                   and    BYTE PTR [edx],dh
1bc:    30 33                   xor    BYTE PTR [ebx],dh
1be:    33 20                   xor    esp,DWORD PTR [eax]
1c0:    33 30                   xor    esi,DWORD PTR [eax]
1c2:    30 41 20                xor    BYTE PTR [ecx+0x20],al
1c5:    33 30                   xor    esi,DWORD PTR [eax]
1c7:    34 31                   xor    al,0x31
1c9:    32 30                   xor    dh,BYTE PTR [eax]
1cb:    20 33                   and    BYTE PTR [ebx],dh
1cd:    33 33                   xor    esi,DWORD PTR [ebx]
1cf:    37                      aaa
1d0:    30 41 20                xor    BYTE PTR [ecx+0x20],al
1d3:    33 33                   xor    esi,DWORD PTR [ebx]
1d5:    33 35 20 32 30 20       xor    esi,DWORD PTR ds:0x20303220
1db:    33 32                   xor    esi,DWORD PTR [edx]
1dd:    20 33                   and    BYTE PTR [ebx],dh
1df:    30 20                   xor    BYTE PTR [eax],ah
1e1:    32 30                   xor    dh,BYTE PTR [eax]
1e3:    33 34 20                xor    esi,DWORD PTR [eax+eiz*1]
1e6:    33 36                   xor    esi,DWORD PTR [esi]
1e8:    32 30                   xor    dh,BYTE PTR [eax]
1ea:    20 33                   and    BYTE PTR [ebx],dh
1ec:    33 33                   xor    esi,DWORD PTR [ebx]
1ee:    36 20 32                and    BYTE PTR ss:[edx],dh
1f1:    30 20                   xor    BYTE PTR [eax],ah
1f3:    33 30                   xor    esi,DWORD PTR [eax]
1f5:    34 31                   xor    al,0x31
1f7:    30                      .byte 0x30
1f8:    41                      inc    ecx