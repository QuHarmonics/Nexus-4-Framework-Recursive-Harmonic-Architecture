0:  30 3a                   xor    BYTE PTR [edx],bh
2:  20 20                   and    BYTE PTR [eax],ah
4:  30 30                   xor    BYTE PTR [eax],dh
6:  20 30                   and    BYTE PTR [eax],dh
8:  30 20                   xor    BYTE PTR [eax],ah
a:  20 20                   and    BYTE PTR [eax],ah
c:  20 20                   and    BYTE PTR [eax],ah
e:  20 20                   and    BYTE PTR [eax],ah
10: 20 20                   and    BYTE PTR [eax],ah
12: 20 20                   and    BYTE PTR [eax],ah
14: 20 20                   and    BYTE PTR [eax],ah
16: 20 20                   and    BYTE PTR [eax],ah
18: 20 20                   and    BYTE PTR [eax],ah
1a: 20 20                   and    BYTE PTR [eax],ah
1c: 61                      popa
1d: 64 64 20 20             fs and BYTE PTR fs:[eax],ah
21: 20 20                   and    BYTE PTR [eax],ah
23: 42                      inc    edx
24: 59                      pop    ecx
25: 54                      push   esp
26: 45                      inc    ebp
27: 20 50 54                and    BYTE PTR [eax+0x54],dl
2a: 52                      push   edx
2b: 20 5b 65                and    BYTE PTR [ebx+0x65],bl
2e: 61                      popa
2f: 78 5d                   js     0x8e
31: 2c 61                   sub    al,0x61
33: 6c                      ins    BYTE PTR es:[edi],dx
34: 0a 32                   or     dh,BYTE PTR [edx]
36: 3a 20                   cmp    ah,BYTE PTR [eax]
38: 20 30                   and    BYTE PTR [eax],dh
3a: 30 20                   xor    BYTE PTR [eax],ah
3c: 30 30                   xor    BYTE PTR [eax],dh
3e: 20 20                   and    BYTE PTR [eax],ah
40: 20 20                   and    BYTE PTR [eax],ah
42: 20 20                   and    BYTE PTR [eax],ah
44: 20 20                   and    BYTE PTR [eax],ah
46: 20 20                   and    BYTE PTR [eax],ah
48: 20 20                   and    BYTE PTR [eax],ah
4a: 20 20                   and    BYTE PTR [eax],ah
4c: 20 20                   and    BYTE PTR [eax],ah
4e: 20 20                   and    BYTE PTR [eax],ah
50: 20 61 64                and    BYTE PTR [ecx+0x64],ah
53: 64 20 20                and    BYTE PTR fs:[eax],ah
56: 20 20                   and    BYTE PTR [eax],ah
58: 42                      inc    edx
59: 59                      pop    ecx
5a: 54                      push   esp
5b: 45                      inc    ebp
5c: 20 50 54                and    BYTE PTR [eax+0x54],dl
5f: 52                      push   edx
60: 20 5b 65                and    BYTE PTR [ebx+0x65],bl
63: 61                      popa
64: 78 5d                   js     0xc3
66: 2c 61                   sub    al,0x61
68: 6c                      ins    BYTE PTR es:[edi],dx
69: 0a 34 3a                or     dh,BYTE PTR [edx+edi*1]
6c: 20 20                   and    BYTE PTR [eax],ah
6e: 66 66 20 20             data16 data16 and BYTE PTR [eax],ah
72: 20 20                   and    BYTE PTR [eax],ah
74: 20 20                   and    BYTE PTR [eax],ah
76: 20 20                   and    BYTE PTR [eax],ah
78: 20 20                   and    BYTE PTR [eax],ah
7a: 20 20                   and    BYTE PTR [eax],ah
7c: 20 20                   and    BYTE PTR [eax],ah
7e: 20 20                   and    BYTE PTR [eax],ah
80: 20 20                   and    BYTE PTR [eax],ah
82: 20 20                   and    BYTE PTR [eax],ah
84: 20 20                   and    BYTE PTR [eax],ah
86: 28 62 61                sub    BYTE PTR [edx+0x61],ah
89: 64 29 0a                sub    DWORD PTR fs:[edx],ecx
8c: 35 3a 20 20 66          xor    eax,0x6620203a
91: 66 20 30                data16 and BYTE PTR [eax],dh
94: 30 20                   xor    BYTE PTR [eax],ah
96: 20 20                   and    BYTE PTR [eax],ah
98: 20 20                   and    BYTE PTR [eax],ah
9a: 20 20                   and    BYTE PTR [eax],ah
9c: 20 20                   and    BYTE PTR [eax],ah
9e: 20 20                   and    BYTE PTR [eax],ah
a0: 20 20                   and    BYTE PTR [eax],ah
a2: 20 20                   and    BYTE PTR [eax],ah
a4: 20 20                   and    BYTE PTR [eax],ah
a6: 20 20                   and    BYTE PTR [eax],ah
a8: 69 6e 63 20 20 20 20    imul   ebp,DWORD PTR [esi+0x63],0x20202020
af: 44                      inc    esp
b0: 57                      push   edi
b1: 4f                      dec    edi
b2: 52                      push   edx
b3: 44                      inc    esp
b4: 20 50 54                and    BYTE PTR [eax+0x54],dl
b7: 52                      push   edx
b8: 20 5b 65                and    BYTE PTR [ebx+0x65],bl
bb: 61                      popa
bc: 78 5d                   js     0x11b
be: 0a 38                   or     bh,BYTE PTR [eax]
c0: 3a 20                   cmp    ah,BYTE PTR [eax]
c2: 20 30                   and    BYTE PTR [eax],dh
c4: 30 20                   xor    BYTE PTR [eax],ah
c6: 30 37                   xor    BYTE PTR [edi],dh
c8: 20 20                   and    BYTE PTR [eax],ah
ca: 20 20                   and    BYTE PTR [eax],ah
cc: 20 20                   and    BYTE PTR [eax],ah
ce: 20 20                   and    BYTE PTR [eax],ah
d0: 20 20                   and    BYTE PTR [eax],ah
d2: 20 20                   and    BYTE PTR [eax],ah
d4: 20 20                   and    BYTE PTR [eax],ah
d6: 20 20                   and    BYTE PTR [eax],ah
d8: 20 20                   and    BYTE PTR [eax],ah
da: 20 61 64                and    BYTE PTR [ecx+0x64],ah
dd: 64 20 20                and    BYTE PTR fs:[eax],ah
e0: 20 20                   and    BYTE PTR [eax],ah
e2: 42                      inc    edx
e3: 59                      pop    ecx
e4: 54                      push   esp
e5: 45                      inc    ebp
e6: 20 50 54                and    BYTE PTR [eax+0x54],dl
e9: 52                      push   edx
ea: 20 5b 65                and    BYTE PTR [ebx+0x65],bl
ed: 64 69 5d 2c 61 6c 0a    imul   ebx,DWORD PTR fs:[ebp+0x2c],0x610a6c61
f4: 61
f5: 3a 20                   cmp    ah,BYTE PTR [eax]
f7: 20 61 38                and    BYTE PTR [ecx+0x38],ah
fa: 20 62 38                and    BYTE PTR [edx+0x38],ah
fd: 20 20                   and    BYTE PTR [eax],ah
ff: 20 20                   and    BYTE PTR [eax],ah
101:    20 20                   and    BYTE PTR [eax],ah
103:    20 20                   and    BYTE PTR [eax],ah
105:    20 20                   and    BYTE PTR [eax],ah
107:    20 20                   and    BYTE PTR [eax],ah
109:    20 20                   and    BYTE PTR [eax],ah
10b:    20 20                   and    BYTE PTR [eax],ah
10d:    20 20                   and    BYTE PTR [eax],ah
10f:    20 74 65 73             and    BYTE PTR [ebp+eiz*2+0x73],dh
113:    74 20                   je     0x135
115:    20 20                   and    BYTE PTR [eax],ah
117:    61                      popa
118:    6c                      ins    BYTE PTR es:[edi],dx
119:    2c 30                   sub    al,0x30
11b:    78 62                   js     0x17f
11d:    38 0a                   cmp    BYTE PTR [edx],cl
11f:    63 3a                   arpl   WORD PTR [edx],di
121:    20 20                   and    BYTE PTR [eax],ah
123:    66 33 20                xor    sp,WORD PTR [eax]
126:    63 36                   arpl   WORD PTR [esi],si
128:    20 20                   and    BYTE PTR [eax],ah
12a:    20 20                   and    BYTE PTR [eax],ah
12c:    20 20                   and    BYTE PTR [eax],ah
12e:    20 20                   and    BYTE PTR [eax],ah
130:    20 20                   and    BYTE PTR [eax],ah
132:    20 20                   and    BYTE PTR [eax],ah
134:    20 20                   and    BYTE PTR [eax],ah
136:    20 20                   and    BYTE PTR [eax],ah
138:    20 20                   and    BYTE PTR [eax],ah
13a:    20 72 65                and    BYTE PTR [edx+0x65],dh
13d:    70 7a                   jo     0x1b9
13f:    20 28                   and    BYTE PTR [eax],ch
141:    62 61 64                bound  esp,QWORD PTR [ecx+0x64]
144:    29 0a                   sub    DWORD PTR [edx],ecx
146:    65 3a 20                cmp    ah,BYTE PTR gs:[eax]
149:    20 63 62                and    BYTE PTR [ebx+0x62],ah
14c:    20 20                   and    BYTE PTR [eax],ah
14e:    20 20                   and    BYTE PTR [eax],ah
150:    20 20                   and    BYTE PTR [eax],ah
152:    20 20                   and    BYTE PTR [eax],ah
154:    20 20                   and    BYTE PTR [eax],ah
156:    20 20                   and    BYTE PTR [eax],ah
158:    20 20                   and    BYTE PTR [eax],ah
15a:    20 20                   and    BYTE PTR [eax],ah
15c:    20 20                   and    BYTE PTR [eax],ah
15e:    20 20                   and    BYTE PTR [eax],ah
160:    20 20                   and    BYTE PTR [eax],ah
162:    72 65                   jb     0x1c9
164:    74 66                   je     0x1cc
166:    0a 66 3a                or     ah,BYTE PTR [esi+0x3a]
169:    20 20                   and    BYTE PTR [eax],ah
16b:    61                      popa
16c:    30 20                   xor    BYTE PTR [eax],ah
16e:    37                      aaa
16f:    62 20                   bound  esp,QWORD PTR [eax]
171:    31 64 20 39             xor    DWORD PTR [eax+eiz*1+0x39],esp
175:    66 20 34 64             data16 and BYTE PTR [esp+eiz*2],dh
179:    20 20                   and    BYTE PTR [eax],ah
17b:    20 20                   and    BYTE PTR [eax],ah
17d:    20 20                   and    BYTE PTR [eax],ah
17f:    20 20                   and    BYTE PTR [eax],ah
181:    20 20                   and    BYTE PTR [eax],ah
183:    6d                      ins    DWORD PTR es:[edi],dx
184:    6f                      outs   dx,DWORD PTR ds:[esi]
185:    76 20                   jbe    0x1a7
187:    20 20                   and    BYTE PTR [eax],ah
189:    20 61 6c                and    BYTE PTR [ecx+0x6c],ah
18c:    2c 64                   sub    al,0x64
18e:    73 3a                   jae    0x1ca
190:    30 78 34                xor    BYTE PTR [eax+0x34],bh
193:    64 39 66 31             cmp    DWORD PTR fs:[esi+0x31],esp
197:    64 37                   fs aaa
199:    62 0a                   bound  ecx,QWORD PTR [edx]
19b:    31 34 3a                xor    DWORD PTR [edx+edi*1],esi
19e:    20 37                   and    BYTE PTR [edi],dh
1a0:    38 20                   cmp    BYTE PTR [eax],ah
1a2:    61                      popa
1a3:    38 20                   cmp    BYTE PTR [eax],ah
1a5:    20 20                   and    BYTE PTR [eax],ah
1a7:    20 20                   and    BYTE PTR [eax],ah
1a9:    20 20                   and    BYTE PTR [eax],ah
1ab:    20 20                   and    BYTE PTR [eax],ah
1ad:    20 20                   and    BYTE PTR [eax],ah
1af:    20 20                   and    BYTE PTR [eax],ah
1b1:    20 20                   and    BYTE PTR [eax],ah
1b3:    20 20                   and    BYTE PTR [eax],ah
1b5:    20 20                   and    BYTE PTR [eax],ah
1b7:    6a 73                   push   0x73
1b9:    20 20                   and    BYTE PTR [eax],ah
1bb:    20 20                   and    BYTE PTR [eax],ah
1bd:    20 30                   and    BYTE PTR [eax],dh
1bf:    78 66                   js     0x227
1c1:    66 66 66 66 66 62 65    data16 data16 data16 data16 bound sp,DWORD PTR [ebp+0xa]
1c8:    0a
1c9:    31 36                   xor    DWORD PTR [esi],esi
1cb:    3a 20                   cmp    ah,BYTE PTR [eax]
1cd:    32 64 20 66             xor    ah,BYTE PTR [eax+eiz*1+0x66]
1d1:    64 20 33                and    BYTE PTR fs:[ebx],dh
1d4:    62 20                   bound  esp,QWORD PTR [eax]
1d6:    37                      aaa
1d7:    37                      aaa
1d8:    20 36                   and    BYTE PTR [esi],dh
1da:    36 20 20                and    BYTE PTR ss:[eax],ah
1dd:    20 20                   and    BYTE PTR [eax],ah
1df:    20 20                   and    BYTE PTR [eax],ah
1e1:    20 20                   and    BYTE PTR [eax],ah
1e3:    20 20                   and    BYTE PTR [eax],ah
1e5:    73 75                   jae    0x25c
1e7:    62 20                   bound  esp,QWORD PTR [eax]
1e9:    20 20                   and    BYTE PTR [eax],ah
1eb:    20 65 61                and    BYTE PTR [ebp+0x61],ah
1ee:    78 2c                   js     0x21c
1f0:    30 78 36                xor    BYTE PTR [eax+0x36],bh
1f3:    36 37                   ss aaa
1f5:    37                      aaa
1f6:    33 62 66                xor    esp,DWORD PTR [edx+0x66]
1f9:    64 0a 31                or     dh,BYTE PTR fs:[ecx]
1fc:    62 3a                   bound  edi,QWORD PTR [edx]
1fe:    20 65 31                and    BYTE PTR [ebp+0x31],ah
201:    20 39                   and    BYTE PTR [ecx],bh
203:    39 20                   cmp    DWORD PTR [eax],esp
205:    20 20                   and    BYTE PTR [eax],ah
207:    20 20                   and    BYTE PTR [eax],ah
209:    20 20                   and    BYTE PTR [eax],ah
20b:    20 20                   and    BYTE PTR [eax],ah
20d:    20 20                   and    BYTE PTR [eax],ah
20f:    20 20                   and    BYTE PTR [eax],ah
211:    20 20                   and    BYTE PTR [eax],ah
213:    20 20                   and    BYTE PTR [eax],ah
215:    20 20                   and    BYTE PTR [eax],ah
217:    6c                      ins    BYTE PTR es:[edi],dx
218:    6f                      outs   dx,DWORD PTR ds:[esi]
219:    6f                      outs   dx,DWORD PTR ds:[esi]
21a:    70 65                   jo     0x281
21c:    20 20                   and    BYTE PTR [eax],ah
21e:    30 78 66                xor    BYTE PTR [eax+0x66],bh
221:    66 66 66 66 66 62 36    data16 data16 data16 data16 bound si,DWORD PTR [esi]
228:    0a 31                   or     dh,BYTE PTR [ecx]
22a:    64 3a 20                cmp    ah,BYTE PTR fs:[eax]
22d:    64 32 20                xor    ah,BYTE PTR fs:[eax]
230:    20 20                   and    BYTE PTR [eax],ah
232:    20 20                   and    BYTE PTR [eax],ah
234:    20 20                   and    BYTE PTR [eax],ah
236:    20 20                   and    BYTE PTR [eax],ah
238:    20 20                   and    BYTE PTR [eax],ah
23a:    20 20                   and    BYTE PTR [eax],ah
23c:    20 20                   and    BYTE PTR [eax],ah
23e:    20 20                   and    BYTE PTR [eax],ah
240:    20 20                   and    BYTE PTR [eax],ah
242:    20 20                   and    BYTE PTR [eax],ah
244:    20 2e                   and    BYTE PTR [esi],ch
246:    62 79 74                bound  edi,QWORD PTR [ecx+0x74]
249:    65 20 30                and    BYTE PTR gs:[eax],dh
24c:    78 64                   js     0x2b2
24e:    32 0a                   xor    cl,BYTE PTR [edx]
250:    31 65 3a                xor    DWORD PTR [ebp+0x3a],esp
253:    20 61 34                and    BYTE PTR [ecx+0x34],ah
256:    20 20                   and    BYTE PTR [eax],ah
258:    20 20                   and    BYTE PTR [eax],ah
25a:    20 20                   and    BYTE PTR [eax],ah
25c:    20 20                   and    BYTE PTR [eax],ah
25e:    20 20                   and    BYTE PTR [eax],ah
260:    20 20                   and    BYTE PTR [eax],ah
262:    20 20                   and    BYTE PTR [eax],ah
264:    20 20                   and    BYTE PTR [eax],ah
266:    20 20                   and    BYTE PTR [eax],ah
268:    20 20                   and    BYTE PTR [eax],ah
26a:    20 20                   and    BYTE PTR [eax],ah
26c:    6d                      ins    DWORD PTR es:[edi],dx
26d:    6f                      outs   dx,DWORD PTR ds:[esi]
26e:    76 73                   jbe    0x2e3
270:    20 20                   and    BYTE PTR [eax],ah
272:    20 42 59                and    BYTE PTR [edx+0x59],al
275:    54                      push   esp
276:    45                      inc    ebp
277:    20 50 54                and    BYTE PTR [eax+0x54],dl
27a:    52                      push   edx
27b:    20 65 73                and    BYTE PTR [ebp+0x73],ah
27e:    3a 5b 65                cmp    bl,BYTE PTR [ebx+0x65]
281:    64 69 5d 2c 42 59 54    imul   ebx,DWORD PTR fs:[ebp+0x2c],0x45545942
288:    45
289:    20 50 54                and    BYTE PTR [eax+0x54],dl
28c:    52                      push   edx
28d:    20 64 73 3a             and    BYTE PTR [ebx+esi*2+0x3a],ah
291:    5b                      pop    ebx
292:    65 73 69                gs jae 0x2fe
295:    5d                      pop    ebp
296:    0a 31                   or     dh,BYTE PTR [ecx]
298:    66 3a 20                data16 cmp ah,BYTE PTR [eax]
29b:    33 31                   xor    esi,DWORD PTR [ecx]
29d:    20 20                   and    BYTE PTR [eax],ah
29f:    20 20                   and    BYTE PTR [eax],ah
2a1:    20 20                   and    BYTE PTR [eax],ah
2a3:    20 20                   and    BYTE PTR [eax],ah
2a5:    20 20                   and    BYTE PTR [eax],ah
2a7:    20 20                   and    BYTE PTR [eax],ah
2a9:    20 20                   and    BYTE PTR [eax],ah
2ab:    20 20                   and    BYTE PTR [eax],ah
2ad:    20 20                   and    BYTE PTR [eax],ah
2af:    20 20                   and    BYTE PTR [eax],ah
2b1:    20 20                   and    BYTE PTR [eax],ah
2b3:    2e 62 79 74             bound  edi,QWORD PTR cs:[ecx+0x74]
2b7:    65 20 30                and    BYTE PTR gs:[eax],dh
2ba:    78 33                   js     0x2ef
2bc:    31                      .byte 0x31