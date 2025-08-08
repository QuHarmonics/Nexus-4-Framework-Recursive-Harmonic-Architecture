0:  61                      popa
1:  63 74 63 74             arpl   WORD PTR [ebx+eiz*2+0x74],si
5:  74 74                   je     0x7b
7:  67 67 63 20             addr16 arpl WORD PTR [bx+si],sp
b:  61                      popa
c:  61                      popa
d:  63 67 61                arpl   WORD PTR [edi+0x61],sp
10: 63 63 63                arpl   WORD PTR [ebx+0x63],sp
13: 63 74 20 63             arpl   WORD PTR [eax+eiz*1+0x63],si
17: 67 74 63                addr16 je 0x7d
1a: 61                      popa
1b: 63 61 61                arpl   WORD PTR [ecx+0x61],sp
1e: 74 61                   je     0x81
20: 20 61 61                and    BYTE PTR [ecx+0x61],ah
23: 67 61                   addr16 popa
25: 74 61                   je     0x88
27: 67 67 67 67 20 67 67    addr16 addr16 addr16 and BYTE PTR [bx+0x67],ah
2e: 63 61 61                arpl   WORD PTR [ecx+0x61],sp
31: 63 74 61 61             arpl   WORD PTR [ecx+eiz*2+0x61],si
35: 61                      popa
36: 20 67 67                and    BYTE PTR [edi+0x67],ah
39: 61                      popa
3a: 61                      popa
3b: 67 63 74 63             arpl   WORD PTR [si+0x63],si
3f: 74 61                   je     0xa2
41: 74 74                   je     0xb7
43: 61                      popa
44: 67 61                   addr16 popa
46: 74 61                   je     0xa9
48: 63 61 67                arpl   WORD PTR [ecx+0x67],sp
4b: 20 67 61                and    BYTE PTR [edi+0x61],ah
4e: 67 63 61 67             arpl   WORD PTR [bx+di+0x67],sp
52: 61                      popa
53: 74 67                   je     0xbc
55: 61                      popa
56: 20 74 61 63             and    BYTE PTR [ecx+eiz*2+0x63],dh
5a: 61                      popa
5b: 67 74 61                addr16 je 0xbf
5e: 74 74                   je     0xd4
60: 61                      popa
61: 20 67 61                and    BYTE PTR [edi+0x61],ah
64: 61                      popa
65: 67 61                   addr16 popa
67: 61                      popa
68: 61                      popa
69: 74 67                   je     0xd2
6b: 61                      popa
6c: 20 67 74                and    BYTE PTR [edi+0x74],ah
6f: 74 74                   je     0xe5
71: 67 63 63 61             arpl   WORD PTR [bp+di+0x61],sp
75: 67 67 20 61 61          addr16 and BYTE PTR [bx+di+0x61],ah
7a: 67 61                   addr16 popa
7c: 74 67                   je     0xe5
7e: 67 61                   addr16 popa
80: 61                      popa
81: 61                      popa
82: 63 63 61                arpl   WORD PTR [ebx+0x61],sp
85: 61                      popa
86: 61                      popa
87: 61                      popa
88: 61                      popa
89: 74 67                   je     0xf2
8b: 61                      popa
8c: 20 74 61 67             and    BYTE PTR [ecx+eiz*2+0x67],dh
90: 67 67 67 67 61          addr16 addr16 addr16 addr16 popa
95: 61                      popa
96: 74 20                   je     0xb8
98: 74 67                   je     0x101
9a: 67 61                   addr16 popa
9c: 67 67 74 74             addr16 addr16 je 0x114
a0: 74 74                   je     0x116
a2: 20 61 74                and    BYTE PTR [ecx+0x74],ah
a5: 63 61 61                arpl   WORD PTR [ecx+0x61],sp
a8: 61                      popa
a9: 67 74 61                addr16 je 0x10d
ac: 61                      popa
ad: 20 67 61                and    BYTE PTR [edi+0x61],ah
b0: 63 61 67                arpl   WORD PTR [ecx+0x67],sp
b3: 74 61                   je     0x116
b5: 74 67                   je     0x11e
b7: 61                      popa
b8: 20 74 63 61             and    BYTE PTR [ebx+eiz*2+0x61],dh
bc: 67 61                   addr16 popa
be: 74 61                   je     0x121
c0: 63 74 63 61             arpl   WORD PTR [ebx+eiz*2+0x61],si
c4: 74 61                   je     0x127
c6: 67 61                   addr16 popa
c8: 61                      popa
c9: 61                      popa
ca: 74 63                   je     0x12f
cc: 74 20                   je     0xee
ce: 67 74 67                addr16 je 0x138
d1: 67 61                   addr16 popa
d3: 63 61 74                arpl   WORD PTR [ecx+0x74],sp
d6: 61                      popa
d7: 61                      popa
d8: 20 61 67                and    BYTE PTR [ecx+0x67],ah
db: 63 74 61 74             arpl   WORD PTR [ecx+eiz*2+0x74],si
df: 61                      popa
e0: 67 67 74 20             addr16 addr16 je 0x104
e4: 61                      popa
e5: 63 61 67                arpl   WORD PTR [ecx+0x67],sp
e8: 74 61                   je     0x14b
ea: 74 74                   je     0x160
ec: 61                      popa
ed: 67 20 74 61             and    BYTE PTR [si+0x61],dh
f1: 67 67 61                addr16 addr16 popa
f4: 63 63 74                arpl   WORD PTR [ebx+0x74],sp
f7: 61                      popa
f8: 63 20                   arpl   WORD PTR [eax],sp
fa: 61                      popa
fb: 63 63 74                arpl   WORD PTR [ebx+0x74],sp
fe: 67 74 63                addr16 je 0x164
101:    61                      popa
102:    61                      popa
103:    63 61 74                arpl   WORD PTR [ecx+0x74],sp
106:    61                      popa
107:    61                      popa
108:    74 74                   je     0x17e
10a:    67 67 61                addr16 addr16 popa
10d:    61                      popa
10e:    20 67 61                and    BYTE PTR [edi+0x61],ah
111:    61                      popa
112:    61                      popa
113:    74 63                   je     0x178
115:    74 67                   je     0x17e
117:    74 74                   je     0x18d
119:    20 67 61                and    BYTE PTR [edi+0x61],ah
11c:    63 74 63 61             arpl   WORD PTR [ebx+eiz*2+0x61],si
120:    67 61                   addr16 popa
122:    74 74                   je     0x198
124:    20 67 67                and    BYTE PTR [edi+0x67],ah
127:    74 74                   je     0x19d
129:    67 63 61 63             arpl   WORD PTR [bx+di+0x63],sp
12d:    74 74                   je     0x1a3
12f:    20 74 61 61             and    BYTE PTR [ecx+eiz*2+0x61],dh
133:    61                      popa
134:    74 74                   je     0x1aa
136:    74 74                   je     0x1ac
138:    6e                      outs   dx,BYTE PTR ds:[esi]
139:    6e                      outs   dx,BYTE PTR ds:[esi]
13a:    20 6e 6e                and    BYTE PTR [esi+0x6e],ch
13d:    6e                      outs   dx,BYTE PTR ds:[esi]
13e:    6e                      outs   dx,BYTE PTR ds:[esi]
13f:    6e                      outs   dx,BYTE PTR ds:[esi]
140:    6e                      outs   dx,BYTE PTR ds:[esi]
141:    6e                      outs   dx,BYTE PTR ds:[esi]
142:    6e                      outs   dx,BYTE PTR ds:[esi]
143:    6e                      outs   dx,BYTE PTR ds:[esi]
144:    6e                      outs   dx,BYTE PTR ds:[esi]
145:    6e                      outs   dx,BYTE PTR ds:[esi]
146:    6e                      outs   dx,BYTE PTR ds:[esi]
147:    6e                      outs   dx,BYTE PTR ds:[esi]
148:    6e                      outs   dx,BYTE PTR ds:[esi]
149:    6e                      outs   dx,BYTE PTR ds:[esi]
14a:    6e                      outs   dx,BYTE PTR ds:[esi]
14b:    6e                      outs   dx,BYTE PTR ds:[esi]
14c:    6e                      outs   dx,BYTE PTR ds:[esi]
14d:    6e                      outs   dx,BYTE PTR ds:[esi]
14e:    6e                      outs   dx,BYTE PTR ds:[esi]
14f:    20 6e 6e                and    BYTE PTR [esi+0x6e],ch
152:    6e                      outs   dx,BYTE PTR ds:[esi]
153:    6e                      outs   dx,BYTE PTR ds:[esi]
154:    6e                      outs   dx,BYTE PTR ds:[esi]
155:    6e                      outs   dx,BYTE PTR ds:[esi]
156:    6e                      outs   dx,BYTE PTR ds:[esi]
157:    6e                      outs   dx,BYTE PTR ds:[esi]
158:    6e                      outs   dx,BYTE PTR ds:[esi]
159:    6e                      outs   dx,BYTE PTR ds:[esi]
15a:    20 6e 6e                and    BYTE PTR [esi+0x6e],ch
15d:    6e                      outs   dx,BYTE PTR ds:[esi]
15e:    6e                      outs   dx,BYTE PTR ds:[esi]
15f:    6e                      outs   dx,BYTE PTR ds:[esi]
160:    6e                      outs   dx,BYTE PTR ds:[esi]
161:    6e                      outs   dx,BYTE PTR ds:[esi]
162:    6e                      outs   dx,BYTE PTR ds:[esi]
163:    6e                      outs   dx,BYTE PTR ds:[esi]
164:    6e                      outs   dx,BYTE PTR ds:[esi]
165:    20 6e 6e                and    BYTE PTR [esi+0x6e],ch
168:    6e                      outs   dx,BYTE PTR ds:[esi]
169:    6e                      outs   dx,BYTE PTR ds:[esi]
16a:    6e                      outs   dx,BYTE PTR ds:[esi]
16b:    6e                      outs   dx,BYTE PTR ds:[esi]
16c:    6e                      outs   dx,BYTE PTR ds:[esi]
16d:    6e                      outs   dx,BYTE PTR ds:[esi]
16e:    6e                      outs   dx,BYTE PTR ds:[esi]
16f:    6e                      outs   dx,BYTE PTR ds:[esi]
170:    20 6e 6e                and    BYTE PTR [esi+0x6e],ch
173:    6e                      outs   dx,BYTE PTR ds:[esi]
174:    6e                      outs   dx,BYTE PTR ds:[esi]
175:    6e                      outs   dx,BYTE PTR ds:[esi]
176:    6e                      outs   dx,BYTE PTR ds:[esi]
177:    6e                      outs   dx,BYTE PTR ds:[esi]
178:    6e                      outs   dx,BYTE PTR ds:[esi]
179:    6e                      outs   dx,BYTE PTR ds:[esi]
17a:    6e                      outs   dx,BYTE PTR ds:[esi]
17b:    20 6e 6e                and    BYTE PTR [esi+0x6e],ch
17e:    6e                      outs   dx,BYTE PTR ds:[esi]
17f:    6e                      outs   dx,BYTE PTR ds:[esi]
180:    6e                      outs   dx,BYTE PTR ds:[esi]
181:    6e                      outs   dx,BYTE PTR ds:[esi]
182:    6e                      outs   dx,BYTE PTR ds:[esi]
183:    6e                      outs   dx,BYTE PTR ds:[esi]
184:    6e                      outs   dx,BYTE PTR ds:[esi]
185:    6e                      outs   dx,BYTE PTR ds:[esi]
186:    6e                      outs   dx,BYTE PTR ds:[esi]
187:    6e                      outs   dx,BYTE PTR ds:[esi]
188:    6e                      outs   dx,BYTE PTR ds:[esi]
189:    6e                      outs   dx,BYTE PTR ds:[esi]
18a:    6e                      outs   dx,BYTE PTR ds:[esi]
18b:    6e                      outs   dx,BYTE PTR ds:[esi]
18c:    6e                      outs   dx,BYTE PTR ds:[esi]
18d:    6e                      outs   dx,BYTE PTR ds:[esi]
18e:    6e                      outs   dx,BYTE PTR ds:[esi]
18f:    6e                      outs   dx,BYTE PTR ds:[esi]
190:    20 6e 6e                and    BYTE PTR [esi+0x6e],ch
193:    6e                      outs   dx,BYTE PTR ds:[esi]
194:    6e                      outs   dx,BYTE PTR ds:[esi]
195:    6e                      outs   dx,BYTE PTR ds:[esi]
196:    6e                      outs   dx,BYTE PTR ds:[esi]
197:    6e                      outs   dx,BYTE PTR ds:[esi]
198:    6e                      outs   dx,BYTE PTR ds:[esi]
199:    6e                      outs   dx,BYTE PTR ds:[esi]
19a:    6e                      outs   dx,BYTE PTR ds:[esi]
19b:    20 6e 6e                and    BYTE PTR [esi+0x6e],ch
19e:    6e                      outs   dx,BYTE PTR ds:[esi]
19f:    6e                      outs   dx,BYTE PTR ds:[esi]
1a0:    6e                      outs   dx,BYTE PTR ds:[esi]
1a1:    6e                      outs   dx,BYTE PTR ds:[esi]
1a2:    6e                      outs   dx,BYTE PTR ds:[esi]
1a3:    6e                      outs   dx,BYTE PTR ds:[esi]
1a4:    6e                      outs   dx,BYTE PTR ds:[esi]
1a5:    6e                      outs   dx,BYTE PTR ds:[esi]
1a6:    20 6e 6e                and    BYTE PTR [esi+0x6e],ch
1a9:    6e                      outs   dx,BYTE PTR ds:[esi]
1aa:    6e                      outs   dx,BYTE PTR ds:[esi]
1ab:    6e                      outs   dx,BYTE PTR ds:[esi]
1ac:    6e                      outs   dx,BYTE PTR ds:[esi]
1ad:    6e                      outs   dx,BYTE PTR ds:[esi]
1ae:    6e                      outs   dx,BYTE PTR ds:[esi]
1af:    6e                      outs   dx,BYTE PTR ds:[esi]
1b0:    61                      popa
1b1:    20 74 67 67             and    BYTE PTR [edi+eiz*2+0x67],dh
1b5:    61                      popa
1b6:    61                      popa
1b7:    61                      popa
1b8:    61                      popa
1b9:    67 67 61                addr16 addr16 popa
1bc:    20 61 67                and    BYTE PTR [ecx+0x67],ah
1bf:    67 67 61                addr16 addr16 popa
1c2:    61                      popa
1c3:    61                      popa
1c4:    61                      popa
1c5:    74 74                   je     0x23b
1c7:    74 63                   je     0x22c
1c9:    61                      popa
1ca:    61                      popa
1cb:    61                      popa
1cc:    61                      popa
1cd:    61                      popa
1ce:    74 74                   je     0x244
1d0:    67 20 67 67             and    BYTE PTR [bx+0x67],ah
1d4:    63 63 74                arpl   WORD PTR [ebx+0x74],sp
1d7:    67 61                   addr16 popa
1d9:    61                      popa
1da:    61                      popa
1db:    61                      popa
1dc:    20 74 63 63             and    BYTE PTR [ebx+eiz*2+0x63],dh
1e0:    61                      popa
1e1:    74 61                   je     0x244
1e3:    63 61 61                arpl   WORD PTR [ecx+0x61],sp
1e6:    74 20                   je     0x208
1e8:    61                      popa
1e9:    63 74 63 63             arpl   WORD PTR [ebx+eiz*2+0x63],si
1ed:    61                      popa
1ee:    67 74 61                addr16 je 0x252
1f1:    74 20                   je     0x213
1f3:    74 74                   je     0x269
1f5:    67 63 63 61             arpl   WORD PTR [bp+di+0x61],sp
1f9:    74 61                   je     0x25c
1fb:    61                      popa
1fc:    61                      popa
1fd:    20 67 61                and    BYTE PTR [edi+0x61],ah
200:    61                      popa
201:    61                      popa
202:    61                      popa
203:    61                      popa
204:    61                      popa
205:    67 61                   addr16 popa
207:    63 61 67                arpl   WORD PTR [ecx+0x67],sp
20a:    74 61                   je     0x26d
20c:    63 74 61 61             arpl   WORD PTR [ecx+eiz*2+0x61],si
210:    61                      popa
211:    74 20                   je     0x233
213:    67 67 61                addr16 addr16 popa
216:    67 61                   addr16 popa
218:    61                      popa
219:    61                      popa
21a:    61                      popa
21b:    74 74                   je     0x291
21d:    20 61 67                and    BYTE PTR [ecx+0x67],ah
220:    74 61                   je     0x283
222:    67 61                   addr16 popa
224:    74 74                   je     0x29a
226:    74 63                   je     0x28b
228:    20 61 67                and    BYTE PTR [ecx+0x67],ah
22b:    61                      popa
22c:    67 61                   addr16 popa
22e:    61                      popa
22f:    63 74 74 61             arpl   WORD PTR [esp+esi*2+0x61],si
233:    20 61 74                and    BYTE PTR [ecx+0x74],ah
236:    61                      popa
237:    61                      popa
238:    67 61                   addr16 popa
23a:    67 61                   addr16 popa
23c:    61                      popa
23d:    63 20                   arpl   WORD PTR [eax],sp
23f:    74 63                   je     0x2a4
241:    61                      popa
242:    61                      popa
243:    67 61                   addr16 popa
245:    63 74 74 63             arpl   WORD PTR [esp+esi*2+0x63],si
249:    74 67                   je     0x2b2
24b:    67 67 61                addr16 addr16 popa
24e:    61                      popa
24f:    67 74 74                addr16 je 0x2c6
252:    63 20                   arpl   WORD PTR [eax],sp
254:    61                      popa
255:    61                      popa
256:    74 74                   je     0x2cc
258:    61                      popa
259:    67 67 61                addr16 addr16 popa
25c:    61                      popa
25d:    74 20                   je     0x27f
25f:    61                      popa
260:    63 63 61                arpl   WORD PTR [ebx+0x61],sp
263:    63 61 74                arpl   WORD PTR [ecx+0x74],sp
266:    63 63 63                arpl   WORD PTR [ebx+0x63],sp
269:    20 67 63                and    BYTE PTR [edi+0x63],ah
26c:    61                      popa
26d:    67 67 67 74 74          addr16 addr16 addr16 je 0x2e6
272:    61                      popa
273:    61                      popa
274:    20 61 61                and    BYTE PTR [ecx+0x61],ah
277:    61                      popa
278:    61                      popa
279:    67 61                   addr16 popa
27b:    61                      popa
27c:    61                      popa
27d:    61                      popa
27e:    61                      popa
27f:    20 61 74                and    BYTE PTR [ecx+0x74],ah
282:    63 61 67                arpl   WORD PTR [ecx+0x67],sp
285:    74 61                   je     0x2e8
287:    61                      popa
288:    63 61 67                arpl   WORD PTR [ecx+0x67],sp
28b:    74 61                   je     0x2ee
28d:    63 74 67 67             arpl   WORD PTR [edi+eiz*2+0x67],si
291:    61                      popa
292:    74 67                   je     0x2fb
294:    20 74 67 67             and    BYTE PTR [edi+eiz*2+0x67],dh
298:    67 74 67                addr16 je 0x302
29b:    61                      popa
29c:    74 67                   je     0x305
29e:    63 20                   arpl   WORD PTR [eax],sp
2a0:    61                      popa
2a1:    74 61                   je     0x304
2a3:    74 74                   je     0x319
2a5:    74 74                   je     0x31b
2a7:    74 63                   je     0x30c
2a9:    61                      popa
2aa:    20 67 74                and    BYTE PTR [edi+0x74],ah
2ad:    74 63                   je     0x312
2af:    63 63 74                arpl   WORD PTR [ebx+0x74],sp
2b2:    74 61                   je     0x315
2b4:    67 20 61 74             and    BYTE PTR [bx+di+0x74],ah
2b8:    67 61                   addr16 popa
2ba:    61                      popa
2bb:    67 61                   addr16 popa
2bd:    63 74 74 20             arpl   WORD PTR [esp+esi*2+0x20],si
2c1:    63 61 67                arpl   WORD PTR [ecx+0x67],sp
2c4:    67 61                   addr16 popa
2c6:    61                      popa
2c7:    67 74 61                addr16 je 0x32b
2ca:    74 61                   je     0x32d
2cc:    63 74 67 63             arpl   WORD PTR [edi+eiz*2+0x63],si
2d0:    61                      popa
2d1:    74 74                   je     0x347
2d3:    74 61                   je     0x336
2d5:    20 63 63                and    BYTE PTR [ebx+0x63],ah
2d8:    61                      popa
2d9:    74 61                   je     0x33c
2db:    63 63 74                arpl   WORD PTR [ebx+0x74],sp
2de:    61                      popa
2df:    67 20 74 61             and    BYTE PTR [si+0x61],dh
2e3:    74 61                   je     0x346
2e5:    61                      popa
2e6:    61                      popa
2e7:    63 61 61                arpl   WORD PTR [ecx+0x61],sp
2ea:    74 20                   je     0x30c
2ec:    67 61                   addr16 popa
2ee:    67 61                   addr16 popa
2f0:    63 61 63                arpl   WORD PTR [ecx+0x63],sp
2f3:    63 61 67                arpl   WORD PTR [ecx+0x67],sp
2f6:    20 67 67                and    BYTE PTR [edi+0x67],ah
2f9:    61                      popa
2fa:    74 74                   je     0x370
2fc:    61                      popa
2fd:    67 61                   addr16 popa
2ff:    74 61                   je     0x362
301:    20 74 63 61             and    BYTE PTR [ebx+eiz*2+0x61],dh
305:    67 74 61                addr16 je 0x369
308:    63 61 61                arpl   WORD PTR [ecx+0x61],sp
30b:    74 67                   je     0x374
30d:    74 67                   je     0x376
30f:    63 74 74 63             arpl   WORD PTR [esp+esi*2+0x63],si
313:    63 61 63                arpl   WORD PTR [ecx+0x63],sp
316:    20 61 67                and    BYTE PTR [ecx+0x67],ah
319:    67 67 61                addr16 addr16 popa
31c:    74 67                   je     0x385
31e:    67 61                   addr16 popa
320:    61                      popa
321:    20 61 67                and    BYTE PTR [ecx+0x67],ah
324:    67 61                   addr16 popa
326:    74 63                   je     0x38b
328:    61                      popa
329:    63 63 61                arpl   WORD PTR [ebx+0x61],sp
32c:    20 67 63                and    BYTE PTR [edi+0x63],ah
32f:    61                      popa
330:    61                      popa
331:    74 61                   je     0x394
333:    74 74                   je     0x3a9
335:    63 63 20                arpl   WORD PTR [ebx+0x20],sp
338:    61                      popa
339:    61                      popa
33a:    61                      popa
33b:    67 74 61                addr16 je 0x39f
33e:    67 63 61 74             arpl   WORD PTR [bx+di+0x74],sp
342:    20 67 61                and    BYTE PTR [edi+0x61],ah
345:    63 61 61                arpl   WORD PTR [ecx+0x61],sp
348:    61                      popa
349:    61                      popa
34a:    61                      popa
34b:    74 63                   je     0x3b0
34d:    74 74                   je     0x3c3
34f:    61                      popa
350:    67 61                   addr16 popa
352:    67 63 63 74             arpl   WORD PTR [bp+di+0x74],sp
356:    74 20                   je     0x378
358:    74 74                   je     0x3ce
35a:    61                      popa
35b:    67 61                   addr16 popa
35d:    61                      popa
35e:    61                      popa
35f:    61                      popa
360:    63 61 20                arpl   WORD PTR [ecx+0x20],sp
363:    61                      popa
364:    61                      popa
365:    61                      popa
366:    74 63                   je     0x3cb
368:    63 61 67                arpl   WORD PTR [ecx+0x67],sp
36b:    61                      popa
36c:    63 20                   arpl   WORD PTR [eax],sp
36e:    61                      popa
36f:    74 61                   je     0x3d2
371:    67 74 74                addr16 je 0x3e8
374:    61                      popa
375:    74 63                   je     0x3da
377:    74 20                   je     0x399
379:    61                      popa
37a:    74 63                   je     0x3df
37c:    61                      popa
37d:    61                      popa
37e:    74 61                   je     0x3e1
380:    63 61 74                arpl   WORD PTR [ecx+0x74],sp
383:    20 67 67                and    BYTE PTR [edi+0x67],ah
386:    61                      popa
387:    74 67                   je     0x3f0
389:    61                      popa
38a:    74 74                   je     0x400
38c:    74 67                   je     0x3f5
38e:    74 61                   je     0x3f1
390:    74 67                   je     0x3f9
392:    74 61                   je     0x3f5
394:    67 67 61                addr16 addr16 popa
397:    74 20                   je     0x3b9
399:    63 74 67 61             arpl   WORD PTR [edi+eiz*2+0x61],si
39d:    63 74 74 61             arpl   WORD PTR [esp+esi*2+0x61],si
3a1:    67 61                   addr16 popa
3a3:    20 61 61                and    BYTE PTR [ecx+0x61],ah
3a6:    74 61                   je     0x409
3a8:    67 67 67 63 61 67       addr16 addr16 arpl WORD PTR [bx+di+0x67],sp
3ae:    20 63 61                and    BYTE PTR [ebx+0x61],ah
3b1:    74 61                   je     0x414
3b3:    67 61                   addr16 popa
3b5:    61                      popa
3b6:    63 61 61                arpl   WORD PTR [ecx+0x61],sp
3b9:    20 61 61                and    BYTE PTR [ecx+0x61],ah
3bc:    61                      popa
3bd:    74 61                   je     0x420
3bf:    67 61                   addr16 popa
3c1:    67 67 61                addr16 addr16 popa
3c4:    20 67 63                and    BYTE PTR [edi+0x63],ah
3c7:    74 67                   je     0x430
3c9:    61                      popa
3ca:    67 61                   addr16 popa
3cc:    63 61 61                arpl   WORD PTR [ecx+0x61],sp
3cf:    63 61 74                arpl   WORD PTR [ecx+0x74],sp
3d2:    63 74 67 74             arpl   WORD PTR [edi+eiz*2+0x74],si
3d6:    74 67                   je     0x43f
3d8:    61                      popa
3d9:    20 67 67                and    BYTE PTR [edi+0x67],ah
3dc:    74 67                   je     0x445
3de:    67 67 67 61             addr16 addr16 addr16 popa
3e2:    63 74 20 74             arpl   WORD PTR [eax+eiz*1+0x74],si
3e6:    61                      popa
3e7:    63 63 61                arpl   WORD PTR [ebx+0x61],sp
3ea:    63 61 63                arpl   WORD PTR [ecx+0x63],sp
3ed:    63 61 20                arpl   WORD PTR [ecx+0x20],sp
3f0:    67 61                   addr16 popa
3f2:    63 61 61                arpl   WORD PTR [ecx+0x61],sp
3f5:    61                      popa
3f6:    61                      popa
3f7:    61                      popa
3f8:    61                      popa
3f9:    63 20                   arpl   WORD PTR [eax],sp
3fb:    61                      popa
3fc:    74 63                   je     0x461
3fe:    61                      popa
3ff:    67 61                   addr16 popa
401:    61                      popa
402:    61                      popa
403:    67 61                   addr16 popa
405:    20 61 63                and    BYTE PTR [ecx+0x63],ah
408:    63 74 63 63             arpl   WORD PTR [ebx+eiz*2+0x63],si
40c:    61                      popa
40d:    74 74                   je     0x483
40f:    63 63 74                arpl   WORD PTR [ebx+0x74],sp
412:    74 74                   je     0x488
414:    67 67 61                addr16 addr16 popa
417:    74 67                   je     0x480
419:    67 20 67 74             and    BYTE PTR [bx+0x74],ah
41d:    74 61                   je     0x480
41f:    74 67                   je     0x488
421:    61                      popa
422:    61                      popa
423:    63 74 20 63             arpl   WORD PTR [eax+eiz*1+0x63],si
427:    63 61 74                arpl   WORD PTR [ecx+0x74],sp
42a:    63 63 74                arpl   WORD PTR [ebx+0x74],sp
42d:    67 61                   addr16 popa
42f:    74 20                   je     0x451
431:    61                      popa
432:    61                      popa
433:    61                      popa
434:    74 67                   je     0x49d
436:    67 61                   addr16 popa
438:    63                      .byte 0x63
439:    61                      popa