...

8:  00 07                   add    BYTE PTR [edi],al
a:  a8 b8                   test   al,0xb8
c:  f3 c6                   repz (bad)
e:  cb                      retf
f:  a0 7b 1d 9f 4d          mov    al,ds:0x4d9f1d7b
14: 78 a8                   js     0xffffffbe
16: 2d fd 3b 77 66          sub    eax,0x66773bfd
1b: e1 99                   loope  0xffffffb6
1d: d2                      .byte 0xd2
1e: a4                      movs   BYTE PTR es:[edi],BYTE PTR ds:[esi]
1f: 31                      .byte 0x31