0:  33 c0                   xor    eax,eax
2:  f7 e2                   mul    edx
4:  c7 45 fc 00 00 00 00    mov    DWORD PTR [ebp-0x4],0x0
b:  ee                      out    dx,al
c:  eb fe                   jmp    0xc
e:  31 db                   xor    ebx,ebx
10: 8b c3                   mov    eax,ebx
12: 20 20                   and    BYTE PTR [eax],ah
14: 75 f6                   jne    0xc