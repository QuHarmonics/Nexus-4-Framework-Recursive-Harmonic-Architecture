5d38f7ab8f18915855a456616a2a9800f249607375f199a75cc32d7abd3528bf

pop ebp
cmp bh, dh
stosd
nop
sbb byte [ecx+0x56a45558], dl
popad
push 0x2a
cwde
add dl, dh
dec ecx
pusha
jae L8a
nop
cdq
cmpsd
pop esp
ret
L8a:
sub eax, 0x2835bd7a
nop


5D38F7AB9018915955A456616A2A9800F2496073059099A75CC32D7ABD352890     from complied


