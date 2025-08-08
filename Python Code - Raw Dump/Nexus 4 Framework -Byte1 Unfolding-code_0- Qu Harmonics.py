
from textwrap import dedent
from capstone import Cs, CS_ARCH_X86, CS_MODE_32

RAW_LISTING = dedent("""
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
""").strip().splitlines()

def listing_to_bytes(lines):
    """
    Extract just the hex-byte tokens, return as real bytes.
    """
    byte_strs = []
    for ln in lines:
        # split on the two spaces after the address field
        try:
            _, hex_chunk, _ = ln.split(None, 2)
        except ValueError:
            # 1-byte rows like "1c: 61                      popa"
            _, hex_chunk = ln.split(None, 1)
        # drop repeating whitespace, keep hex pairs
        parts = hex_chunk.strip().split()
        byte_strs.extend(parts)
    return bytes.fromhex("".join(byte_strs))

code_bytes = listing_to_bytes(RAW_LISTING)

# sanity-check disassembly
md = Cs(CS_ARCH_X86, CS_MODE_32)
print("---- quick re-disassembly check ----")
for i, ins in enumerate(md.disasm(code_bytes, 0)):
    print(f"{ins.address:04x}:\t{ins.mnemonic}\t{ins.op_str}")

# optionally write to disk if you want to feed Ghidra/r2
open("snippet.bin", "wb").write(code_bytes)
print("\nWrote raw bytes to snippet.bin")
