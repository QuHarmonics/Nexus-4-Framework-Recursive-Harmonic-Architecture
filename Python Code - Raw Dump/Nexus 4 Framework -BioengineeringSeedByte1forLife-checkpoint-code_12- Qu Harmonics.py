from capstone import *

# Initialize Capstone for x86 32-bit mode
md = Cs(CS_ARCH_X86, CS_MODE_32)

def format_nasm():
    print("Enter raw hex bytes (e.g., '9f803d3dde5f5f6a7b'):")
    raw_hex = input("> ").strip()
    try:
        code = bytes.fromhex(raw_hex)
    except ValueError:
        print("Invalid input. Please enter valid hexadecimal bytes.")
        return

    print("\nNASM-Compatible Disassembly:")
    print("section .text")
    print("    global _start")
    print("\n_start:")
    for i in md.disasm(code, 0x1000):  # Disassemble starting at 0x1000
        mnemonic = i.mnemonic
        op_str = i.op_str
        print(f"    {mnemonic} {op_str}")

if __name__ == "__main__":
    format_nasm()
