class Emulator:
    def __init__(self):
        # Simplified registers (32-bit)
        self.reg = {
            'eax': 0,
            'ebx': 0,
            'ecx': 0,
            'edx': 0,
            'esi': 0,
            'edi': 0,
            'ebp': 0,
            'esp': 0x1000  # let's assume stack starts at address 0x1000
        }
        # We'll use a simple list to represent the stack.
        self.stack = []
        # Simulated memory as a dictionary: address -> value (byte or dword)
        self.memory = {}
        # Flags (for CMP, etc.)
        self.flags = {
            'ZF': 0,  # Zero flag
            'SF': 0   # Sign flag
        }
        # For simplicity, we'll simulate segment registers as always valid.
        self.es = self.memory  # For STOS, we treat ES as memory

    def push(self, value):
        # push a 32-bit value onto our stack
        self.stack.append(value)
        self.reg['esp'] -= 4

    def pop(self):
        if not self.stack:
            raise Exception("Stack underflow")
        self.reg['esp'] += 4
        return self.stack.pop()

    def execute(self):
        # Here we “simulate” our instruction stream.
        # Instruction at 0: "pop ebp"
        self.reg['ebp'] = self.pop()
        print("POP EBP -> EBP =", hex(self.reg['ebp']))

        # Instruction at 1: "cmp bh, dh"
        # We need to extract the 8-bit registers from EBX and EDX.
        # Let’s assume reg['ebx'] and reg['edx'] are already set.
        bh = (self.reg['ebx'] >> 8) & 0xFF
        dh = (self.reg['edx'] >> 8) & 0xFF
        result = bh - dh
        self.flags['ZF'] = int(result == 0)
        self.flags['SF'] = int((result & 0x80) != 0)
        print("CMP BH, DH -> BH =", hex(bh), "DH =", hex(dh), "ZF =", self.flags['ZF'])

        # Instruction at 3: "stos DWORD PTR es:[edi], eax"
        # Write EAX to memory at address in EDI (simulate ES:EDI = EAX)
        self.memory[self.reg['edi']] = self.reg['eax']
        print("STOS: Stored EAX =", hex(self.reg['eax']), "at address", hex(self.reg['edi']))

        # Instruction at 4: Byte 0x8f is not a valid instruction; skip it.
        print("Skipping invalid byte at address 4.")

        # Instruction at 5: "sbb BYTE PTR [ecx+0x56a45558], dl"
        addr = self.reg['ecx'] + 0x56A45558
        # For simulation, read a byte from memory (default 0 if not set)
        mem_val = self.memory.get(addr, 0)
        # For SBB, assume the carry flag (CF) is 0 for simplicity.
        cf = 0
        # Get DL (lower 8 bits of EDX)
        dl = self.reg['edx'] & 0xFF
        new_val = (mem_val - dl - cf) & 0xFF
        self.memory[addr] = new_val
        print(f"SBB: Memory[{hex(addr)}] changed from {hex(mem_val)} to {hex(new_val)} using DL = {hex(dl)}")

        # Instruction at B: "popa"
        # In 32-bit, POPA restores registers in the order:
        # EDI, ESI, EBP, ESP (skipped), EBX, EDX, ECX, EAX.
        # We simulate by popping 8 values from our stack into dummy registers.
        # For our simulation, assume the registers get restored to known values.
        if len(self.stack) >= 7:
            self.reg['edi'] = self.pop()
            self.reg['esi'] = self.pop()
            self.reg['ebp'] = self.pop()
            # Skip restoration of ESP as it's our stack pointer.
            self.reg['ebx'] = self.pop()
            self.reg['edx'] = self.pop()
            self.reg['ecx'] = self.pop()
            self.reg['eax'] = self.pop()
            print("POPA: Registers restored from stack.")
        else:
            print("POPA: Not enough values on stack; skipping.")

        # Instruction at C: "push 0x2a"
        self.push(0x2A)
        print("PUSH 0x2A")

        # Instruction at E: "cwde"
        # Convert word to dword: sign-extend AX to EAX.
        ax = self.reg['eax'] & 0xFFFF
        if ax & 0x8000:
            self.reg['eax'] = ax | 0xFFFF0000
        else:
            self.reg['eax'] = ax
        print("CWDE: EAX =", hex(self.reg['eax']))

        # Instruction at F: "add dl, dh"
        dl = self.reg['edx'] & 0xFF
        dh = (self.reg['edx'] >> 8) & 0xFF
        dl = (dl + dh) & 0xFF
        self.reg['edx'] = (self.reg['edx'] & 0xFFFFFF00) | dl
        print("ADD DL, DH -> DL =", hex(dl))

        # Instruction at 11: "dec ecx"
        self.reg['ecx'] = (self.reg['ecx'] - 1) & 0xFFFFFFFF
        print("DEC ECX -> ECX =", hex(self.reg['ecx']))

        # Instruction at 12: "pusha"
        # Push all registers (we simulate by saving them in order)
        saved_regs = (self.reg['eax'], self.reg['ecx'], self.reg['edx'],
                      self.reg['ebx'], self.reg['esp'], self.reg['ebp'],
                      self.reg['esi'], self.reg['edi'])
        for r in saved_regs:
            self.push(r)
        print("PUSHA: All registers pushed onto stack.")

        # Instruction at 13: "jae 0x8a"
        # For simulation, if ZF is 1 (or unsigned comparison shows above or equal), jump.
        # Our current flags from the earlier cmp: if ZF == 1, then condition is met.
        if self.flags['ZF'] == 1:
            print("JAE condition met; would jump to address 0x8a.")
        else:
            print("JAE condition not met; continue sequentially.")

        # Instruction at 15: "icebp"
        # ICEBP is typically used as a debugger breakpoint.
        print("ICEBP encountered - simulation breakpoint (no operation).")

        # Instruction at 16: "cdq"
        # CDQ sign-extends EAX into EDX.
        if self.reg['eax'] & 0x80000000:
            self.reg['edx'] = 0xFFFFFFFF
        else:
            self.reg['edx'] = 0
        print("CDQ: EDX =", hex(self.reg['edx']))

        # Instruction at 17: "cmps DWORD PTR ds:[esi], DWORD PTR es:[edi]"
        # For simulation, compare a DWORD at memory[esi] with memory[edi]
        val_esi = self.memory.get(self.reg['esi'], 0)
        val_edi = self.memory.get(self.reg['edi'], 0)
        result = val_esi - val_edi
        self.flags['ZF'] = int(result == 0)
        print("CMPS: Compared memory at ESI and EDI, result =", result)

        # Instruction at 18: "pop esp"
        # Pop a value from the stack and place it in ESP.
        self.reg['esp'] = self.pop()
        print("POP ESP -> ESP =", hex(self.reg['esp']))

        # Instruction at 19: "ret"
        # For simulation, ret pops a return address from the stack and “jumps” to it.
        ret_addr = self.pop()
        print("RET -> Would jump to address", hex(ret_addr))

        # Instruction at 1A: "sub eax, 0x2835bd7a"
        self.reg['eax'] = (self.reg['eax'] - 0x2835bd7a) & 0xFFFFFFFF
        print("SUB EAX, 0x2835bd7a -> EAX =", hex(self.reg['eax']))

        # Instruction at 1F: raw byte 0xbf; we treat it as data.
        print("Byte 0xbf encountered (treated as data).")

# Example usage:
emu = Emulator()
# Initialize registers for simulation; these values are arbitrary.
emu.reg['eax'] = 0x12345678
emu.reg['ebx'] = 0x87654321
emu.reg['ecx'] = 0x100
emu.reg['edx'] = 0x0A0B0C0D
emu.reg['esi'] = 0x200
emu.reg['edi'] = 0x300
emu.reg['ebp'] = 0x400

# Prepopulate stack with some values for POP instructions.
emu.stack = [0xDEADBEEF, 0xCAFEBABE, 0xFEEDFACE, 0xB16B00B5, 0x0BADF00D, 0x0D15EA5E, 0x12345678, 0x9ABCDEF0]

emu.execute()
