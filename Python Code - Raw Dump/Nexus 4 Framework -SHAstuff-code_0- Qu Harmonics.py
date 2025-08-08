import hashlib
import sys
from typing import List, Tuple, Optional

class ResonanceDifferentiator:
    """
    A tool to compute SHA-256 resonance metrics for a list of inputs.
    It hashes each input, converts the hex digest to ASCII-hex, reverses nibbles,
    and computes the delta to an optional target hash, reporting harmonicity.
    """
    def __init__(self, target_hash: Optional[str] = None):
        self.target_hash = target_hash.lower() if target_hash else None
        self.target_ascii_hex: Optional[str] = None
        if self.target_hash:
            self.target_ascii_hex = self.ascii_hexify(self.target_hash)

    @staticmethod
    def sha256_hexdigest(text: str) -> str:
        return hashlib.sha256(text.encode('utf-8')).hexdigest()

    @staticmethod
    def ascii_hexify(hex_str: str) -> str:
        return ''.join(format(ord(c), '02x') for c in hex_str)

    @staticmethod
    def reverse_nibbles(ascii_hex: str) -> str:
        return ascii_hex[::-1]

    @staticmethod
    def count_trailing_zeros(hex_str: str) -> int:
        # Count number of '0' characters at end of hex string
        return len(hex_str) - len(hex_str.rstrip('0'))

    def compute_metrics(self, text: str) -> dict:
        hexdigest = self.sha256_hexdigest(text)
        ascii_hex = self.ascii_hexify(hexdigest)
        reversed_nibs = self.reverse_nibbles(ascii_hex)
        zeros = self.count_trailing_zeros(reversed_nibs)
        metrics = {
            'input': text,
            'hash': hexdigest,
            'reversed_ascii_hex': reversed_nibs,
            'trailing_zeros': zeros,
        }
        if self.target_ascii_hex:
            # compute absolute numerical delta
            try:
                delta = abs(int(reversed_nibs, 16) - int(self.target_ascii_hex, 16))
                metrics['delta'] = delta
            except ValueError:
                metrics['delta'] = None
        return metrics

    def analyze_list(self, inputs: List[str]) -> List[dict]:
        return [self.compute_metrics(item) for item in inputs]


def main():
    if len(sys.argv) < 2:
        print("Usage: python resonance_differentiator.py <input1> [input2 ...] [--target <hash>]")
        sys.exit(1)

    args = sys.argv[1:]
    target = None
    if '--target' in args:
        idx = args.index('--target')
        if idx + 1 < len(args):
            target = args[idx + 1]
            args = args[:idx]

    rd = ResonanceDifferentiator(target)
    results = rd.analyze_list(args)
    for r in results:
        print(r)

if __name__ == '__main__':
    main()
