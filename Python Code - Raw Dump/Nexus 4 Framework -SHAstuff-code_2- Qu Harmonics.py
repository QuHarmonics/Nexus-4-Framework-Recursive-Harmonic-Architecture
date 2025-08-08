import hashlib
import logging

logging.basicConfig(level=logging.INFO)

class ResonanceDifferentiator:
    """
    A SHA-harmonic mapping engine that processes inputs, computes SHA-256 hashes,
    converts them to ASCII-hex, reverses nibbles, computes deltas, and logs harmonicity metrics.
    """

    def __init__(self, target_hash: str = None):
        """
        Initialize with an optional target hash.
        :param target_hash: SHA-256 hex digest string (64 hex chars)
        """
        self.target_hash = target_hash
        if target_hash:
            self._target_rev_int = self._reverse_and_int(target_hash)
        else:
            self._target_rev_int = None

    @staticmethod
    def _sha256_hexdigest(data: str) -> str:
        return hashlib.sha256(data.encode('utf-8')).hexdigest()

    @staticmethod
    def _ascii_hexify(hex_str: str) -> str:
        # Convert each hex char into its two-digit ASCII hex code
        return ''.join(format(ord(c), '02x') for c in hex_str)

    @staticmethod
    def _reverse_nibbles(ascii_hex: str) -> str:
        return ascii_hex[::-1]

    def _reverse_and_int(self, hex_digest: str) -> int:
        ascii_hex = self._ascii_hexify(hex_digest)
        rev = self._reverse_nibbles(ascii_hex)
        return int(rev, 16)

    @staticmethod
    def _count_trailing_zeros(num: int) -> int:
        # Count trailing zero bits in binary representation
        if num == 0:
            return 0
        count = 0
        while (num & 1) == 0:
            num >>= 1
            count += 1
        return count

    def process(self, inputs: list):
        """
        Process a list of input strings, computing and logging harmonicity metrics.
        """
        for data in inputs:
            digest = self._sha256_hexdigest(data)
            ascii_hex = self._ascii_hexify(digest)
            rev_int = self._reverse_and_int(digest)

            metrics = {'input': data, 'sha256': digest}

            if self._target_rev_int is not None:
                delta = abs(rev_int - self._target_rev_int)
                zeros = self._count_trailing_zeros(delta)
                is_power16 = (delta != 0 and delta & (delta - 1) == 0 and delta.bit_length() % 4 == 0)
                metrics.update({'delta': delta, 'trailing_zero_bits': zeros, 'power_of_16': is_power16})

            logging.info("Harmonicity: %s", metrics)

        return

# Example usage
if __name__ == '__main__':
    # Initialize with an optional target hash to compare against
    target = None  # or set to a 64-character SHA-256 hex digest
    rd = ResonanceDifferentiator(target)
    rd.process(["Hello", "hello", "World"])
