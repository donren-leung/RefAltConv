import pysam
import sys
from dataclasses import dataclass
from enum import Enum, auto
from Bio.Seq import reverse_complement
from typing import Optional

class StrandDirection(Enum):
    FORWARD = auto()
    REVERSE = auto()
    BOTH = auto()
    UNKNOWN = auto()

# Variants extracted from vcf file
@dataclass(frozen=True)
class Variant():
    chrom: str
    position: int
    ref_allele: str
    alt_allele: str

    strand_direction: StrandDirection = StrandDirection.UNKNOWN

# Variants with ref/alt allele and context window before and after
class VariantContext():
    variant: Variant
    context_length: int
    before: str
    after: str

    def __init__(self, ref_genome: pysam.FastaFile, variant: Variant,
                 context_length: int, *, exit_on_mismatch: bool = False) -> None:
        self.variant = variant
        self.context_length = context_length

        ref_length = len(variant.ref_allele)
        start = variant.position - context_length
        end = variant.position + ref_length + context_length - 1
        
        region = f"{variant.chrom}:{start}-{end}"
        sequence = ref_genome.fetch(region=region).upper()

        self.before = sequence[:context_length]
        try:
            assert(variant.ref_allele == sequence[context_length:context_length + ref_length])
        except AssertionError:
            print(f"Variant reference allele ({variant.ref_allele}) did not match the provided reference genome "
                f"({sequence[context_length:context_length + ref_length]}) at {variant.chrom}:{variant.position}"
                "... skipping or exiting", file=sys.stderr)            
            
            raise ValueError

        self.after = sequence[context_length + ref_length:]

    # TODO: find a way to not have to remake strings for context window
    # TODO: assert context_length is long enough for window_length
    def ref_sequence_fwd(self, window_length: Optional[int]=None) -> str:
        if window_length is not None:
            window_length -= 1
            return f"{self.before[-window_length:]}{self.ref_allele}{self.after[:window_length]}"
        return f"{self.before}{self.variant.ref_allele}{self.after}"
    
    def alt_sequence_fwd(self, window_length: Optional[int]=None) -> str:
        if window_length is not None:
            window_length -= 1
            return f"{self.before[-window_length:]}{self.alt_allele}{self.after[:window_length]}"
        return f"{self.before}{self.variant.alt_allele}{self.after}"
    
    def ref_sequence_rev(self, window_length: Optional[int]=None) -> str:
        return reverse_complement(self.ref_sequence_fwd(window_length))

    def alt_sequence_rev(self, window_length: Optional[int]=None) -> str:
        return reverse_complement(self.alt_sequence_fwd(window_length))

    def debug(self) -> tuple:
        tup = (self.before, self.variant.ref_allele, self.after, self.variant.alt_allele)
        return tup, [len(elem) for elem in tup]
    
    def __repr__(self) -> str:
        return str(self.debug())
    

def main():
    pass

if __name__ == '__main__':
    main()
