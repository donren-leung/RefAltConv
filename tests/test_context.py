import pytest
from varconv.variants import Variant, VariantContext
import pysam
import Bio.Seq

from unittest.mock import MagicMock

# Source: hg38 ucsc genome browser
# See test_context_large for verification of these mocks
class fake_pysam_fastafile(dict):
    mocks = {
        #AK5, context window 7
        'chr1:77297573-77297587' :  'GGAAAGGGTACTCAG',
        'chr1:77297572-77297587' : 'TGGAAAGGGTACTCAG',
        'chr1:77297573-77297588' :  'GGAAAGGGTACTCAGA',
        #AK5, context window 9
        'chr1:77297571-77297589' :  'GTGGAAAGGGTACTCAGAG',
        'chr1:77297570-77297589' : 'AGTGGAAAGGGTACTCAGAG',
        'chr1:77297571-77297590' :  'GTGGAAAGGGTACTCAGAGT',
        
        'chr1:1393386-1393408'   : 'CTCAACTCACCATGAGGTCTGGA',
        'chr1:42189293-42189307' : 'CAACACTCACCTTGA',

        # Reverse direction genes, context window 7
        # PHP4
        'chr1:5873231-5873245'   : 'TGTCCTCCGTTGCCC',
        # CCNL2
        'chr1:1393386-1393408'   : 'CTCAACTCACCATGAGGTCTGGA',
        # FOXJ3
        'chr1:42189293-42189307' : 'CAACACTCACCTTGA'
    }

    @staticmethod
    def fetch(region: str) -> str:
        try:
            return __class__.mocks[region]
        except KeyError:
            raise KeyError(f"Mock does not support {region=}")

@pytest.fixture
def ref_genome():
    mock_fasta_file = MagicMock(spec=pysam.FastaFile)
    mock_fasta_file.fetch.side_effect = fake_pysam_fastafile.fetch
    yield mock_fasta_file

def test_fake_pysam(ref_genome: pysam.FastaFile):
    def test(variant: Variant, sequence: str, context_length: int):
        variant_context = VariantContext(ref_genome, variant, context_length)
        assert(variant_context.ref_sequence_fwd() == sequence)

    CONTEXT_LENGTH = 7
    # Variant() will fetch the region chr1:77297573-77297587 using pysam.FastaFile.fetch()
    test(Variant("chr1", 77297580, "G", "T"), 'GGAAAGGGTACTCAG', CONTEXT_LENGTH)
    with pytest.raises(ValueError):
        test(Variant("chr1", 77297580, "A", "T"), 'GGAAAGGGTACTCAG', CONTEXT_LENGTH)

def test_VariantContext(ref_genome: pysam.FastaFile):
    def test(variant: Variant, context_length: int, expected: tuple):
        variant_context = VariantContext(ref_genome, variant, context_length)
        assert(variant_context.debug() == expected)

    CONTEXT_LENGTH = 7
    test(Variant("chr1", 77297580, "G", "T"), CONTEXT_LENGTH, (('GGAAAGG', 'G', 'TACTCAG', 'T'), [7, 1, 7, 1]))
    test(Variant("chr1", 77297580, "G", "GA"), CONTEXT_LENGTH, (('GGAAAGG', 'G', 'TACTCAG', 'GA'), [7, 1, 7, 2]))
    test(Variant("chr1", 77297579, "GG", "TA"), CONTEXT_LENGTH, (('TGGAAAG', 'GG', 'TACTCAG', 'TA'), [7, 2, 7, 2]))
    test(Variant("chr1", 77297580, "GT", "A"), CONTEXT_LENGTH, (('GGAAAGG', 'GT', 'ACTCAGA', 'A'), [7, 2, 7, 1]))

    test(Variant("chr1", 1393393, "CACCATGAG", "C"), CONTEXT_LENGTH, (('CTCAACT', 'CACCATGAG', 'GTCTGGA', 'C'), [7, 9, 7, 1]))
    test(Variant("chr1", 42189300, "C", "CA"), CONTEXT_LENGTH, (('CAACACT', 'C', 'ACCTTGA', 'CA'), [7, 1, 7, 2]))

    CONTEXT_LENGTH = 9
    test(Variant("chr1", 77297580, "G", "T"), CONTEXT_LENGTH, (('GTGGAAAGG', 'G', 'TACTCAGAG', 'T'),  [9, 1, 9, 1]))
    test(Variant("chr1", 77297580, "G", "GA"), CONTEXT_LENGTH, (('GTGGAAAGG', 'G', 'TACTCAGAG', 'GA'),  [9, 1, 9, 2]))
    test(Variant("chr1", 77297579, "GG", "TA"), CONTEXT_LENGTH, (('AGTGGAAAG', 'GG', 'TACTCAGAG', 'TA'), [9, 2, 9, 2]))
    test(Variant("chr1", 77297580, "GT", "A"), CONTEXT_LENGTH, (('GTGGAAAGG', 'GT', 'ACTCAGAGT', 'A'), [9, 2, 9, 1]))

def test_fwd_sequences(ref_genome: pysam.FastaFile):
    from dataclasses import dataclass
    @dataclass
    class Expected(dict):
        ref_seq_fwd: str
        alt_seq_fwd: str

    def test_forward(variant: Variant, context_length: int, expected: Expected):
        variant_context = VariantContext(ref_genome, variant, context_length)
        assert(variant_context.ref_sequence_fwd() == expected.ref_seq_fwd)
        assert(variant_context.alt_sequence_fwd() == expected.alt_seq_fwd)

    CONTEXT_LENGTH = 7
    test_forward(Variant("chr1", 77297580, "G", "T"), CONTEXT_LENGTH,
             Expected('GGAAAGG' + 'G' + 'TACTCAG', 'GGAAAGG' + 'T' + 'TACTCAG'))
    test_forward(Variant("chr1", 77297580, "G", "GA"), CONTEXT_LENGTH,
             Expected('GGAAAGG' + 'G' + 'TACTCAG','GGAAAGG' + 'GA' + 'TACTCAG'))
    test_forward(Variant("chr1", 77297579, "GG", "TA"), CONTEXT_LENGTH,
             Expected('TGGAAAG' + 'GG' + 'TACTCAG','TGGAAAG' + 'TA' + 'TACTCAG'))
    test_forward(Variant("chr1", 77297580, "GT", "A"), CONTEXT_LENGTH,
             Expected('GGAAAGG' + 'GT' + 'ACTCAGA','GGAAAGG' + 'A' + 'ACTCAGA'))

def test_reverse_sequences(ref_genome: pysam.FastaFile):
    from dataclasses import dataclass
    @dataclass
    class Expected(dict):
        ref_seq_rev: str
        alt_seq_rev: str

    def test_reverse(variant: Variant, context_length: int, expected: Expected):
        variant_context = VariantContext(ref_genome, variant, context_length)
        assert(variant_context.ref_sequence_rev() == Bio.Seq.reverse_complement(expected.ref_seq_rev))
        assert(variant_context.alt_sequence_rev() == Bio.Seq.reverse_complement(expected.alt_seq_rev))

    CONTEXT_LENGTH = 7
    # Reverse strand, supply Expected() as forward strand sequence, sourced from ucsc DNA browser
    test_reverse(Variant("chr1", 5873238, "C", "T"), CONTEXT_LENGTH,
             Expected('TGTCCTC' + 'C' + 'GTTGCCC', 'TGTCCTC' + 'T' + 'GTTGCCC'))
    test_reverse(Variant("chr1", 1393393, "CACCATGAG", "C"), CONTEXT_LENGTH,
             Expected('CTCAACT' + 'CACCATGAG' + 'GTCTGGA', 'CTCAACT' + 'C' + 'GTCTGGA'))
    test_reverse(Variant("chr1", 42189300, "C", "CA"), CONTEXT_LENGTH,
             Expected('CAACACT' + 'C' + 'ACCTTGA', 'CAACACT' + 'CA' + 'ACCTTGA'))
