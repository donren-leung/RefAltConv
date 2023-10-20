import pytest
import pysam
from tests.test_context import fake_pysam_fastafile

@pytest.fixture
def ref_genome():
    return pysam.FastaFile('working_files/hg38.fa')

@pytest.mark.largefile
def test_mock_with_real_ref(ref_genome: pysam.FastaFile):
    for region, sequence in fake_pysam_fastafile.mocks.items():
        true_sequence = ref_genome.fetch(region=region)
        assert(fake_pysam_fastafile.fetch(region) == true_sequence)
        assert(sequence == true_sequence)
