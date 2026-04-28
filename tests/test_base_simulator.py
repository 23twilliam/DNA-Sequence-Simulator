import pytest
from src.simulators.base_simulator import BaseSimulator
from pathlib import Path
import math

class Simulator(BaseSimulator):
    def _precompute(self, sequence): pass
    def _get_worker_args(self, batch, sequence_header, sequence): return (batch, sequence_header, sequence)
    def _get_worker_fn(self):
        def worker(batch, sequence_header, sequence):
            return [{
                "id": f"{sequence_header}_read_{i + 1}",
                "description": "test read",
                "sequence": sequence[:10],
                "qualities": [30] * 10,
                "true_start": 1,
                "true_end": 10,
                "strand": "F",
                "errors": []
            } for i in batch]
        return worker
    def get_platform_name(self): return "Test"
    def generate_quality_scores(self, read): return [30] * len(read)
    def inject_errors(self, seq, quals): return seq, quals, []

TEST_DIR = Path(__file__).resolve().parent
FILES_DIR = TEST_DIR / "test_files"

@pytest.fixture
def sim():
    profile = FILES_DIR / "test_profile.txt"
    return Simulator(read_length=100, coverage=30, profile_path=profile, num_workers=1)

def test_read_length_exceeds_sequence(sim):
    with pytest.raises(ValueError):
        sim.generate_reads(sequence_header="test_seq", sequence="ATGC")

def test_read_structure_complete(sim):
    sequence = "ATGC" * 25
    reads = sim.generate_reads(sequence_header="test_seq", sequence=sequence)
    expected_keys = {'id', 'description', 'sequence', 'qualities', 'true_start', 'true_end', 'strand', 'errors'}
    for read in reads:
        assert expected_keys.issubset(read.keys())

def test_simulate_from_file(sim):
    class MockIO:
        def read_sequences(_):
            sequence = "ATGC" * 25
            return {"test_seq": sequence}

    reads = sim.simulate_from_file(fasta_file="fake.fasta", io_functions=MockIO)
    print(reads)
    assert len(reads) > 0

########################################################################################################################

# calculate_num_reads tests
def test_num_reads_basic(sim):
    # (seq_len * coverage) / read_length = (3000 * 30) / 100 = 900
    assert sim.calculate_num_reads(3000) == 900

def test_num_reads_truncates_not_rounds(sim):
    # Should use ceil
    result = sim.calculate_num_reads(101)
    assert isinstance(result, int)
    assert result == math.ceil((101 * 30) / 100)

def test_num_reads_coverage_effect():
    profile = FILES_DIR / "test_profile.txt"
    sim_high = Simulator(read_length=100, coverage=60, profile_path=profile)
    sim_low = Simulator(read_length=100, coverage=30, profile_path=profile)
    assert sim_high.calculate_num_reads(1000) == 2 * sim_low.calculate_num_reads(1000)

########################################################################################################################

# reverse_complement tests
def test_reverse_complement_normal(sim):
    assert sim.reverse_complement("ATCG") == "CGAT"

def test_reverse_complement_all_same(sim):
    assert sim.reverse_complement("AAAA") == "TTTT"

def test_reverse_complement_palindrome(sim):
    # ATAT reverse complement is ATAT
    assert sim.reverse_complement("ATAT") == "ATAT"

def test_reverse_complement_n(sim):
    # N should still map to N
    assert sim.reverse_complement("ATNG") == "CNAT"

def test_reverse_complement_double(sim):
    # Applying twice should give return the original sequence
    seq = "ATCGATCG"
    assert sim.reverse_complement(sim.reverse_complement(seq)) == seq

########################################################################################################################

# calculate_gc_content tests
def test_gc_content_all_gc(sim):
    assert sim.calculate_gc_content("GCGCGC") == 1.0

def test_gc_content_no_gc(sim):
    assert sim.calculate_gc_content("ATATAT") == 0.0

def test_gc_content_half(sim):
    assert sim.calculate_gc_content("ATGC") == pytest.approx(0.5)

def test_quality_to_error_probability(sim):
    # Q30 should give 0.001
    assert sim.quality_to_error_probability(30) == pytest.approx(0.001)
    # Q20 should give 0.01
    assert sim.quality_to_error_probability(20) == pytest.approx(0.01)
    # Q10 should give 0.1
    assert sim.quality_to_error_probability(10) == pytest.approx(0.1)