import pytest
from src.simulators.pacbio.simulator import _simulate_batch

def make_test_profile(read_length: int):
    """
    Create a minimal test profile.
    Every position and every base returns the same quality score.
    """
    profile = {}
    for pos in range(read_length):
        profile[pos] = {}
        for base in ['A', 'T', 'G', 'C']:
            profile[pos][base] = {
                "scores": [30],
                "cumulative": [1000],
                "total": 1000
            }
    return profile

def reverse_complement(sequence: str) -> str:
    """
    Return the reverse complement of sequence. (Borrowed from BaseSimulator)
    """
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(complement[base] for base in reversed(sequence))

@pytest.fixture
def batch_args():
    sequence = "ACGT" * 50  # 200bp
    read_length = 100
    return dict(
        batch_indices=list(range(5)),
        sequence_header="test_header",
        sequence=sequence,
        mean_read_length=read_length,
        coverage=10,
        cv=0.2,
        platform_name="test",
        profile=make_test_profile(read_length),
        homopolymer_pos={},
        reverse_complement_fn=reverse_complement,
    )

def test_batch_produces_reads(batch_args):
    reads = _simulate_batch(**batch_args)
    assert len(reads) > 0

def test_batch_read_ids_are_sequential(batch_args):
    reads = _simulate_batch(**batch_args)
    ids = list(r["id"] for r in reads)
    assert ids == [f"test_header_read_{i + 1}" for i in range(5)]

def test_batch_start_and_end_within_sequence_bounds(batch_args):
    reads = _simulate_batch(**batch_args)
    seq_len = len(batch_args["sequence"])

    for read in reads:
        assert 1 <= read["true_start"] <= seq_len
        assert read["true_end"] <= seq_len
        assert read["true_start"] < read["true_end"]

def test_quality_and_sequence_are_same_length(batch_args):
    # more important for pacbio than it was for illumina
    reads = _simulate_batch(**batch_args)
    for read in reads:
        assert len(read["sequence"]) == len(read["qualities"])

def test_read_structure_complete(batch_args):
    reads = _simulate_batch(**batch_args)
    expected_keys = {'id', 'description', 'sequence', 'qualities', 'true_start', 'true_end', 'errors'}
    for read in reads:
        assert expected_keys.issubset(read.keys())

def test_read_lengths_vary(batch_args):
    batch_args["batch_indices"] = list(range(50)) # do more reads to guarantee more variance
    reads = _simulate_batch(**batch_args)
    lengths = [r["true_end"] - r["true_start"] + 1 for r in reads]
    assert len(set(lengths)) > 1