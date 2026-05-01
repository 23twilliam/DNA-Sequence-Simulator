import pytest
from src.simulators.base_simulator import BaseSimulator
from pathlib import Path
import math

TEST_DIR = Path(__file__).resolve().parent
FILES_DIR = TEST_DIR / "test_files"
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

@pytest.fixture
def simulator():
    profile = FILES_DIR / "test_profile.txt"
    return Simulator(read_length=100, coverage=30, profile_path=profile, num_workers=1)

def test_homopolymer_detects_run(simulator):
    # run of 7 starting at position 0
    result = simulator.find_homopolymer_runs("AAAAAAATCG")
    assert 0 in result
    assert result[0] == 7

def test_homopolymer_min_length(simulator):
    # Default min_length = 3, so runs of 2 should not appear
    result = simulator.find_homopolymer_runs("AATCG")
    assert 0 not in result

def test_homopolymer_all_positions_in_run_marked(simulator):
    # Run of 5 A's all positions should be in dict
    result = simulator.find_homopolymer_runs("AAAAATCG")
    start = 5
    for i in range(5):
        assert result[i] == start - i

def test_homopolymer_multiple_runs(simulator):
    result = simulator.find_homopolymer_runs("AAAATTTTT")
    assert 0 in result
    assert result[0] == 4
    assert 4 in result
    assert result[4] == 5

def test_homopolymer_no_runs(simulator):
    result = simulator.find_homopolymer_runs("ATCGATCG")
    assert len(result) == 0

def test_homopolymer_override_min_length(simulator):
    # With min length = 5 a run of 4 should not appear in results
    result = simulator.find_homopolymer_runs("AAAATCG", min_length=5)
    assert 0 not in result

def test_homopolymer_returns_dict(simulator):
    assert isinstance(simulator.find_homopolymer_runs("AATCG"), dict)
