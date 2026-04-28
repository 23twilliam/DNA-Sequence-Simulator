import pytest
from src.simulators.illumina.simulator import IlluminaSimulator
from pathlib import Path

TEST_DIR = Path(__file__).resolve().parent
FILES_DIR = TEST_DIR / "test_files"


@pytest.fixture
def simulator():
    profile = FILES_DIR / "test_profile.txt"
    return IlluminaSimulator(read_length=100, coverage=30.0, profile_filename=str(profile))

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
