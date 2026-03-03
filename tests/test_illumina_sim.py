import pytest
from src.simulators.illumina.simulator import IlluminaSimulator
from pathlib import Path

TEST_DIR = Path(__file__).resolve().parent
FILES_DIR = TEST_DIR / "test_files"

@pytest.fixture
def simulator():
    profile = FILES_DIR / "art_illumina_profile_hiseq.txt"
    return IlluminaSimulator(read_length=100, coverage=30.0, profile_filename=str(profile))

def test_something(simulator):
    assert simulator.read_length == 100