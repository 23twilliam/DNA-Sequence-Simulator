import pytest
from src.simulators.illumina.simulator import IlluminaSimulator
from pathlib import Path

TEST_DIR = Path(__file__).resolve().parent
FILES_DIR = TEST_DIR / "test_files"


@pytest.fixture
def simulator():
	profile = FILES_DIR / "test_profile.txt"
	return IlluminaSimulator(read_length=100, coverage=30.0, profile_filename=str(profile))


def test_something(simulator):
	assert simulator.read_length == 100


@pytest.fixture
def profile():
	profile = FILES_DIR / "test_profile.txt"
	simulator = IlluminaSimulator(read_length=100, coverage=30.0, profile_filename=str(profile))
	result = simulator.load_quality_profile(profile)
	return result


def test_profile_exists(profile):
	assert len(profile[0]["A"]) > 0


def test_scores_parse_correctly(profile):
	assert profile[0]["A"]["scores"] == [10, 20, 30, 40]


def test_cumulative_count(profile):
	assert profile[0]["."]["cumulative"] == [4, 12, 24, 40]
	assert profile[0]["A"]["cumulative"] == [1, 3, 6, 10]


def test_total(profile):
	assert profile[0]["."]["total"] == 40
	assert profile[0]["A"]["total"] == 10


def test_position_uniqueness(profile):
	assert profile[0]["A"] is not profile[1]["A"]


@pytest.fixture
def qualities():
	profile = FILES_DIR / "art_illumina_profile_hiseq.txt"
	gene = "CTCAAAAGTCTAGAGCCACCGTCCAGGGAGCAGGTAGCTGCTGGGCTCCGGGGACACTTTGCGTTCGGGC"
	simulator = IlluminaSimulator(read_length=100, coverage=30.0, profile_filename=str(profile))
	result = simulator.generate_quality_scores(gene)
	return result

def test_quality_length(qualities):
	assert len(qualities) == 70