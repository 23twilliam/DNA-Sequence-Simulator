import pytest
from src.simulators.illumina.simulator import IlluminaSimulator
from pathlib import Path

TEST_DIR = Path(__file__).resolve().parent
FILES_DIR = TEST_DIR / "test_files"


@pytest.fixture
def simulator():
	profile = FILES_DIR / "test_profile.txt"
	gene = "CTCAAAAGTCTAGAGCCACCGTCCAGGGAGCAGGTAGCTGCTGGGCTCCGGGGACACTTTGCGTTCGGGC"
	return IlluminaSimulator(read_length=100, coverage=30.0, profile_filename=str(profile), sequence=gene)


def test_something(simulator):
	assert simulator.read_length == 100


@pytest.fixture
def profile():
	profile = FILES_DIR / "test_profile.txt"
	gene = "CTCAAAAGTCTAGAGCCACCGTCCAGGGAGCAGGTAGCTGCTGGGCTCCGGGGACACTTTGCGTTCGGGC"
	simulator = IlluminaSimulator(read_length=100, coverage=30.0, profile_filename=str(profile), sequence=gene)
	result = simulator.load_quality_profile(str(profile))
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
	simulator = IlluminaSimulator(read_length=100, coverage=30.0, profile_filename=str(profile), sequence=gene)
	result = simulator.generate_quality_scores(gene)
	return result

def test_quality_length(qualities):
	assert len(qualities) == 70

@pytest.fixture
def ggc_simulator():
	profile = FILES_DIR / "art_illumina_profile_hiseq.txt"
	gene = "AAAAA" + "GGC" + "A" * 100
	return IlluminaSimulator(read_length=50, coverage=10.0, profile_filename=str(profile), sequence=gene)

@pytest.fixture
def ir_simulator():
	profile = FILES_DIR / "art_illumina_profile_hiseq.txt"
	forward = "TTACGGTAA"
	reverse_c = "TTACCGTAA"
	gene = "AAAAA" + forward + "A" * 20 + reverse_c + "A" * 100
	return IlluminaSimulator(read_length=50, coverage=10.0, profile_filename=str(profile), sequence=gene)

@pytest.fixture
def combined_simulator():
	profile = FILES_DIR / "art_illumina_profile_hiseq.txt"
	forward = "TTACGGTAA"
	reverse_c = "TTACCGTAA"
	gene = "AAAGGC" + forward + "A" * 20 + reverse_c + "A" * 100
	return IlluminaSimulator(read_length=50, coverage=10.0, profile_filename=str(profile), sequence=gene)


def test_ggc_trigger(ggc_simulator):
	# GGC at position 5, impact area should cover 5-44
	assert 5 in ggc_simulator.sse_trigger_map
	assert 44 in ggc_simulator.sse_trigger_map
	assert 4 not in ggc_simulator.sse_trigger_map  # before trigger


def test_ir_trigger(ir_simulator):
	# Inverted repeats start at position 5, impact are should cover 5-44
	assert 5 in ir_simulator.sse_trigger_map
	assert 44 in ir_simulator.sse_trigger_map
	assert 4 not in ir_simulator.sse_trigger_map  # before trigger


def test_combined_trigger(combined_simulator):
	# GGC at position 3, IR at position 6
	assert 3 in combined_simulator.sse_trigger_map
	assert 6 in combined_simulator.sse_trigger_map
	assert 2 not in combined_simulator.sse_trigger_map  # before trigger
