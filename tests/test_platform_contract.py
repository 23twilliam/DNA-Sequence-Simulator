import pytest
import os
from src.simulators.illumina.simulator import IlluminaSimulator
from src.simulators.pacbio.simulator import PacBioSimulator
from pathlib import Path

TEST_DIR = Path(__file__).resolve().parent
FILES_DIR = TEST_DIR / "test_files"

@pytest.fixture
def illumina_simulator():
    profile = FILES_DIR / "ecoli_hiseqR1.txt"
    if not os.path.exists(profile):
        pytest.fail("Illumina profile file not found")
    return IlluminaSimulator(
        read_length=100, coverage=1, profile_filename=str(profile)
    )


@pytest.fixture
def pacbio_simulator():
    profile = FILES_DIR / "ecoli_hiseqR1.txt"  # whatever your PacBio profile is named
    if not os.path.exists(profile):
        pytest.fail("PacBio profile file not found")
    return PacBioSimulator(
        mean_read_length=100, coverage=1, profile_filename=str(profile)
    )

@pytest.fixture
def short_sequence():
    # Small realistic sequence with some GGC motifs and homopolymers
    return "ATCGGGCAAAAAATCGATCGATCG" * 20  # 480bp (enough for 100bp reads)

def test_pipeline_produces_reads(illumina_simulator, short_sequence):
    reads = illumina_simulator.generate_reads(sequence_header="test_seq", sequence=short_sequence)
    assert len(reads) > 0

def test_read_structure_complete(illumina_simulator, short_sequence):
    reads = illumina_simulator.generate_reads(sequence_header="test_seq", sequence=short_sequence)
    expected_keys = {'id', 'description', 'sequence', 'qualities', 'strand', 'true_start', 'true_end', 'errors'}
    for read in reads:
        assert expected_keys.issubset(read.keys())

def test_read_length_within_expected_range(illumina_simulator, short_sequence):
    reads = illumina_simulator.generate_reads(sequence_header="test_seq", sequence=short_sequence)
    for read in reads:
        # reads can be trimmed by prephasing so we allow the reads to be in the range 97 - 100
        assert 97 <= len(read["sequence"]) <= 100

def test_quality_and_sequence_are_same_length(illumina_simulator, short_sequence):
    reads = illumina_simulator.generate_reads(sequence_header="test_seq", sequence=short_sequence)
    for read in reads:
        assert len(read["sequence"]) == len(read["qualities"])

def test_coordinates_within_sequence_bounds(illumina_simulator, short_sequence):
    reads = illumina_simulator.generate_reads(sequence_header="test_seq", sequence=short_sequence)
    seq_len = len(short_sequence)
    for read in reads:
        assert 1 <= read["true_start"] <= seq_len
        assert read["true_end"] <= seq_len
        assert read["true_start"] < read["true_end"]

def test_coverage_produces_expected_read_count(illumina_simulator):
    seq = "ATCG" * 1000 # 4000bp
    illumina_simulator.coverage = 10
    reads = illumina_simulator.generate_reads(sequence_header="test_seq", sequence=seq)
    expected = int((4000 * 10) / 100)
    assert len(reads) == expected

def test_read_length_exceeds_sequence_raises_error(illumina_simulator):
    with pytest.raises(ValueError):
        illumina_simulator.generate_reads(sequence_header="test_seq", sequence="ATCG") # Too short (< 100)

def test_all_bases_valid(illumina_simulator, short_sequence):
    reads = illumina_simulator.generate_reads(sequence_header="test_seq", sequence=short_sequence)
    valid_bases = set('ATCGN')
    for read in reads:
        assert set(read["sequence"]).issubset(valid_bases)

def test_quality_scores_in_phred_range(illumina_simulator, short_sequence):
    reads = illumina_simulator.generate_reads(sequence_header="test_seq", sequence=short_sequence)
    for read in reads:
        for q in read["qualities"]:
            assert 2 <= q <= 40 # Our floor is 2 and illumina max is 40


# ====================================================================== #
# PacBio contract tests
# ====================================================================== #

def test_pacbio_pipeline_produces_reads(pacbio_simulator, short_sequence):
    reads = pacbio_simulator.generate_reads("test_seq", short_sequence)
    assert len(reads) > 0


def test_pacbio_read_structure_complete(pacbio_simulator, short_sequence):
    reads = pacbio_simulator.generate_reads("test_seq", short_sequence)
    expected_keys = {'id', 'description', 'sequence', 'qualities', 'strand', 'true_start', 'true_end', 'errors'}
    for read in reads:
        assert expected_keys.issubset(read.keys())


def test_pacbio_quality_and_sequence_same_length(pacbio_simulator, short_sequence):
    reads = pacbio_simulator.generate_reads("test_seq", short_sequence)
    for read in reads:
        assert len(read["sequence"]) == len(read["qualities"])


def test_pacbio_coordinates_within_sequence_bounds(pacbio_simulator, short_sequence):
    reads = pacbio_simulator.generate_reads("test_seq", short_sequence)
    seq_len = len(short_sequence)
    for read in reads:
        assert 1 <= read["true_start"] <= seq_len
        assert read["true_end"] <= seq_len
        assert read["true_start"] < read["true_end"]


def test_pacbio_coverage_produces_expected_read_count(pacbio_simulator):
    seq = "ATCG" * 1000
    pacbio_simulator.coverage = 10
    reads = pacbio_simulator.generate_reads("test_seq", seq)
    expected = int((4000 * 10) / 100)
    assert len(reads) == expected


def test_pacbio_read_length_exceeds_sequence_raises(pacbio_simulator):
    with pytest.raises(ValueError):
        pacbio_simulator.generate_reads("test_seq", "ATCG")


def test_pacbio_all_bases_valid(pacbio_simulator, short_sequence):
    reads = pacbio_simulator.generate_reads("test_seq", short_sequence)
    valid_bases = set('ATCGN')
    for read in reads:
        assert set(read["sequence"]).issubset(valid_bases)