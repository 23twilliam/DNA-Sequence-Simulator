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

def test_empty_sequence_no_triggers(simulator):
    trigger_map = simulator.build_sse_trigger_map("AAAAAAAAAAAAAAAAAAAAAA")
    assert len(trigger_map) == 0

def test_trigger_map_returns_dict(simulator):
    trigger_map = simulator.build_sse_trigger_map("ATCGATCGATCG")
    assert isinstance(trigger_map, dict)

def test_ggc_marks_forward_positions(simulator):
    # GGC at position 0 should mark forward positions 1-40
    seq = "GGC" + "A" * 50
    trigger_map = simulator.build_sse_trigger_map(seq)
    assert trigger_map.get(1) == 'F'
    assert trigger_map.get(40) == 'F'

def test_ggc_does_not_mark_before_trigger(simulator):
    seq = "A" * 10 + "GGC" + "A" * 50
    trigger_map = simulator.build_sse_trigger_map(seq)
    # Positions before index 10 should not be marked F from this trigger
    assert trigger_map.get(0) != 'F'
    assert trigger_map.get(9) != 'F'

def test_ggc_impact_area_boundary(simulator):
    # Impact area is 40 bases, so position 41 should not be marked
    seq = "GGC" + "A" * 60
    trigger_map = simulator.build_sse_trigger_map(seq)
    assert 40 in trigger_map  # position 40 should be marked
    assert 41 not in trigger_map  # position 41 shouldn't

def test_ggc_marks_reverse_positions(simulator):
    # GCC is GGC on reverse strand
    seq = "A" * 50 + "GCC" + "A" * 10
    trigger_map = simulator.build_sse_trigger_map(seq)
    # Should have some R strand marks BEFORE the GCC
    r_positions = [pos for pos, strand in trigger_map.items() if strand == 'R']
    assert len(r_positions) > 0

def test_ggc_and_gcc_both_present(simulator):
    seq = "GGC" + "A" * 50 + "GCC"
    trigger_map = simulator.build_sse_trigger_map(seq)
    f_positions = [pos for pos, strand in trigger_map.items() if strand == 'F']
    r_positions = [pos for pos, strand in trigger_map.items() if strand == 'R']
    assert len(f_positions) > 0
    assert len(r_positions) > 0

########################################################################################################################

def test_inverted_repeat_detected(simulator):
    # kmer and its reverse complement within window
    kmer = "ATCGATCGA"  # 9bp (min to be considered)
    rc_kmer = simulator.reverse_complement(kmer)
    # Place kmer at 0, reverse complement within 40bp window
    seq = kmer + "A" * 5 + rc_kmer + "A" * 100
    trigger_map = simulator.build_sse_trigger_map(seq)
    assert len(trigger_map) > 0

def test_inverted_repeat_too_far_not_detected(simulator):
    # reverse complement kmer placed beyond 40bp window should not trigger
    kmer = "ATCGATCGA"
    rc_kmer = simulator.reverse_complement(kmer)
    # Place kmer at 0, reverse complement 50bp away (beyond window)
    seq = kmer + "A" * 50 + rc_kmer + "A" * 100
    trigger_map = simulator.build_sse_trigger_map(seq)

    seq_far = kmer + "A" * 100 + rc_kmer + "A" * 100
    trigger_map_far = simulator.build_sse_trigger_map(seq_far)
    # Far placement should result in fewer or different trigger positions
    assert len(trigger_map_far) <= len(trigger_map)
