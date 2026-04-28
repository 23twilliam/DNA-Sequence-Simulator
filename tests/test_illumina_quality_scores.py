import pytest
from src.simulators.illumina.simulator import generate_quality_scores, IlluminaSimulator
from pathlib import Path

TEST_DIR = Path(__file__).resolve().parent
FILES_DIR = TEST_DIR / "test_files"

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

def test_quality_scores_correct_length():
    profile = make_test_profile(read_length=12)
    seq = "ATCGATCGATCG"
    scores = generate_quality_scores(
        seq, start_pos=0, end_pos=len(seq) ,strand='F', is_dephased=False,
        profile=profile, sse_trigger_map={}, homopolymer_pos={}
    )
    assert len(scores) == len(seq)

def test_quality_scores_minimum_floor():
    # Scores should never go below 2
    profile = make_test_profile(read_length=12)
    seq = "ATCGATCGATCG"
    scores = generate_quality_scores(
        seq, start_pos=0, end_pos=len(seq), strand='F', is_dephased=False,
        profile=profile,
        sse_trigger_map={i: (i, 'F') for i in range(12)},
        homopolymer_pos={i: 5 for i in range(12)}
    )
    assert all(score >= 2 for score in scores)

def test_sse_penalty_reduces_quality():
    profile = make_test_profile(read_length=12)
    seq = "ATCGATCGATCG"

    scores_no_see = generate_quality_scores(
        seq, start_pos=0, end_pos=len(seq),  strand='F', is_dephased=True, profile=profile, sse_trigger_map={}, homopolymer_pos= {}
    )
    scores_with_see = generate_quality_scores(
        seq, start_pos=0, end_pos=len(seq), strand='F', is_dephased=True, profile=profile,
        sse_trigger_map={i: (i, 'F') for i in range(12)},
        homopolymer_pos={}
    )
    # SSE should reduce scores
    assert sum(scores_with_see) < sum(scores_no_see)

def test_homopolymer_penalty_reduces_quality():
    profile = make_test_profile(read_length=12)
    seq = "ATCGATCGATCG"

    scores_no_homopolymer = generate_quality_scores(
        seq, start_pos=0, end_pos=len(seq),  strand='F', is_dephased=True, profile=profile, sse_trigger_map={}, homopolymer_pos= {}
    )
    scores_with_homopolymer = generate_quality_scores(
        seq, start_pos=0, end_pos=len(seq), strand='F', is_dephased=True, profile=profile,
        sse_trigger_map={},
        homopolymer_pos={i: 5 for i in range (12)}
    )
    # Homopolymer should reduce scores
    assert sum(scores_with_homopolymer) < sum(scores_no_homopolymer)

def test_strand_specificity_sse():
    # SSE marked 'R' should not penalise 'F' strand reads
    profile = make_test_profile(read_length=12)
    seq = "ATCGATCGATCG"

    scores_f = generate_quality_scores(
        seq, start_pos=0, end_pos=len(seq), strand='F', is_dephased=True, profile=profile,
        sse_trigger_map={i: (i, 'R') for i in range(12)},
        homopolymer_pos={})
    scores_no_sse = generate_quality_scores(
        seq, start_pos=0, end_pos=len(seq), strand='F', is_dephased=True, profile=profile, sse_trigger_map={}, homopolymer_pos={}
    )

    # R-strand SSE should not affect F-strand quality
    assert scores_f == scores_no_sse

def test_load_profile_raises_on_malformed_line(tmp_path):
    profile_file = tmp_path / "bad_profile.txt"
    profile_file.write_text("A\tNOT_AN_INT\t30\nA\t0\t100\n")

    with pytest.raises(ValueError, match="Improper profile entry at line 1"):
        IlluminaSimulator.load_quality_profile(self = None, filename= str(profile_file))


def test_load_profile_valid_file(tmp_path):
    profile_file = tmp_path / "profile.txt"
    profile_file.write_text("A\t0\t30\t31\nA\t0\t100\t200\n")

    profile = IlluminaSimulator.load_quality_profile(self = None, filename= str(profile_file))
    assert profile[0]["A"]["total"] == 200

def test_profile(tmp_path):
    profile_file = FILES_DIR / "test_profile.txt"
    profile = IlluminaSimulator.load_quality_profile(self = None, filename= str(profile_file))
    print(profile)
    assert profile[0]["A"]["total"] == 4