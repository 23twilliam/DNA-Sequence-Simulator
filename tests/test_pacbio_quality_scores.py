import pytest
from src.simulators.pacbio.simulator import generate_quality_scores, PacBioSimulator
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
        seq, start_pos=0, end_pos=len(seq) ,strand='F',
        profile=profile, homopolymer_pos={}
    )
    assert len(scores) == len(seq)

def test_quality_scores_minimum_floor():
    # Scores should never go below 2
    profile = make_test_profile(read_length=12)
    seq = "ATCGATCGATCG"
    scores = generate_quality_scores(
        seq, start_pos=0, end_pos=len(seq), strand='F',
        profile=profile,
        homopolymer_pos={i: 50 for i in range(12)}
    )
    assert all(score >= 7 for score in scores)

def test_missing_profile_falls_back_to_baseline():
    # should default to 38 when no profile
    seq = "ATCGATCGATCG"
    scores = generate_quality_scores(
        seq, start_pos=0, end_pos=len(seq),  strand='F', profile={}, homopolymer_pos= {}
    )
    assert all(score == 38 for score in scores)
