import pytest
from src.simulators.pacbio.simulator import inject_errors

def test_injection_errors_output_types():
    seq = "ATCGATCGATCG"
    qualities = [30] * 12

    result_seq, result_qualities, errors = inject_errors(
        sequence=seq, quality_scores=qualities, start_pos=0, homopolymer_pos={}
    )
    assert isinstance(result_seq, str)
    assert isinstance(result_qualities, list)
    assert isinstance(errors, list)

def test_injection_errors_high_quality_produces_few_errors():
    # Q40 = 0.01% error rate, very few errors expected
    seq = "ATCGATCGAT" * 10
    qualities = [40] * 100
    _, _, errors = inject_errors(
        sequence=seq, quality_scores=qualities, start_pos=0, homopolymer_pos={}
    )
    # With 0.01% error rate, expect very few errors in 100bp (Should be 1 on avg but go higher so test is extremely unlikely to fall over)
    assert len(errors) < 5

def test_injection_errors_low_quality_produces_more_errors():
    seq = "ATCGATCGAT" * 10 # 100 bp
    high_quals = [40] * 100
    low_quals = [10] * 100  # Q10 = 10% error rate

    high_q_errors = []
    low_q_errors = []

    # Run many times so that we can remove issues with randomness
    for _ in range(100):
        _, _, errors_h = inject_errors(sequence=seq, quality_scores=high_quals, start_pos=0, homopolymer_pos={})
        _, _, errors_l = inject_errors(sequence=seq, quality_scores=low_quals, start_pos=0, homopolymer_pos={})
        high_q_errors.append(len(errors_h))
        low_q_errors.append(len(errors_l))

    assert sum(low_q_errors) > sum(high_q_errors)


def test_error_record_structure():
    # Try to force error using a low score
    seq = "ATCGATCGATCG"
    qualities = [2] * 12  # Low score to make error very likely
    errors_found = []
    for _ in range(100):
        _, _, errors = inject_errors(
            sequence=seq, quality_scores=qualities, start_pos=0, homopolymer_pos={}
        )
        errors_found.extend(errors)
        if errors_found:
            break

    if errors_found:
        err = errors_found[0]
        assert 'position' in err
        assert 'type' in err
        assert 'mode' in err
        assert 'original_base' in err
        assert 'mutated_base' in err
        assert 'quality' in err
        assert err['type'] in ['substitution', 'insertion', 'deletion']
        assert err['mode'] in ['background', 'lagging', 'prephasing']

def test_output_length_matches_indel_count():
    seq = "ATCGATCG" * 5
    qualities = [2] * 40

    for _ in range(40):
        result_seq, result_quals, errors = inject_errors(
            sequence=seq, quality_scores=qualities, start_pos=0, homopolymer_pos={}
        )
    insertions = sum(1 for error in errors if error['type'] == 'insertion')
    deletions = sum(1 for error in errors if error['type'] == 'deletion')
    expected_len = len(seq) + insertions - deletions

    assert len(result_seq) == expected_len
    assert len(result_quals) == expected_len

def test_output_length_stays_same_with_homopolymers():
    # Same test but when homopolymer indels are involved
    seq = "AAAAAAAAAA"
    qualities = [2] * 10
    homopolymer_pos = {i: 10 - i for i in range(10)}

    for _ in range(100):
        result_seq, result_quals, errors = inject_errors(
            sequence=seq, quality_scores=qualities, start_pos=0, homopolymer_pos=homopolymer_pos
        )

        insertions = sum(1 for error in errors if error['type'] == 'insertion')
        deletions = sum(1 for error in errors if error['type'] == 'deletion')
        expected_len = len(seq) + insertions - deletions

        assert len(result_seq) == expected_len
        assert len(result_quals) == expected_len

def test_all_three_background_error_types_appear():
    # With no homopolymers, all three types should appear
    seq = "ATCGATCG" * 5
    qualities = [2] * 40

    types_seen = set()
    for _ in range(100):
        _, _, errors = inject_errors(
            sequence=seq, quality_scores=qualities, start_pos=0, homopolymer_pos={}
        )
        for e in errors:
            types_seen.add(e['type'])

    assert {'substitution', 'insertion', 'deletion'}.issubset(types_seen)

def test_substitution_mutated_base_different():
    # Same test but when homopolymer indels are involved
    seq = "AAAAAAAAAA"
    qualities = [2] * 10
    homopolymer_pos = {i: 10 - i for i in range(10)}

    for _ in range(100):
        _, _, errors = inject_errors(
            sequence=seq, quality_scores=qualities, start_pos=0, homopolymer_pos=homopolymer_pos
        )

        for error in errors:
            if error['type'] == 'substitution':
                assert error['mutated_base'] != error['original_base']

def test_homopolymer_mostly_deletions():
    seq = "AAAAAAAAAA"
    qualities = [2] * 10
    homopolymer_pos = {i: 10 - i for i in range(10)}

    deletions = 0
    insertions = 0
    for _ in range(10000):
        _, _, errors = inject_errors(
            sequence=seq, quality_scores=qualities, start_pos=0, homopolymer_pos=homopolymer_pos
        )
        for error in errors:
            if error['mode'] == 'homopolymer':
                if error['type'] == 'deletion':
                    deletions += 1
                elif error['type'] == 'insertion':
                    insertions += 1

    total_indels = deletions + insertions
    assert total_indels > 0
    assert deletions / total_indels == pytest.approx(0.95, abs=0.01) # some leeway

def test_homopolymer_min_length():
    seq = "AAAAAAAAAA"
    qualities = [2] * 10

    # mark positions with run lengths below 4
    homopolymer_pos = {i: 3 for i in range(10)}

    homopolymer_triggers = 0
    for _ in range(100):
        _, _, errors = inject_errors(
            sequence=seq, quality_scores=qualities, start_pos=0, homopolymer_pos=homopolymer_pos
        )
        homopolymer_triggers += sum(1 for error in errors if error['mode'] == 'homopolymer')

    assert homopolymer_triggers == 0

def test_homopolymer_uses_genome_coordinates():
    seq = "AAAAAAAAAA"
    qualities = [2] * 10
    # Mark genome positions 100-109 as homopolymer, read start at 100
    homopolymer_pos = {100 + i: 10 - i for i in range(10)}

    homopolymer_triggers = 0
    for _ in range(100):
        _, _, errors = inject_errors(
            sequence=seq, quality_scores=qualities, start_pos=100, homopolymer_pos=homopolymer_pos
        )
        homopolymer_triggers += sum(1 for error in errors if error['mode'] == 'homopolymer')

    assert homopolymer_triggers > 0

def test_background_error_ratio_approx():
    seq = "ATCGATCG" * 25
    qualities = [2] * 200

    deletions = 0
    insertions = 0
    substitutions = 0

    for _ in range(10000):
        _, _, errors = inject_errors(
            sequence=seq, quality_scores=qualities, start_pos=0, homopolymer_pos={}
        )
        for error in errors:
            if error['type'] == 'deletion':
                deletions += 1
            elif error['type'] == 'insertion':
                insertions += 1
            elif error['type'] == 'substitution':
                substitutions += 1

    total_indels = deletions + insertions + substitutions

    assert deletions / total_indels == pytest.approx(0.33, abs=0.01)
    assert insertions / total_indels == pytest.approx(0.45, abs=0.01)
    assert substitutions / total_indels == pytest.approx(0.22, abs=0.01)