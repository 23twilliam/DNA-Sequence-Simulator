import pytest
from src.simulators.illumina.simulator import inject_errors, weighted_substitution

def test_injection_errors_output_types():
    seq = "ATCGATCGATCG"
    qualities = [30] * 12

    result_seq, result_qualities, errors = inject_errors(
        sequence=seq, quality_scores=qualities, start_pos=0, strand='F', is_dephased=False, sse_trigger_map={}, homopolymer_pos={}
    )
    assert isinstance(result_seq, str)
    assert isinstance(result_qualities, list)
    assert isinstance(errors, list)

def test_injection_errors_high_quality_produces_few_errors():
    # Q40 = 0.01% error rate, very few errors expected
    seq = "ATCGATCGAT" * 10
    qualities = [40] * 100
    _, _, errors = inject_errors(
        sequence=seq, quality_scores=qualities, start_pos=0, strand='F', is_dephased=False, sse_trigger_map={}, homopolymer_pos={}
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
        _, _, errors_h = inject_errors(sequence=seq, quality_scores=high_quals, start_pos=0, strand='F', is_dephased=False, sse_trigger_map={}, homopolymer_pos={})
        _, _, errors_l = inject_errors(sequence=seq, quality_scores=low_quals, start_pos=0, strand='F', is_dephased=False, sse_trigger_map={}, homopolymer_pos={})
        high_q_errors.append(len(errors_h))
        low_q_errors.append(len(errors_l))

    assert sum(low_q_errors) > sum(high_q_errors)

def test_error_record_structure():
    # Try to force error using a low score
    seq = "ATCGATCGATCG"
    qualities = [2] * 12 # Low score to make error very likely
    errors_found = []
    for _ in range(100):
        _, _, errors = inject_errors(
        sequence=seq, quality_scores=qualities, start_pos=0, strand='F', is_dephased=False, sse_trigger_map={}, homopolymer_pos={}
        )
        errors_found.extend(errors)
        if errors_found:
            break

    if errors_found:
        err = errors_found[0]
        assert 'position' in err
        assert 'type' in err
        assert 'mode' in err
        assert 'original' in err
        assert 'mutated' in err
        assert 'quality' in err
        assert err['type'] == 'substitution'
        assert err['mode'] in ['background', 'lagging', 'prephasing']

def test_error_mode_background_when_not_depahsed():
    # When is_dephased is false, all errors must be background errors
    seq = "ATCGATCGATCG"
    qualities = [2] * 12 # Low score to make error very likely

    errors_found = []
    for _ in range(100):
        _, _, errors = inject_errors(
            sequence=seq, quality_scores=qualities, start_pos=0, strand='F', is_dephased=False,
            sse_trigger_map={i: 'F' for i in range(12)},
            homopolymer_pos={i: 5 for i in range(12)}
        )
        errors_found.extend(errors)

    for error in errors_found:
        assert error['mode'] == 'background'

def test_depahsed_can_produce_lagging_errors():
    seq = "ATCGATCGATCG"
    qualities = [2] * 12

    modes_seen = set()
    for _ in range(500):
        _, _, errors = inject_errors(
            sequence=seq, quality_scores=qualities, start_pos=0, strand='F', is_dephased=True,
            sse_trigger_map={i: (i, 'F') for i in range(12)},
            homopolymer_pos={}
        )
        for error in errors:
            modes_seen.add(error['mode'])

    # With SSE triggers and dephasing, lagging mode should appear
    assert 'lagging' in modes_seen

def test_prephasing_trim_occurs():
    # When 3+ prephasing errors, sequence is trimmed, this is hard to force so we run many times
    seq = "ATATATATATAT" # alternate bases so that next_base != original_base is always true (except first base)
    qualities = [2] * 12

    trimmed_count = 0
    for _ in range(100):
        result_seq, result_quals, errors = inject_errors(
            sequence=seq, quality_scores=qualities, start_pos=0, strand='F', is_dephased=True,
            sse_trigger_map={},
            # All positions in a homopolymer run of length 4 to trigger prephasing
            homopolymer_pos={i: 4 for i in range(12)},
        )
        if len(result_seq) < len(seq):
            trimmed_count += 1

    # Given Q2 quality and run_length=4 (p=0.84 prephasing),
    # expect ~50% chance of prephasing per position.
    # With 12 bp, getting 3+ prephasing errors is almost certain.
    # Trim should occur in the vast majority of runs. Let's say over 90% to allow for some margin of error
    assert trimmed_count > 90

def test_mutated_base_differs_from_original():
    seq = "AAAAAAAAAA"
    qualities = [2] * 10

    for _ in range(100):
        _, _, errors = inject_errors(
            sequence=seq, quality_scores=qualities, start_pos=0, strand='F', is_dephased=True,sse_trigger_map={},homopolymer_pos={}
        )
        for error in errors:
            assert error['mutated'] != error['original']

########################################################################################################################

# Weighted Substitution Tests

def test_weighted_substitution_never_return_original():
    for base in ['A', 'C', 'G', 'T']:
        for _ in range(100):
            result = weighted_substitution(base)
            assert result != base

def test_weighted_substitution_returns_valid_base():
    for base in ['A', 'C', 'G', 'T']:
        result = weighted_substitution(base)
        assert result in ['A', 'C', 'G', 'T']

def test_weighted_substitution_distribution():
    # e.g. C->A should be most common for C (weight 0.13 vs 0.11 T, 0.04 G)
    # Run many times to see pattern
    results = [weighted_substitution('C') for _ in range(10000)]
    a_count = results.count('A')
    g_count = results.count('G')
    t_count = results.count('T')

    # A should be most common for C (weight 0.13), T second (weight 0.11), G third (weight 0.04)
    assert  a_count > t_count > g_count