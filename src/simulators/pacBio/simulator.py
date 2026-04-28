import bisect
import math
import os
from ..base_simulator import BaseSimulator
from typing import List, Dict, Tuple
import random
import time

# empirical PacBio error ratios (ono2022)
BACKGROUND_WEIGHTS = [0.22, 0.45, 0.33]
BACKGROUND_TYPES = ["substitution", "insertion", "deletion"]

# homopolymer indel ratio 95% deletion 5% insertion (karst 2021)
HOMOPOLYMER_DEL_PROB = 0.95
HOMOPOLYMER_MIN_LEN = 4

# TODO verify quality floor from lit
PACBIO_QUALITY_FLOOR = 2

HOMOPOLYMER_DECAY_CONST = 2 # TODO check (same as illumina)

# defaault log-normal CV (coeff of variance)
DEFAULT_CV = None
DEFAULT_MEAN_READ_LEN = None

def generate_quality_scores(
        read_seq: str, start_pos: int, end_pos: int, strand: str,
        profile: Dict, homopolymer_pos: Dict[int, int]) -> List[int]:
    qualities = []

    for idx, base in enumerate(read_seq):
        if strand == 'F':
            genome_pos = start_pos + idx
        else:
            # For reverse reads, index 0 is the end of the segment
            genome_pos = end_pos - 1 - idx

        # Added some robustness for when a profile is not complete in some way (found this by using N)
        if idx not in profile or base not in profile[idx]:
            scores, cumulative, total = [], [], 0
        else:
            scores = profile[idx][base]["scores"]
            cumulative = profile[idx][base]["cumulative"]
            total = profile[idx][base]["total"]

        if total == 0:
            q_baseline = 38
        else:
            r = random.uniform(0, total - 1)
            score_idx = min(bisect.bisect_left(cumulative, r), len(scores) - 1)
            q_baseline = scores[score_idx]

        q = q_baseline

        # Homopolymer same decay form as illumina need to TODO check const
        run_length = homopolymer_pos.get(genome_pos, 0)
        if run_length > 1:
            penalty = (run_length - 1) * HOMOPOLYMER_DECAY_CONST
            q = q - penalty

        qualities.append(max(PACBIO_QUALITY_FLOOR, q))
    return qualities

def inject_errors(
        sequence: str, quality_scores: List[int], start_pos: int,homopolymer_pos: Dict[int, int]) \
        -> Tuple[str, List[int], List[Dict]]:

    out_seq = []
    out_qual = []
    errors_injected = []
    # WE have to append a lot in this because of insertions and deletions, let elegant than Illumina

    for i, base in enumerate(sequence):
        error_probability = 10 ** (-quality_scores[i] / 10.0)

        if random.random() >= error_probability:
            out_seq.append(base)
            out_qual.append(quality_scores[i])
            continue

        genome_pos = start_pos + i
        run_length = homopolymer_pos.get(genome_pos, 0)

        if run_length >= HOMOPOLYMER_MIN_LEN:
            # homopolymer indels
            if random.random() < HOMOPOLYMER_DEL_PROB:
                # Deletion, omit base and quality
                errors_injected.append(
                {
                    'position': i,
                    'type': 'deletion',
                    'mode': 'homopolymer',
                    'original_base': base,
                    'mutated': None,
                    'quality': quality_scores[i]
                })
            else:
                # insertion, append extra copy of original
                out_seq.append(base)
                out_qual.append(quality_scores[i])
                out_seq.append(base)
                out_qual.append(quality_scores[i])
                errors_injected.append(
                    {
                        'position': i,
                        'type': 'insertion',
                        'mode': 'homopolymer',
                        'original_base': base,
                        'mutated': base,
                        'quality': quality_scores[i]
                    })
        else:
            # Background errors 22:45:33 subs:ins:del
            error_type = random.choices(BACKGROUND_TYPES, weights = BACKGROUND_WEIGHTS, k=1)[0]

            if error_type == 'substitution':
                alternative_bases = [b for b in ['A', 'C', 'G', 'T'] if b != base]
                mutated = random.choice(alternative_bases)
                out_seq.append(mutated)
                out_qual.append(quality_scores[i])
                errors_injected.append(
                    {
                        'position': i,
                        'type': 'substitution',
                        'mode': 'background',
                        'original_base': base,
                        'mutated': mutated,
                        'quality': quality_scores[i]
                    })
            elif error_type == 'insertion':
                inserted = random.choice(['A', 'T', 'G', 'C'])
                out_seq.append(inserted)
                out_qual.append(quality_scores[i])
                out_seq.append(base)
                out_qual.append(quality_scores[i])
                errors_injected.append(
                    {
                        'position': i,
                        'type': 'insertion',
                        'mode': 'background',
                        'original_base': base,
                        'mutated': inserted,
                        'quality': quality_scores[i]
                    })
            else:
                # Deletion omit base
                errors_injected.append(
                    {
                        'position': i,
                        'type': 'deletion',
                        'mode': 'background',
                        'original_base': base,
                        'mutated': None,
                        'quality': quality_scores[i]
                    })

    return ''.join(out_seq), out_qual, errors_injected

def _simulate_batch(
    batch_indices,
    sequence_header,
    sequence,
    mean_read_length,
    cv,
    coverage,
    platform_name,
    profile,
    homopolymer_pos,
    reverse_complement_fn):

    seq_len = len(sequence)
    reads = []

    # derived log-normal parameters from mean and CV
    s = cv * mean_read_length
    mu = math.log(mean_read_length ** 2 / math.sqrt(mean_read_length ** 2 + s ** 2))
    sigma2 = math.log(1 + (s ** 2 / mean_read_length ** 2))
    sigma = math.sqrt(sigma2)

    for i in batch_indices:
        # Sample read length from log-normal
        read_length = max(1, int(random.lognormvariate(mu, sigma)))
        read_length = min(read_length, seq_len)

        start_pos = random.randint(0, seq_len - read_length)
        end_pos = start_pos + read_length
        strand = random.choice(['F', 'R'])

        read_seq = sequence[start_pos:end_pos]
        if strand == 'R':
            read_seq = reverse_complement_fn(read_seq)

        quality_scores = generate_quality_scores(
            read_seq, start_pos, end_pos, strand,
            profile, homopolymer_pos
        )
        read_seq_with_errors, final_qualities, errors_injected = inject_errors(
            read_seq, quality_scores, start_pos, homopolymer_pos
        )

        reads.append({
            "id": f"{sequence_header}_read_{i + 1}",
            "description": (
                f"simulated_read platform={platform_name} "
                f"start={start_pos + 1} end={end_pos} "
                f"length={len(read_seq_with_errors)} coverage={coverage}"
            ),
            "sequence": read_seq_with_errors,
            "qualities": final_qualities,
            "true_start": start_pos + 1,
            "true_end": end_pos,
            "strand": strand,
            "errors": errors_injected,
        })

    return reads

class PacBioSimulator(BaseSimulator):
    def __init__(self, mean_read_length: int, coverage: float, profile_filename: str,
                 cv: float = DEFAULT_CV, num_workers: int = os.cpu_count(), batch_size: int = 5000):
        """
        Initialise simulator with the PacBio profile.

        :param:
            mean_read_length: Target length of reads to generate
            coverage:         Desired coverage depth
            profile_filename: File containing the quality profile
            num_workers:      Number of worker processes (default: os.cpu_count())
            batch_size:       Reads per worker batch (default: 5000)
        """
        super().__init__(mean_read_length, coverage, profile_filename, num_workers, batch_size)
        self.mean_read_length = mean_read_length
        self.cv = cv

    def get_platform_name(self) -> str:
        return "PacBio HiFi"

    def _precompute(self, sequence):
        start = time.time()
        self.homopolymer_pos = self.find_homopolymer_runs(
            sequence, min_length=HOMOPOLYMER_MIN_LEN
        )
        print(f"PacBio Precomputation took {time.time() - start:.2f}s")

    def _get_worker_fn(self):
        return _simulate_batch

    def _get_worker_args(self, batch, sequence_header, sequence) -> tuple:
        return (
            batch,
            sequence_header,
            sequence,
            self.mean_read_length,
            self.cv,
            self.coverage,
            self.get_platform_name(),
            self.profile,
            self.homopolymer_pos,
            self.reverse_complement,
        )

    def generate_quality_scores(self, read_seq, start_pos, end_pos, strand):
        return generate_quality_scores(
            read_seq, start_pos, end_pos, strand, self.profile, self.homopolymer_pos
        )

    def inject_errors(self, sequence, quality_scores):
        return inject_errors(
            sequence, quality_scores, self.homopolymer_pos, 0
        )