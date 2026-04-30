import bisect
import math
import os
from ..base_simulator import BaseSimulator
from typing import List, Dict, Tuple
import random
import time

# Park et al. Fig 1A
STR_WEIGHTS_BY_LENGTH = {
    1: [0.05, 0.17, 0.78], # sub ins del
    2: [0.24, 0.30, 0.46],
    3: [0.42, 0.27, 0.31]
}

BACKGROUND_WEIGHTS = [0.31, 0.30, 0.38]
BACKGROUND_TYPES = ["substitution", "insertion", "deletion"]

# A <-> G and C <-> T are 3-5 times more likely, using 4 as an average
SUBST_CHOICES: Dict[str, Tuple[List[str], List[float]]] = {
    "A": (["C", "G", "T"], [1, 4, 1]),
    "C": (["A", "G", "T"], [4, 1, 1]),
    "G": (["A", "C", "T"], [1, 1, 4]),
    "T": (["A", "C", "G"], [1, 4, 1]),
    "N": (["N"], [1.0])
}

DEFAULT_CV = 0.2
DEFAULT_MEAN_READ_LEN = 20_000

# TODO verify quality floor from lit
NANOPORE_QUALITY_FLOOR = 2

HOMOPOLYMER_DECAY_CONST = 2 # TODO check (same as illumina)
HOMOPOLYMER_MIN_LEN = 3

END_PENALTY_LEN = 100
END_PENALTY_MAX = 10  # TODO check vs Delahaye

# Park et al.
STR_ACCURATE_QUALITY = {
    1: 7,   # 19% per base error at 1bp
    2: 11,  # 7% per base error at 2bp
    3: 14   # 4% per base error at 3bp
}

STR_OVERREPORTED_QUALITY = {
    1: 75,   # Median 70-80 fig 6A
    2: 90,   # Basecaller pins at fixed Q90
    3: 75,   # Median 70-80 fig 6A
}

def generate_quality_scores(
        read_seq: str, start_pos: int, end_pos: int, strand: str,
        profile: Dict, corrected_q_map: Dict[int, int]) -> List[int]:
    qualities = []
    read_len = len(read_seq)

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

        # end of read quality penalty delahaye et al.
        if idx < END_PENALTY_LEN:
            q -= int((END_PENALTY_LEN - idx) / END_PENALTY_LEN * END_PENALTY_MAX)
        elif idx >= read_len - END_PENALTY_LEN:
            distance_from_end = idx - (read_len - END_PENALTY_LEN)
            q -= int(distance_from_end / END_PENALTY_LEN * END_PENALTY_MAX)

        #STR override, replace baseline with corrected (low) quality scores
        if genome_pos in corrected_q_map:
            q = corrected_q_map[genome_pos]

        qualities.append(max(NANOPORE_QUALITY_FLOOR, q))
    return qualities

def inject_errors(
        sequence: str, quality_scores: List[int], start_pos: int, repeat_map: Dict[int, Dict], overreported_q_map: Dict[int, int]) -> Tuple[str, List[int], List[Dict]]:

    out_seq = []
    out_qual = []
    errors_injected = []
    # WE have to append a lot in this because of insertions and deletions, less elegant than Illumina
    out_genome_pos = []

    for i, base in enumerate(sequence):
        genome_pos = start_pos + i
        error_probability = 10 ** (-quality_scores[i] / 10.0)

        if random.random() >= error_probability:
            out_seq.append(base)
            out_qual.append(quality_scores[i])
            out_genome_pos.append(genome_pos)
            continue

        genome_pos = start_pos + i
        str_entry = repeat_map.get(genome_pos)
        if str_entry:
            weights = STR_WEIGHTS_BY_LENGTH[str_entry["motif_len"]]
            mode = "str"
        else:
            weights = BACKGROUND_WEIGHTS
            mode = "background"

        error_type = random.choices(BACKGROUND_TYPES, weights=weights, k=1)[0]

        if error_type == 'substitution':
            if mode == 'background':
                mutated_base = weighted_substitution(base)
            else:
                alternative_bases = [b for b in ['A', 'C', 'G', 'T'] if b != base]
                mutated_base = random.choice(alternative_bases)
            out_seq.append(mutated_base)
            out_qual.append(quality_scores[i])
            out_genome_pos.append(genome_pos)
            errors_injected.append(
                {
                    'position': i,
                    'type': 'substitution',
                    'mode': mode,
                    'original_base': base,
                    'mutated_base': mutated_base,
                    'quality': quality_scores[i]
                })
        elif error_type == 'insertion':
            # STR insertions duplicate current base, background gets random
            if mode == 'str':
                inserted = base
            else:
                alternative_bases = [b for b in ['A', 'C', 'G', 'T']]
                inserted = random.choice(alternative_bases)
            out_seq.append(inserted)
            out_qual.append(quality_scores[i])
            out_genome_pos.append(genome_pos)
            out_seq.append(base)
            out_qual.append(quality_scores[i])
            out_genome_pos.append(genome_pos)
            errors_injected.append(
                {
                    'position': i,
                    'type': 'insertion',
                    'mode': mode,
                    'original_base': base,
                    'mutated_base': inserted,
                    'quality': quality_scores[i]
                })
        else:
            # Deletion omit base
            errors_injected.append(
                {
                    'position': i,
                    'type': 'deletion',
                    'mode': mode,
                    'original_base': base,
                    'mutated_base': None,
                    'quality': quality_scores[i]
                })

    # Sqap STR positions to the over reported quality score for the output
    final_qualities = [
        overreported_q_map.get(gp, q)
        for q, gp in zip(out_qual, out_genome_pos) # TODO make more robust zip can cause error
    ]

    return ''.join(out_seq), final_qualities, errors_injected

def _simulate_batch(
    batch_indices,
    sequence_header,
    sequence,
    mean_read_length,
    cv,
    coverage,
    platform_name,
    profile,
    repeat_map,
    quality_override_map,
    quality_overreport_map,
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
            profile, quality_override_map
        )
        read_seq_with_errors, final_qualities, errors_injected = inject_errors(
            read_seq, quality_scores, start_pos, repeat_map, quality_overreport_map
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

def weighted_substitution(original: str) -> str:
    keys, weights = SUBST_CHOICES[original]
    return random.choices(keys, weights=weights, k=1)[0]

class NanoporeSimulator(BaseSimulator):
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
        return "Oxford Nanopore"

    def _precompute(self, sequence):
        start = time.time()
        homopolymers = self.find_homopolymer_runs(
            sequence, min_length=HOMOPOLYMER_MIN_LEN
        )
        # convert from dict[int, int] so we can merge properly
        homopolymers = {
            pos: {"motif": sequence[pos], "motif_len":1, "run_len": run_len}
            for pos, run_len in homopolymers.items()
        }
        dinucs = self.find_repeats(sequence, 2)
        trinucs = self.find_repeats(sequence, 3)
        self.repeat_map = self.merge_repeat_maps(homopolymers, dinucs, trinucs)
        self.quality_override_map = self.build_corrected_quality_map(self.repeat_map)
        self.quality_overreport_map = self.build_overreported_quality_map(self.repeat_map)
        print(f"Nanopore Precomputation took {time.time() - start:.2f}s")

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
            self.repeat_map,
            self.quality_override_map,
            self.quality_overreport_map,
            self.reverse_complement,
        )

    def find_repeats(self, sequence: str, motif_len: int, min_repeats: int = 2) -> Dict[int, Dict]:
        """Find tandem repeats of given motif length. Returns {pos: {motif, length, run_len}}"""
        result = {}
        seq_len = len(sequence)
        i = 0
        while i <= seq_len - motif_len * min_repeats:
            motif = sequence[i:i + motif_len]
            # Count consecutive non overlapping occurrences
            count = 1
            j = i + motif_len
            while j + motif_len <= seq_len and sequence[j:j + motif_len] == motif:
                count += 1
                j += motif_len
            if count >= min_repeats:
                run_len = count * motif_len
                for k in range(i, i + run_len):
                    result[k] = {"motif": motif, "motif_len": motif_len, "run_len": run_len}
                i += run_len # skip past
            else:
                i += 1
        return result

    def merge_repeat_maps(self, *maps: Dict[int, Dict]) -> Dict[int, Dict]:
        merged = {}
        for m in maps:
            for pos, entry in m.items():
                if pos not in merged or entry["motif_len"] > merged[pos]["motif_len"]:
                    merged[pos] = entry
        return merged

    def build_corrected_quality_map(self, repeat_map: Dict[int, Dict]) -> Dict[int, int]:
        """Map of position : low quality (reflecting true error rate at STRs)"""
        return {
            pos: STR_ACCURATE_QUALITY[entry["motif_len"]]
            for pos, entry in repeat_map.items()
        }

    def build_overreported_quality_map(self, repeat_map: Dict[int, Dict]) -> Dict[int, int]:
        """Map of position : high quality (reflecting base caller over reporting at STRs)"""
        return {
            pos: STR_OVERREPORTED_QUALITY[entry["motif_len"]]
            for pos, entry in repeat_map.items()
        }