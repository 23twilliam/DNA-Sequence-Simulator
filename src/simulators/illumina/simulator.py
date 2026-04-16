import bisect
from ..base_simulator import BaseSimulator
from typing import List, Dict, Tuple
import random
from collections import defaultdict
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

SUBST_CHOICES: Dict[str, Tuple[List[str], List[float]]] = {
    "A": (["C", "G", "T"], [0.04, 0.07, 0.07]),
    "C": (["A", "G", "T"], [0.13, 0.04, 0.11]),
    "G": (["A", "C", "T"], [0.07, 0.04, 0.11]),
    "T": (["A", "C", "G"], [0.07, 0.04, 0.07]),
    "N": (["N"], [1.0]),
}

PREPHASE_PROB = {1: 0.64, 2: 0.77, 3: 0.80, 4: 0.84}
PREPHASE_MAX = 0.86
LAGGING_PROB = 0.565 # avg of 49-64%
DEPHASING_RATE = 1.0 #0.0515
HOMOPOLYMER_DECAY = 2  # TODO Double check research for better decay rate
SSE_LOW, SSE_HIGH = 5, 12


def inject_errors(
        sequence: str, quality_scores: List[int], start_pos: int, strand: str, is_dephased: bool,
        sse_trigger_map: Dict[int, Tuple[int, str]], homopolymer_pos: Dict[int, int]) \
        -> Tuple[str, List[int], List[Dict]]:
    read_list = list(sequence)
    errors_injected = []

    for i in range(len(read_list)):
        error_probability = 10 ** (-quality_scores[i] / 10.0)
        if random.random() >= error_probability:
            continue

        original_base = read_list[i]
        mutated_base = None
        mode = "background"

        if is_dephased:
            sse_data = sse_trigger_map.get(start_pos + i)
            # Mode A Lagging
            if sse_data and sse_data[1] == strand:  # [1] is the index of the strand in the tuple
                if random.random() < LAGGING_PROB:
                    if i >= 2:
                        candidate = sequence[i - random.choice([1, 2])]
                    elif i == 1:
                        candidate = sequence[0]
                    else:
                        candidate = None
                    if candidate and candidate != original_base:
                        mutated_base = candidate
                        mode = "lagging"

            # Mode B Pre-phasing
            if mutated_base is None:
                run_len = homopolymer_pos.get(start_pos + i, 0)
                if run_len > 0 and i + 1 < len(sequence):
                    p = PREPHASE_PROB.get(run_len, PREPHASE_MAX)
                    next_base = sequence[i + 1]
                    if next_base != original_base and random.random() < p:
                        mutated_base = next_base
                        mode = "prephasing"

        # Mode C Background
        if mutated_base is None:
            mutated_base = weighted_substitution(original_base)
            mode = "background"

        read_list[i] = mutated_base
        errors_injected.append({
            'position': i,
            'type': 'substitution',
            'mode': mode,
            'original': original_base,
            'mutated': mutated_base,
            'quality': quality_scores[i],
        })

    read_seq_with_errors = ''.join(read_list)
    quality_out = list(quality_scores)

    if is_dephased:
        prephase_count = sum(1 for e in errors_injected if e["mode"] == "prephasing")
        if prephase_count >= 3:
            trim = random.randint(1, 3)
            read_seq_with_errors = read_seq_with_errors[:-trim]
            quality_out = quality_out[:-trim]

    return read_seq_with_errors, quality_out, errors_injected

def weighted_substitution(original: str) -> str:
    keys, weights = SUBST_CHOICES[original]
    return random.choices(keys, weights=weights, k=1)[0]


def generate_quality_scores(
        read_seq: str, start_pos: int, end_pos: int, strand: str, is_dephased: bool,
        profile, sse_trigger_map: Dict[int, str], homopolymer_pos: Dict[int, int]) -> List[int]:
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

        # SSE Penalty
        sse_data = sse_trigger_map.get(genome_pos)  # gets (sse_start_pos, strand)
        if is_dephased and sse_data:
            sse_start_pos, trigger_strand = sse_data

            if trigger_strand == strand:
                # get distance from the start of the SSE
                distance = genome_pos - sse_start_pos

                # Calculate rephasing progress (how close we are to not being in trigger boundary)
                # (0 at start, 1 at 40bp) basically it is fully unaffected 40bp away from the trigger. cannot be negative
                # so we use abs()
                rephase_progress = abs(distance / 40.0)
                # Get value based on weighted random score and how far away it is from trigger
                q_min = 2
                q = int(q_min + (q_baseline - q_min) * rephase_progress) # Cast to int as Q scores can only be an integer


        # Homopolymer
        run_length = homopolymer_pos.get(genome_pos, 0)
        if run_length > 1:
            penalty = (run_length - 1) * HOMOPOLYMER_DECAY
            q = q - penalty
        # Qualities should never go below 2 (Illumina only reports 2 or above)
        qualities.append(max(2, q))
    return qualities


def _simulate_batch(
    batch_indices,
    sequence_header,
    sequence,
    read_length,
    coverage,
    platform_name,
    profile,
    sse_trigger_map,
    homopolymer_pos,
    reverse_complement_fn):

    seq_len = len(sequence)
    reads = []

    for i in batch_indices:
        start_pos = random.randint(0, seq_len - read_length)
        end_pos = start_pos + read_length
        strand = random.choice(['F', 'R'])
        is_dephased = random.random() < DEPHASING_RATE

        read_seq = sequence[start_pos:end_pos]
        if strand == 'R':
            read_seq = reverse_complement_fn(read_seq)

        quality_scores = generate_quality_scores(
            read_seq, start_pos, end_pos, strand, is_dephased,
            profile, sse_trigger_map, homopolymer_pos
        )
        read_seq_with_errors, final_qualities, errors_injected = inject_errors(
            read_seq, quality_scores, start_pos, strand, is_dephased,
            sse_trigger_map, homopolymer_pos
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
            "errors": errors_injected,
        })

    return reads


class IlluminaSimulator(BaseSimulator):
    def __init__(self, read_length: int, coverage: float, profile_filename: str,
                 num_workers: int = 8, batch_size: int = 5000):
        """
        Initialise simulator with the Illumina profile.

        :param:
            read_length:      Length of reads to generate
            coverage:         Desired coverage depth
            profile_filename: File containing the quality profile
            num_workers:      Number of worker processes (default: os.cpu_count())
            batch_size:       Reads per worker batch (default: 5000)
        """
        super().__init__(read_length, coverage)
        self.profile = self.load_quality_profile(profile_filename)
        self.sse_trigger_map = None
        self.num_workers = num_workers
        self.batch_size = batch_size

    def get_platform_name(self) -> str:
        return "Illumina Simulator"

    # https://stackoverflow.com/questions/4664850/how-to-find-all-occurrences-of-a-substring

    def build_sse_trigger_map(self, sequence: str) -> Dict[int, str]:
        trigger_map = {}
        seq_len = len(sequence)
        IMPACT_ZONE = 40 # bp

        def mark_area(sse_start_pos: int, strand: str):
            """
            Marks the 40bp zone starting AFTER the trigger.
            """
            if strand == 'F':
                for j in range(sse_start_pos, min(sse_start_pos + IMPACT_ZONE, seq_len)):
                    trigger_map[j] = (sse_start_pos, 'F')
            else:
                for j in range(max(0, sse_start_pos - IMPACT_ZONE), sse_start_pos):
                    trigger_map[j] = (sse_start_pos, 'R')

        i = sequence.find("GGC")
        while i != -1:
            # GGC ends at i + 3 so the SSE position starts at i + 3
            mark_area(i + 3, 'F')
            i = sequence.find("GGC", i + 1)

        # GGC on forward = GCC on reverse.
        i = sequence.find("GCC")
        while i != -1:
            mark_area(i, 'R')
            i = sequence.find("GCC", i + 1)


        i = sequence.find("GGT")
        while i != -1:
            # GGC ends at i + 3 so the SSE position starts at i + 3
            mark_area(i + 3, 'F')
            i = sequence.find("GGT", i + 1)

        # GGT on forward = ACC on reverse.
        i = sequence.find("ACC")
        while i != -1:
            mark_area(i, 'R')
            i = sequence.find("ACC", i + 1)


        i = sequence.find("GGG")
        while i != -1:
            # GGC ends at i + 3 so the SSE position starts at i + 3
            mark_area(i + 3, 'F')
            i = sequence.find("GGG", i + 1)

        # GGG on forward = CCC on reverse.
        i = sequence.find("CCC")
        while i != -1:
            mark_area(i, 'R')
            i = sequence.find("CCC", i + 1)


        i = sequence.find("CGG")
        while i != -1:
            # GGC ends at i + 3 so the SSE position starts at i + 3
            mark_area(i + 3, 'F')
            i = sequence.find("GGG", i + 1)

        # CGG on forward = CCG on reverse.
        i = sequence.find("CCG")
        while i != -1:
            mark_area(i, 'R')
            i = sequence.find("CCC", i + 1)

        """
                Identify positions in the sequence where inverted repeats exist within a sliding window. Marks the impact area downstream of each trigger.
                Criteria for a valid trigger:
                - Repeat must be longer than 8bp (k=9 minimum)
                - Reverse complement must appear within +-40bp,
                 from the current position (40bp total window)
                - Errors are in the 40 bases downstream of the trigger position
            k: K-mer length (must be > 8bp per trigger)
            window: Maximum distance to search for reverse complement
        """
        k, window = 9, 40
        # Build list of all kmers
        kmer_index = defaultdict(list)
        for i in range(seq_len - k):
            kmer_index[sequence[i:i + k]].append(i)
        # Check each position for a nearby reverse compliment
        for i in range(seq_len - k):
            kmer = sequence[i:i + k]
            reverse_comp = self.reverse_complement(kmer)
            for j in kmer_index.get(reverse_comp, []):
                if i + k <= j <= i + k + window:
                    # For the forward trigger the SSE starts after the k-mer (i+k)
                    mark_area(i + k, 'F')
                    # For the reverse trigger the SSE starts before the k-mer (j)
                    mark_area(j, 'R')
                    break
        return trigger_map

    def generate_reads(self, sequence_header: str, sequence: str):
        start = time.time()
        self.sse_trigger_map = self.build_sse_trigger_map(sequence)
        print(f"SSE map built in {time.time() - start:.2f}s")
        self.homopolymer_pos = self.find_homopolymer_runs(sequence)

        sequence_len = len(sequence)
        if self.read_length > sequence_len:
            raise ValueError('Read length exceeds sequence length')

        num_reads = self.calculate_num_reads(sequence_len)

        all_reads = list(range(num_reads))
        batches = [
            all_reads[i:i + self.batch_size]
            for i in range(0, num_reads, self.batch_size)
        ]

        worker_args = [
            (
                batch,
                sequence_header,
                sequence,
                self.read_length,
                self.coverage,
                self.get_platform_name(),
                self.profile,
                self.sse_trigger_map,
                self.homopolymer_pos,
                self.reverse_complement,
            )
            for batch in batches
        ]

        reads = []
        with ProcessPoolExecutor(max_workers=self.num_workers) as pool:
            futures = [pool.submit(_simulate_batch, *args) for args in worker_args]
            for future in as_completed(futures):
                reads.extend(future.result())

        reads.sort(key=lambda r: int(r["id"].rsplit("_", 1)[-1]))
        return reads

    def find_homopolymer_runs(self, sequence: str, min_length: int = 3) -> Dict[int, int]:
        """
            Identify homopolymer runs in the sequence.

            :param:
                sequence:   Input genome sequence
                min_length: Minimum run length to record (default 3)

            :return:
                Dict - start position;run length
                e.g. {12: 6, 40: 5} means a run of 6 at pos 12, run of 5 at pos 40
            """
        homopolymer_pos = {}
        i = 0
        while i < len(sequence):
            base = sequence[i]
            run_start = i
            while i < len(sequence) and sequence[i] == base:
                i += 1
            run_length = i - run_start
            if run_length >= min_length:
                for p in range(run_start, run_start + run_length):
                    homopolymer_pos[p] = run_length

        return homopolymer_pos

    def load_quality_profile(self, filename: str) -> dict[int, dict[str, dict]]:
        data = defaultdict(lambda: defaultdict(dict))
        with open(filename, "r") as f:
            lines = f.readlines()

        for i in range(0, len(lines), 2):
            if i+1 >= len(lines):
                break

            score_line = lines[i].strip().split("\t")
            freq_line = lines[i+1].strip().split("\t")

            if len(score_line) < 2 or len(freq_line) < 2:
                continue

            try:
                base = score_line[0]
                position = int(score_line[1])
                scores = list(map(int, score_line[2:]))
                frequencies = list(map(int, freq_line[2:]))
            except (ValueError):
                raise ValueError(f"Improper profile entry at line {i + 1}")

            cumulative = []
            running_total = 0

            for freq in frequencies:
                running_total += freq
                cumulative.append(running_total)

            try:
                total = cumulative[-1]
            except IndexError:
                total = 0

            data[position][base] = {
                "scores": scores,
                "cumulative": cumulative,
                "total": total
            }

        return {pos: dict(bases) for pos, bases in data.items()}

    def generate_quality_scores(self, read_seq, start_pos, end_pos, strand, is_dephased):
        return generate_quality_scores(
            read_seq, start_pos, end_pos, strand, is_dephased,
            self.profile, self.sse_trigger_map, self.homopolymer_pos
        )

    def inject_errors(self, sequence, quality_scores, start_pos, strand, is_dephased):
        return inject_errors(
            sequence, quality_scores, start_pos, strand, is_dephased,
            self.sse_trigger_map, self.homopolymer_pos
        )