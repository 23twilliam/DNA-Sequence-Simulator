import bisect
from ..base_simulator import BaseSimulator
from typing import List, Dict, Tuple
import random
from collections import defaultdict
import time

SUBST_CHOICES: Dict[str, Tuple[List[str], List[float]]] = {
    "A": (["C", "G", "T"], [0.04, 0.07, 0.07]),
    "C": (["A", "G", "T"], [0.13, 0.04, 0.11]),
    "G": (["A", "C", "T"], [0.07, 0.04, 0.11]),
    "T": (["A", "C", "G"], [0.07, 0.04, 0.07]),
}

PREPHASE_PROB = {1: 0.64, 2: 0.71, 3: 0.79, 4: 0.86}
PREPHASE_MAX = 0.86
LAGGING_PROB = 0.56
DEPHASING_RATE = 0.0515
HOMOPOLYMER_DECAY = 0.5  # TODO Double check research for better decay rate
SSE_LOW, SSE_HIGH = 3, 9


class IlluminaSimulator(BaseSimulator):
    def __init__(self, read_length: int, coverage: float, profile_filename: str, error_rate: float = 0.001, sequence: str = None):
        """
        Initialise simulator with the Illumina profile.

        :param:
            read_length: Length of reads to generate
            coverage: Desired coverage depth
            profile_filename: file containing the quality profile
            error_rate: Base error rate (may be overridden by platform)
        """
        super().__init__(read_length, coverage, error_rate) # initialise BaseSimulator
        self.sequence = sequence
        self.profile = self.load_quality_profile(profile_filename)
        self.gc_content = None
        self.sse_trigger_map = None
        self.current_strand = 'F'
        self.current_is_dephased = False

    def get_platform_name(self) -> str:
        return "Illumina Simulator"

    # https://stackoverflow.com/questions/4664850/how-to-find-all-occurrences-of-a-substring

    def build_sse_trigger_map(self, sequence: str, impacted_area: int = 40) -> Dict[int, str]:
        trigger_map = {}
        seq_len = len(sequence)

        def mark_area(start: int, strand: str):
            if strand == 'F':
                for j in range(start + 1, min(start + impacted_area + 1, seq_len)):
                    trigger_map[j] = 'F'
            else:
                for j in range(max(0, start - impacted_area), start):
                    trigger_map[j] = 'R'

        i = sequence.find("GGC")
        while i != -1:
            mark_area(i, 'F')
            i = sequence.find("GGC", i + 1)

        # GCC = GGC on reverse strand
        i = sequence.find("GCC")
        while i != -1:
            mark_area(i + 2, 'R')
            i = sequence.find("GCC", i + 1)

        """
                Identify positions in the sequence where inverted repeats exist within a sliding window. Marks the impact area downstream of each trigger.
                Criteria for a valid trigger:
                - Repeat must be longer than 8bp (k=9 minimum)
                - Reverse complement must appear within +-40bp, instead of looking forward and backward, we just look forward
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
                    mark_area(i, 'F')
                    mark_area(j, 'R')
                    break

        return trigger_map

    def generate_reads(self, sequence_header: str, sequence: str):
        self.gc_content = self.calculate_gc_content(sequence)
        start = time.time()
        self.sse_trigger_map = self.build_sse_trigger_map(sequence)
        print(f"SSE map built in {time.time() - start:.2f}s")
        self.homopolymer_pos = self.find_homopolymer_runs(sequence)

        # loop here is identical to BaseSimulator.generate_reads except
        # for the three lines that set per-read state
        sequence_len = len(sequence)
        if self.read_length > sequence_len:
            raise ValueError('Read length exceeds sequence length')

        num_reads = self.calculate_num_reads(sequence_len)
        reads = []

        for i in range(num_reads):
            start_pos = self.sample_position(sequence_len)
            end_pos = start_pos + self.read_length

            self.current_strand = random.choice(['F', 'R'])
            self.current_is_dephased = random.random() < DEPHASING_RATE

            read_seq = sequence[start_pos:end_pos]
            if self.current_strand == 'R':
                read_seq = self.reverse_complement(read_seq)

            quality_scores = self.generate_quality_scores(read_seq, start_pos)
            read_seq_with_errors, final_qualities, errors_injected = self.inject_errors(
                read_seq, quality_scores, start_pos
            )

            reads.append({
                "id": f"{sequence_header}_read_{i + 1}",
                "description": (
                    f"simulated_read platform={self.get_platform_name()} "
                    f"start={start_pos + 1} end={end_pos} "
                    f"length={len(read_seq_with_errors)} coverage={self.coverage}"
                ),
                "sequence": read_seq_with_errors,
                "qualities": final_qualities,
                "true_start": start_pos + 1,
                "true_end": end_pos,
                "errors": errors_injected,
            })

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

    def weighted_substitution(self, original: str) -> str:
        keys, weights = SUBST_CHOICES[original]
        return random.choices(keys, weights=weights, k=1)[0]

    def load_quality_profile(self, filename: str) -> defaultdict[int, defaultdict[str, dict]]:
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
            except (ValueError, IndexError):
                continue

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

        return data

    def generate_quality_scores(self, read_seq: str, start_pos: int) -> List[int]:
        qualities = []

        for idx, base in enumerate(read_seq):
            genome_pos = start_pos + idx
            scores = self.profile[idx][base]["scores"]
            cumulative = self.profile[idx][base]["cumulative"]
            total = self.profile[idx][base]["total"]

            r = random.uniform(0, total)
            score_idx = min(bisect.bisect_left(cumulative, r), len(scores) - 1)
            q = scores[score_idx]

            # SSE Penalty
            trigger_strand = self.sse_trigger_map.get(genome_pos)
            if self.current_is_dephased and trigger_strand == self.current_strand:
                q = q - random.randint(SSE_LOW, SSE_HIGH)

            # Homopolymer
            run_length = self.homopolymer_pos.get(genome_pos, 0)
            if run_length > 1:
                q = q - int((run_length - 1) * HOMOPOLYMER_DECAY)

            qualities.append(max(2, q))

        return qualities

    def inject_errors(self, sequence: str, quality_scores: List[int], start_pos: int) -> Tuple[str, List[int], List[Dict]]:
        read_list = list(sequence)
        errors_injected = []

        for i in range(len(read_list)):
            error_probability = self.quality_to_error_probability(quality_scores[i])
            if random.random() >= error_probability:
                continue

            original_base = read_list[i]
            mutated_base = None
            mode = "background"

            if self.current_is_dephased:

                # Mode A — Lagging
                if self.sse_trigger_map.get(start_pos + i) == self.current_strand:
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

                # Mode B — Pre-phasing
                if mutated_base is None:
                    run_len = self.homopolymer_pos.get(start_pos + i, 0)
                    if run_len > 0 and i + 1 < len(sequence):
                        p = PREPHASE_PROB.get(run_len, PREPHASE_MAX)
                        next_base = sequence[i + 1]
                        if next_base != original_base and random.random() < p:
                            mutated_base = next_base
                            mode = "prephasing"

            # Mode C — Background
            if mutated_base is None:
                mutated_base = self.weighted_substitution(original_base)
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

        if self.current_is_dephased:
            prephase_count = sum(1 for e in errors_injected if e["mode"] == "prephasing")
            if prephase_count >= 3:
                trim = random.randint(1, 3)
                read_seq_with_errors = read_seq_with_errors[:-trim]
                quality_out = quality_out[:-trim]

        return read_seq_with_errors, quality_out, errors_injected