import bisect
from ..base_simulator import BaseSimulator
from typing import List, Dict, Tuple
import random
from collections import defaultdict

class IlluminaSimulator(BaseSimulator):
    def __init__(self, read_length: int, coverage: float, profile_filename: str, error_rate: float = 0.001, sequence: str = None):
        """
        Initialise simulator with the Illumina profile.

        Args:
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

    def get_platform_name(self) -> str:
        return "Illumina Simulator"

    # https://stackoverflow.com/questions/4664850/how-to-find-all-occurrences-of-a-substring

    def build_sse_trigger_map(self, sequence: str, impacted_area: int = 40) -> set:
        ggc_positions = set()
        for pos in self.find_all_ggc(sequence):
            for i in range(pos, pos + impacted_area):
                ggc_positions.add(i)

        ir_positions = self.find_inverted_repeats(sequence)
        return ggc_positions | ir_positions

    def find_inverted_repeats(self, sequence: str, k: int = 9, window: int = 40, impacted_area: int = 40) -> set:
        """
        Identify positions in the sequence where inverted repeats exist within a sliding window. Marks the impact area downstream of each trigger.
        Criteria for a valid trigger:
        - Repeat must be longer than 8bp (k=9 minimum)
        - Reverse complement must appear within +-40bp, instead of looking forward and backward, we just look forward
         from the current position (40bp total window)
        - Errors are in the 40 bases downstream of the trigger position
        :param:
            sequence: Input genome sequence
            k: K-mer length (must be > 8bp per trigger)
            window: Maximum distance to search for reverse complement
            impacted_area: Number of bases downstream to mark

        :return:Set of sequence positions flagged as active
        """
        trigger_positions = set()

        for i in range(len(sequence) - k):
            kmer = sequence[i:i + k]
            reverse_comp = self.reverse_complement(kmer)

            # Search only forward within the window the reverse complement must appear downstream of the sequence
            search_window = sequence[i + k: i + k + window]

            if reverse_comp in search_window:
                for pos in range(i, i + impacted_area):
                    trigger_positions.add(pos)

        return trigger_positions

    def find_all_ggc(self, sequence: str) -> List[int]:
        """
        Returns all the position of the patter "GGC" in the sequence.
        :param: sequence: the original genome sequence
        :return: list of indexes of the sequence where GGC can be found
        """
        i = sequence.find("GGC")
        positions = []
        while i != -1:
            positions.append(i)
            i = sequence.find("GGC", i+1)
        return positions

    def generate_reads(self, sequence_header: str, sequence: str):
        self.gc_content = self.calculate_gc_content(sequence)
        self.sse_trigger_map = self.build_sse_trigger_map(sequence)
        return super().generate_reads(sequence_header, sequence)

    # TODO Implement
    def find_homopolymer_runs(self, sequence: str) -> List[int]:
        pass

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


    def generate_quality_scores(self, read_seq: str) -> List[int]:
        # would have like to use a fenwick tree here, however as we never update the frequencies it doesnt make sense,
        # a binary search is more effective
        qualities = []
        # TODO fix length issue
        # this will break if the profile file has an amount of positions for base that is lower than the read length
        for idx, base in enumerate(read_seq):
            scores = self.profile[idx][base]["scores"]
            cumulative = self.profile[idx][base]["cumulative"]
            total = self.profile[idx][base]["total"]

            r = random.uniform(0, total)
            score_idx = bisect.bisect_left(cumulative, r)
            qualities.append(scores[score_idx])

        return qualities

    def inject_errors(self, sequence: str, quality_scores: List[int]) -> Tuple[str, List[int], List[Dict]]:
        """
        Inject substitution errors based on quality scores.
        Lower quality -> higher error probability.
        """
        read_list = list(sequence)
        errors_injected = []

        for i in range(len(read_list)):
            # Use quality score to determine error probability
            error_probability = self.quality_to_error_probability(quality_scores[i])

            if random.random() < error_probability:
                original_base = read_list[i]
                mutated_base = random.choice([b for b in self.bases if b != original_base])
                read_list[i] = mutated_base

                errors_injected.append({
                    'position': i,
                    'type': 'substitution',
                    'original': original_base,
                    'mutated': mutated_base,
                    'quality': quality_scores[i]  # Include quality in ground truth
                })

        read_seq_with_errors = ''.join(read_list)

        # For substitutions only, quality list length doesn't change
        # But for indels (in other platforms), we need to adjust
        return read_seq_with_errors, quality_scores, errors_injected