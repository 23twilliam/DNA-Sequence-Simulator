import bisect
from ..base_simulator import BaseSimulator
from typing import List, Dict, Tuple
import random
from collections import defaultdict

class IlluminaSimulator(BaseSimulator):
    def __init__(self, read_length: int, coverage: float, profile_filename: str, error_rate: float = 0.001):
        """
        Initialise simulator with the Illumina profile.

        Args:
            read_length: Length of reads to generate
            coverage: Desired coverage depth
            profile_filename: file containing the quality profile
            error_rate: Base error rate (may be overridden by platform)
        """
        super().__init__(read_length, coverage, error_rate) # initialise BaseSimulator
        self.profile = self.load_quality_profile(profile_filename)

    def get_platform_name(self) -> str:
        return "Illumina Simulator"

    # https://stackoverflow.com/questions/4664850/how-to-find-all-occurrences-of-a-substring

    def find_all_ggc(self, sequence: str) -> List[int]:
        """
        Returns all the position of the patter "GGC" in the sequence.
        :param sequence:
        :return: list of indexes of the sequence where GGC can be found
        """
        i = sequence.find("GGC")
        positions = []
        while i != -1:
            positions.append(i)
            i = sequence.find("GGC", i+1)
        return positions

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