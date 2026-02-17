from ..base_simulator import BaseSimulator
from typing import List, Dict, Tuple
import random


class RandomSimulator(BaseSimulator):

    def get_platform_name(self) -> str:
        return "Random Simulator"

    def generate_quality_scores(self, read_length: int) -> List[int]:
        """Generate quality scores: 10% at Q20-Q30, 90% at Q30-Q40"""
        qualities = [
            random.randint(20, 30) if random.random() < 0.1 else random.randint(30, 40)
            for _ in range(read_length)
        ]
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