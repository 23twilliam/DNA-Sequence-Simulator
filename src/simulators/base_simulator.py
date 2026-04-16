from abc import ABC, abstractmethod
from typing import List, Dict, Tuple
import random


class BaseSimulator(ABC):
    """
    Abstract base class for all sequencing platform simulators.
    """

    def __init__(self, read_length: int, coverage: float):
        """
        Initialise simulator with basic parameters.

        Args:
            read_length: Length of reads to generate
            coverage: Desired coverage depth
        """
        self.read_length = read_length
        self.coverage = coverage
        self.bases = ['A', 'T', 'G', 'C']
        self.sequence = None

    @abstractmethod
    def get_platform_name(self) -> str:
        """Return platform name (e.g. 'Illumina', 'PacBio', 'Nanopore')"""
        pass

    @abstractmethod
    def generate_quality_scores(self, read: str) -> List[int]:
        """
        Generate quality scores for a read.
        Quality scores should be generated first because they will influence error injection probability.

        Returns:
            List of quality scores (Phred scale)
        """
        pass

    @abstractmethod
    def inject_errors(self, sequence: str, quality_scores: List[int]) -> Tuple[str, List[int], List[Dict]]:
        """
        Inject sequencing errors into read based on quality scores.
        Errors should be informed by quality scores
        - Lower quality positions -> higher error probability
        - Insertions/deletions will change sequence length -> quality list must adjust

        Args:
            sequence: Original read sequence
            quality_scores: Quality scores that guide error probability

        Returns:
            Tuple of:
            - sequence_with_errors: Modified sequence
            - adjusted_qualities: Quality scores
            - errors_injected: List of errors for ground truth
        """
        pass

    def calculate_num_reads(self, sequence_length: int) -> int:
        """
        Calculate number of reads needed for desired coverage.
        """
        return int((sequence_length * self.coverage) / self.read_length)

    def calculate_gc_content(self, sequence: str) -> float:
        """
        Calculate the relative gc content of the sequence.
        """
        gc_count = sequence.count("G") + sequence.count("C")
        gc_proportion = gc_count / len(sequence)

        return gc_proportion

    def reverse_complement(self, sequence: str) -> str:
        """
        Return the reverse complement of sequence.
        """
        # TODO bitwise implementation?
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
        return ''.join(complement[base] for base in reversed(sequence))

    def sample_position(self, sequence_length: int) -> int:
        """
        Sample a random starting position in the genome.
        """
        return random.randint(0, sequence_length - self.read_length)

    def quality_to_error_probability(self, quality_score: int) -> float:
        """
        Convert Phred quality score to error probability.
        Q = -10 * log10(P)  =>  P = 10^(-Q/10)
        """
        return 10 ** (-quality_score / 10)

    def generate_reads(self, sequence_header: str, sequence: str) -> List[Dict]:
        """
        Main read generation method.
        Quality scores -> Errors -> Output
        """
        self.sequence = sequence
        sequence_len = len(sequence)

        if self.read_length > sequence_len:
            raise ValueError('Read length exceeds sequence length')

        num_reads = self.calculate_num_reads(sequence_len)
        reads = []

        for i in range(num_reads):
            start_pos = self.sample_position(sequence_len)
            end_pos = start_pos + self.read_length

            # Extract original sequence
            read_seq = sequence[start_pos:end_pos]

            # STEP 1: Generate quality scores
            quality_scores = self.generate_quality_scores(read_seq)

            # STEP 2: Use quality scores to guide error injection
            read_seq_with_errors, final_qualities, errors_injected = self.inject_errors(
                read_seq,
                quality_scores
            )

            # Create read record
            read_record = {
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
                "errors": errors_injected
            }

            reads.append(read_record)

        return reads

    def simulate_from_file(self, fasta_file: str, io_functions) -> List[Dict]:
        """
        Generate reads from all sequences in a FASTA file.

        Args:
            fasta_file: Path to input FASTA file
            io_functions: io_functions module
        Returns:
            List of all generated reads
        """
        sequences = io_functions.read_sequences(fasta_file)
        all_reads = []

        for sequence_header, sequence in sequences.items():
            reads = self.generate_reads(sequence_header, sequence)
            all_reads.extend(reads)

        return all_reads