import math
import os
from abc import ABC, abstractmethod
from typing import List, Dict, Tuple
from collections import defaultdict
import random
from concurrent.futures import ProcessPoolExecutor, as_completed

class BaseSimulator(ABC):
    """
    Abstract base class for all sequencing platform simulators.
    """

    def __init__(self, read_length: int, coverage: float, profile_path: str, num_workers: int = os.cpu_count(), batch_size: int = 5000):
        """
        Initialise simulator with basic parameters.

        Args:
            read_length: Length of reads to generate
            coverage: Desired coverage depth
            profile_filename: File containing the quality profile
            num_workers:      Number of worker processes (default: os.cpu_count())
            batch_size:       Reads per worker batch (default: 5000)
        """
        self.read_length = read_length
        self.coverage = coverage
        self.bases = ['A', 'T', 'G', 'C']
        self.sequence = None
        self.num_workers = num_workers
        self.batch_size = batch_size
        self.profile = self.load_quality_profile(profile_path)

    @abstractmethod
    def _precompute(self, sequence):
        """Build platform-specific precomputed structures (e.g. homopolymer maps or SSE triggers) before batch dispatch."""
        pass

    @abstractmethod
    def _get_worker_fn(self):
        """Return the module level batch worker function for this simulator."""
        pass

    def _get_worker_args(self, batch, sequence_header, sequence) -> tuple:
        """Return the args tuple to pass the worker for this batch."""
        pass

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
        return math.ceil((sequence_length * self.coverage) / self.read_length)

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
                    homopolymer_pos[p] = run_length - (p - run_start)

        return homopolymer_pos

    def load_quality_profile(self, filename: str) -> Dict:
        """
        Load empirical quality profile from file
        Profile stores cumulative frequencies
        :param filename: name of profile file
        :return: dict with the structure of [pos][bas] = {scores, cumulative, total]
        """
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
                cumulative = list(map(int, freq_line[2:]))
            except (ValueError):
                raise ValueError(f"Improper profile entry at line {i + 1}")

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

    # def generate_reads(self, sequence_header: str, sequence: str) -> List[Dict]:
    #     """
    #     Main read generation method.
    #     Quality scores -> Errors -> Output
    #     """
    #     self.sequence = sequence
    #     sequence_len = len(sequence)
    #
    #     if self.read_length > sequence_len:
    #         raise ValueError('Read length exceeds sequence length')
    #
    #     num_reads = self.calculate_num_reads(sequence_len)
    #     reads = []
    #
    #     for i in range(num_reads):
    #         start_pos = self.sample_position(sequence_len)
    #         end_pos = start_pos + self.read_length
    #         read_seq = sequence[start_pos:end_pos]
    #
    #         strand = random.choice(['F', 'R'])
    #         if strand == 'R':
    #             read_seq = self.reverse_complement(read_seq)
    #
    #         # STEP 1: Generate quality scores
    #         quality_scores = self.generate_quality_scores(read_seq)
    #
    #         # STEP 2: Use quality scores to guide error injection
    #         read_seq_with_errors, final_qualities, errors_injected = self.inject_errors(
    #             read_seq,
    #             quality_scores
    #         )
    #
    #         # Create read record
    #         read_record = {
    #             "id": f"{sequence_header}_read_{i + 1}",
    #             "description": (
    #                 f"simulated_read platform={self.get_platform_name()} "
    #                 f"start={start_pos + 1} end={end_pos} "
    #                 f"length={len(read_seq_with_errors)} coverage={self.coverage}"
    #             ),
    #             "sequence": read_seq_with_errors,
    #             "qualities": final_qualities,
    #             "true_start": start_pos + 1,
    #             "true_end": end_pos,
    #             "errors": errors_injected
    #         }
    #
    #         reads.append(read_record)
    #
    #     return reads

    def generate_reads(self, sequence_header: str, sequence: str) -> List[Dict]:
        sequence_len = len(sequence)

        if self.read_length > sequence_len:
            raise ValueError('Read length exceeds sequence length')

        # Platform specific precomputation
        self._precompute(sequence)

        num_reads = self.calculate_num_reads(sequence_len)
        all_reads = list(range(num_reads))
        batches = [
            all_reads[i:i + self.batch_size]
            for i in range(0, num_reads, self.batch_size)
        ]

        worker_fn = self._get_worker_fn()
        worker_args = [
            self._get_worker_args(batch, sequence_header, sequence)
            for batch in batches
        ]

        reads = []
        with ProcessPoolExecutor(max_workers=self.num_workers) as pool:
            futures = [pool.submit(worker_fn, *args) for args in worker_args]
            for future in as_completed(futures):
                reads.extend(future.result())

        reads.sort(key=lambda r: int(r["id"].rsplit("_", 1)[-1]))
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