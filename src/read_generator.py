from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import random


# Generate random reads from sequence with simulated sequencing errors and quality scores
def generate_reads(sequence, read_length=150, coverage=30, error_rate=0.001):

    sequence_len = len(sequence)

    num_reads = int((sequence_len * coverage) / read_length)

    reads = []
    bases = ['A', 'T', 'G', 'C']

    for i in range(num_reads):
        start_pos = random.randint(0, sequence_len - read_length)
        read_seq = sequence[start_pos:start_pos + read_length]
        # TODO make this work on cases where read length is > sequence_len

        # Intorducing sequencing errors
        read_list = list(read_seq)
        for j in range(len(read_list)):
            if random.random() < error_rate:
                read_list[j] = random.choice([b for b in bases if b != read_list[j]])

        read_seq_with_errors = ''.join(read_list)

        qualities = []
        # Assign quality scores 10% at Q20-Q30, 90% at Q30-Q40
        for k in range(len(read_seq_with_errors)):
            if random.random() < 0.1:
                quality = random.randint(20, 30)
            else:
                quality = random.randint(30, 40)
            qualities.append(quality)

        # Create SeqRecord with quality scores
        read_record = SeqRecord(
            Seq(read_seq_with_errors),
            id=f"read_{i + 1}",
            description=f"simulated_read start={start_pos} length={len(read_seq_with_errors)}",
            letter_annotations={"phred_quality": qualities}
        )

        reads.append(read_record)

    return reads

def read_sequence(filename):
    # This currently gets just the first sequence
    with open(filename) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            return str(record.seq)


if __name__ == "__main__":
    example_input = "..\\Example_FASTA\\gene.fna"

    example_sequence = read_sequence(example_input)
    reads = generate_reads(example_sequence)

    with open("simulated_reads.fastq", "w") as output_handle:
        SeqIO.write(reads, output_handle, "fastq")
