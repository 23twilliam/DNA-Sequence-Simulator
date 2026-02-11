import random
import io_functions

# Generate random reads from sequence with simulated sequencing errors and quality scores
def generate_reads(sequence, read_length, coverage, error_rate=0.001):

    sequence_len = len(sequence)

    num_reads = int((sequence_len * coverage) / read_length)

    reads = []
    bases = ['A', 'T', 'G', 'C']

    for i in range(num_reads):
        if read_length > sequence_len:
            raise ValueError('Read length exceeds sequence length')
        start_pos = random.randint(0, sequence_len - read_length)
        read_seq = sequence[start_pos:start_pos + read_length]

        # Introducing sequencing errors
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

        # Creates a dictionary with id, description, sequence, quality scores
        read_record = {
            "id": f"read_{i + 1}",
            "description": f"simulated_read start={start_pos} length={len(read_seq_with_errors)} coverage={coverage}",
            "sequence": read_seq_with_errors,
            "qualities": qualities
        }

        reads.append(read_record)

    return reads

def run_example(read_length, coverage):
    example_input = "../Example_FASTA/gene.fna"

    example_sequences = io_functions.read_sequences(example_input)
    reads = generate_reads(example_sequences[1][1], read_length, coverage)

    io_functions.write_fastq(reads, "simulated_reads.fastq")