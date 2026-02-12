import random
import io_functions

# Generate random reads from sequence with simulated sequencing errors and quality scores
def generate_reads(sequence_header, sequence, read_length, coverage, error_rate=0.001):

    sequence_len = len(sequence)

    if read_length == 0:
        raise ValueError('Read length must be greater than 0')
    if read_length > sequence_len:
        raise ValueError('Read length exceeds sequence length')

    num_reads = int((sequence_len * coverage) / read_length)

    reads = []
    bases = ['A', 'T', 'G', 'C']

    for i in range(num_reads):
        start_pos = random.randint(0, sequence_len - read_length)
        read_seq = sequence[start_pos:start_pos + read_length]

        # Introducing sequencing errors
        read_list = list(read_seq)
        for j in range(len(read_list)):
            if random.random() < error_rate:
                read_list[j] = random.choice([b for b in bases if b != read_list[j]])

        read_seq_with_errors = ''.join(read_list)


        # Find way to do this more effectively using list comprehension
        # https://www.geeksforgeeks.org/python/generating-random-number-list-in-python/
        # Assign quality scores 10% at Q20-Q30, 90% at Q30-Q40
        qualities = [
            random.randint(20, 30) if random.random() < 0.1 else random.randint(30,40)
            for _ in range(len(read_seq_with_errors))
        ]

        # Creates a dictionary with id, description, sequence, quality scores
        read_record = {
             "id": f"{sequence_header}_read_{i + 1}",
            "description": f"simulated_read start={start_pos} length={len(read_seq_with_errors)} coverage={coverage}",
            "sequence": read_seq_with_errors,
            "qualities": qualities
        }

        reads.append(read_record)

    return reads

def run_example(read_length, coverage):
    example_input = "../Example_FASTA/gene.fna"

    example_sequences = io_functions.read_sequences(example_input)
    reads = []
    for sequence_header, sequence in example_sequences.items():
        reads.extend(generate_reads(sequence_header, sequence, read_length, coverage))

    io_functions.write_fastq(reads, "simulated_reads.fastq")