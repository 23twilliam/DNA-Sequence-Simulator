import random

# Generate random reads from sequence with simulated sequencing errors and quality scores
def generate_reads(sequence, read_length=150, coverage=30, error_rate=0.001):

    sequence_len = len(sequence)

    num_reads = int((sequence_len * coverage) / read_length)

    reads = []
    bases = ['A', 'T', 'G', 'C']

    for i in range(num_reads):
        if read_length > sequence_len:
            raise ValueError('Read length exceeds sequence length')
        start_pos = random.randint(0, sequence_len - read_length)
        read_seq = sequence[start_pos:start_pos + read_length]

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

        # Creates a dictionary with id, description, sequence, quality scores
        read_record = {
            "id": f"read_{i + 1}",
            "description": f"simulated_read start={start_pos} length={len(read_seq_with_errors)}",
            "sequence": read_seq_with_errors,
            "qualities": qualities
        }

        reads.append(read_record)

    return reads

def write_fastq(reads, filename):
    with open(filename, "w") as f:
        for read in reads:
            f.write(f"@{read['id']} {read['description']}\n")
            f.write(f"{read['sequence']}\n")
            f.write("+\n")
            # we add 33 to phred scores and then convert to a char to get ASCII representation of the Phred score
            # https://help.basespace.illumina.com/files-used-by-basespace/quality-scores
            f.write("".join(chr(q + 33) for q in read["qualities"]) + "\n")

def read_sequences(filename):
    sequences = []
    header = None
    current_sequence = [] # originally used "" however appending to a string is much slower than a list in python
    # This is as += for a string, creates a new strinbg, copies old string into it, appends the new line, gets rid of
    # old string O(n^2) as it is in a loop. Whereas appending to a list is O(1) so O(N) for this loop.

    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                # If we already have a header then we have reached the end of a sequence, which we must add it to the
                # list of sequences
                if header:
                    sequences.append((header, ''.join(current_sequence)))
                current_header = line[1:]  # Remove the '>' character
                # Reset the current sequence to be blank
                header = line
                current_sequence = []
            else:
                current_sequence.append(line)
        if header:
            sequences.append((header, ''.join(current_sequence)))

    return sequences


if __name__ == "__main__":
    example_input = "..\\Example_FASTA\\gene.fna"

    example_sequences = read_sequences(example_input)
    reads = generate_reads(example_sequences[1][1])

    write_fastq(reads, "simulated_reads.fastq")
