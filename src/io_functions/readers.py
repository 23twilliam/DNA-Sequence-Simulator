def read_sequences(fasta_file,):
    sequences = {}
    header = None
    current_sequence = []  # originally used "" however appending to a string is much slower than a list in python
    # This is as += for a string, creates a new string, copies old string into it, appends the new line, gets rid of
    # old string O(n^2) as it is in a loop. Whereas appending to a list is O(1) so O(N) for this loop.
    lines = open(fasta_file).readlines()
    for line in lines:
        line = line.strip()
        if line.startswith('>'):
            # If we already have a header then we have reached the end of a sequence, which we must add it to the
            # list of sequences
            if header:
                sequences[header] = ''.join(current_sequence)
            # Reset the current sequence to be blank
            header = line
            current_sequence = []
        else:
            current_sequence.append(line)

    sequences[header] = ''.join(current_sequence)
    return sequences