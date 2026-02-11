def read_sequences(filename):
    sequences = []
    header = None
    current_sequence = [] # originally used "" however appending to a string is much slower than a list in python
    # This is as += for a string, creates a new string, copies old string into it, appends the new line, gets rid of
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