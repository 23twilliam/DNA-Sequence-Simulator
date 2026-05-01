import sys
import gzip
from collections import defaultdict

# Making sure that the correct number of arguments are given by the user
if len(sys.argv) != 2:
    print("Usage: python fastqReadAvg.py FastqFile")
    sys.exit(1)

# Getting the filename to the file of interest
infile = sys.argv[1]

# Open normally or via gzip depending on extension
if infile.endswith(".gz"):
    input_file = gzip.open(infile, "rt")
else:
    input_file = open(infile, "r")

# Output file (same behaviour as Perl: basename + .txt)
out_file = infile + ".txt"
out = open(out_file, "w")

# Configuration constants
NUM_SCORES = 71        # possible quality scores
BASES = ["A", "T", "G", "C", "N"]
ASCII_OFFSET = 33      # Phred+33

# Data storage
# data[quality][position][base] = count
data = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))

max_length = 0

# Main FASTQ parsing loop

# FASTQ structure:
# line 0: @header
# line 1: bases
# line 2: +
# line 3: quality scores

for i, line in enumerate(input_file):
    line = line.rstrip()

    # Sequence line
    if i % 4 == 1:
        bases = list(line)
        max_length = max(max_length, len(bases))

    # Quality line
    elif i % 4 == 3:
        scores = list(line)
        if len(bases) != len(scores):
            print("Error in zipping bases and scores, their lengths do not match")
            sys.exit(1)
        for pos, (base, score_char) in enumerate(zip(bases, scores)):
            qual = ord(score_char) - ASCII_OFFSET # ord converts the ASCII character into their number ie ! becomes 33
            if base not in BASES:
                base = "N"

            data[qual][pos][base] += 1

input_file.close() # have to close because we couldnt use with earlier


# Write the longest read length
out.write(f"Longest Length Read: {max_length}\n")

# Column headings
for base in BASES:
    for pos in range(max_length):
        out.write(f"{base}{pos}\t")
out.write("\n")

# Build full matrix and averages

# final_matrix[quality][column]
final_matrix = []
avg_numerators = []
avg_denominators = []

total_columns = max_length * len(BASES)

for qual in range(NUM_SCORES):

    row = []
    for base in BASES:
        for pos in range(max_length):
            value = data[qual][pos][base]
            row.append(value)

    final_matrix.append(row)

# Initialise average tracking
avg_numerators = [0] * total_columns
avg_denominators = [0] * total_columns

# Calculate averages
for qual in range(NUM_SCORES):
    for col in range(total_columns):
        count = final_matrix[qual][col]
        avg_numerators[col] += qual * count
        avg_denominators[col] += count


# Write data matrix + averages
for qual in range(NUM_SCORES + 1):

    if qual == NUM_SCORES:
        # Write average row
        for col in range(total_columns):
            num = avg_numerators[col]
            denom = avg_denominators[col]

            if num != 0 and denom != 0:
                out.write(f"{int(num / denom)}\t")
            else:
                out.write("0\t")
    else:
        for value in final_matrix[qual]:
            out.write(f"{value}\t")

    out.write("\n")

# Write numerator and denominator rows
for row in (avg_numerators, avg_denominators):
    for value in row:
        out.write(f"{value}\t")
    out.write("\n")

out.close()