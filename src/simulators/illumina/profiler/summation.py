import sys
from collections import defaultdict


# Argument handling
if len(sys.argv) < 4:
    print("Usage: python summation.py fastqReadAvg_file1 fastqReadAvg_file2 [fastqReadAvg_fileN] outputFile")
    sys.exit(1)

input_files = sys.argv[1:-1]
outfile = sys.argv[-1]

NUM_QUAL = 71
BASES = ["A", "T", "G", "C", "N"]
NUM_BASES = 5

# combined[quality][base][position] = frequency
# https://stackoverflow.com/questions/8419401/python-defaultdict-and-lambda
# https://www.devcuriosity.com/manual/details/python-default-dict/
# Use nested defaultdict to avoid key error when trying to access missing keys automatically creates dicts and
# initialises counts to 0. We also dont need to check the existence of key before adding :
# combined[qual][base][pos] += value
combined = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))

max_length = 0

# Process each input file
for infile in input_files:

    with open(infile) as f:
        lines = f.readlines()

    # Quality score rows are lines 2–72
    quality_lines = lines[2:2 + NUM_QUAL]

    for qual, line in enumerate(quality_lines):

        values = list(map(int, line.strip().split("\t")))
        current_length = len(values) // NUM_BASES

        max_length = max(max_length, current_length)

        for b_index, base in enumerate(BASES):
            start = b_index * current_length
            end = start + current_length

            for pos, value in enumerate(values[start:end]):
                # when we are doing assignment like this we would have to check things in a regular dict. we use default
                # dict so that we do not need to check and create each level
                combined[qual][base][pos] += value

# Write output
with open(outfile, "a") as out:

    # Write header
    for base in BASES:
        for pos in range(max_length):
            out.write(f"{base}{pos}\t")
    out.write("\n")

    # Write combined frequencies
    for qual in range(NUM_QUAL):
        for base in BASES:
            for pos in range(max_length):
                # defaultdict returns 0 for missing keys (if there were any)
                out.write(f"{combined[qual][base][pos]}\t")
        out.write("\n")