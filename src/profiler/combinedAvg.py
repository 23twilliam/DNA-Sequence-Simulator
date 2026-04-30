import sys

if len(sys.argv) != 2:
    print("Usage: python combinedAvg.py FreqFile")
    sys.exit(1)

infile = sys.argv[1]

with open(infile, "r") as file:
    lines = file.readlines()

NUM_SCORES = 71
data_lines = lines[2:2 + NUM_SCORES]  # skip header
all_scores = list(range(NUM_SCORES))

avg_score = [[], []]  # [weighted_sums, frequency_sums]

count = 0

for line in data_lines:
    temp = line.strip().split("\t")

    # Convert to integers
    temp = [int(x) for x in temp]

    # Initialise storage on first iteration
    if count == 0:
        avg_score[0] = [0] * len(temp)
        avg_score[1] = [0] * len(temp)

    for pos, value in enumerate(temp):
        avg_score[0][pos] += value * all_scores[count]
        avg_score[1][pos] += value

    count += 1

outfile = infile

with open(outfile, "a") as out:
    # Write out the numerators
    out.write("\t".join(str(x) for x in avg_score[0]) + "\t\n")
    out.write("\n")
    # Write out the denominators
    out.write("\t".join(str(x) for x in avg_score[1]) + "\t\n")
    out.write("\n")

    # Write out the average quality score for each position
    if not len(avg_score[0]) == len(avg_score[1]):
        print(f"Length mismatch: numerators={len(avg_score[0])}, denominators={len(avg_score[1])}")
        sys.exit(1)

    averages = []
    for numerator, denominator in zip(avg_score[0], avg_score[1]):
        if numerator == 0 or denominator == 0:
            averages.append("0")
        else:
            averages.append(str(numerator / denominator))

    out.write("\t".join(averages))