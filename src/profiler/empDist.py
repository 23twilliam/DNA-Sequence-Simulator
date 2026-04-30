import sys
import os

if len(sys.argv) != 3:
    print("Usage: python empDist.py freqFile outputFile")
    sys.exit(1)

infile = sys.argv[1]
outfile = sys.argv[2]

NUM_QUAL = 71
nucleotides = [".", "A", "T", "G", "C", "N"]

# Create storage: data[nucleotide][quality][position]
data = {nuc: [] for nuc in nucleotides}

with open(infile) as f:
    for count, line in enumerate(f):
        if 2<= count <= NUM_QUAL + 1:
            values = list(map(int, line.strip().split("\t"))) # strip line into its elements eg 1   1 becomes 1,1. then
            # map them to an int. then put that into a list.

            read_length = len(values) // 5

            # Slice into nucleotide blocks
            blocks = {
                "A": values[0:read_length],
                "T": values[read_length:2*read_length],
                "G": values[2*read_length:3*read_length],
                "C": values[3*read_length:4*read_length],
                "N": values[4*read_length:5*read_length],
            }

            # Store per-quality data
            for nuc in ["A", "T", "G", "C", "N"]:
                data[nuc].append(blocks[nuc])

            # Build ALL by summing nucleotides
            # zip should be fine as all are the same length (should all be read_length len). but check anyway just in case
            if len(blocks["A"]) != len(blocks["T"]) != len(blocks["G"]) != len(blocks["C"]) != len(blocks["N"]):
                print("Nucleotide blocks have mismatching lengths")
                sys.exit(1)

            all_block = [
                sum(vals) for vals in zip(
                    blocks["A"],
                    blocks["T"],
                    blocks["G"],
                    blocks["C"],
                    blocks["N"],
                )
            ]
            data["."].append(all_block)

read_length = len(data["A"][0])

out_path = os.path.join("profiles", os.path.basename(outfile))
os.makedirs("profiles", exist_ok=True)
with open(out_path, "a") as out:
    # loop over the possible nucleotides
    for nuc in nucleotides:
        # loop over the possible positions
        for pos in range(read_length):

            # Get non-zero quality scores
            non_zero_quals = [
                qual for qual in range(NUM_QUAL)
                if data[nuc][qual][pos] != 0
            ]

            # Print quality scores
            out.write(f"{nuc}\t{pos}\t")
            out.write("\t".join(map(str, non_zero_quals)) + "\n")

            # display the cumulative frequency for each quality score
            cumulative = 0
            out.write(f"{nuc}\t{pos}\t")
            for qual in non_zero_quals:
                cumulative += data[nuc][qual][pos]
                out.write(str(cumulative) + "\t")

            out.write("\n")