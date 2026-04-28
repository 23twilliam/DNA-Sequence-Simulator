def write_fastq(reads, filename):
    lines = []
    for read in reads:
        lines.append(f"@{read['id']} {read['description']}")
        lines.append(read["sequence"])
        lines.append("+")
        # we add 33 to phred scores and then convert to a char to get ASCII representation of the Phred score
        # https://help.basespace.illumina.com/files-used-by-basespace/quality-scores
        lines.append("".join(chr(q + 33) for q in read["qualities"]))

    with open(filename, "w") as f:
        f.write("\n".join(lines))

def write_fasta(reads, filename):
    lines = []
    for read in reads:
        # no @ on fasta output
        lines.append(f"{read['id']} {read['description']}")
        lines.append(read["sequence"])

    with open(filename, "w") as f:
        f.write("\n".join(lines))