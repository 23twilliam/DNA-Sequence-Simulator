def write_fastq(reads, filename):
    with open(filename, "w") as f:
        for read in reads:
            f.write(f"@{read['id']} {read['description']}\n")
            f.write(f"{read['sequence']}\n")
            f.write("+\n")
            # we add 33 to phred scores and then convert to a char to get ASCII representation of the Phred score
            # https://help.basespace.illumina.com/files-used-by-basespace/quality-scores
            f.write("".join(chr(q + 33) for q in read["qualities"]) + "\n")