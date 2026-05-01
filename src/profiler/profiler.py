import argparse
import os
import subprocess
import sys
from pathlib import Path

# Run fastqReadAvg on each file
def process_fastq(fq):
    output_file = fq.with_suffix(fq.suffix + ".txt")

    if output_file.exists():
        print(f"Skipping {fq.name} (already processed)")
        return output_file

    print(f"Processing {fq.name}")

    result = subprocess.run([sys.executable, "fastqReadAvg.py", str(fq)],capture_output=True)

    if result.returncode != 0:
        print(result.stderr.decode())
        raise RuntimeError(f"Error processing {fq}")

    return output_file

# Helper to run final three scripts
def build_profile(file_list, output_prefix):

    if not file_list:
        return

    freq_file = f"{output_prefix}.freq.txt"

    if len(file_list) == 1:
        os.rename(file_list[0], freq_file)
    else:
        subprocess.run([sys.executable, "summation.py", *map(str, file_list), freq_file],check=True)

    subprocess.run([sys.executable, "combinedAvg.py", freq_file], check=True)
    subprocess.run([sys.executable, "empDist.py", freq_file, f"{output_prefix}.txt"],check=True)

    print(f"Created {output_prefix}.txt")

    os.remove(freq_file)
    for f in file_list:
        if f.exists():
            os.remove(f)

# Argument parsing
parser = argparse.ArgumentParser(
    description="Create an Illumina read quality profile from multiple FASTQ files"
)

parser.add_argument("output_profile_name")
parser.add_argument("input_fastq_dir")
parser.add_argument("fastq_extension")
args = parser.parse_args()

out_prefix = args.output_profile_name
input_dir = Path(args.input_fastq_dir)
ext = args.fastq_extension

if not input_dir.exists():
    sys.exit(f"Error: directory {input_dir} does not exist")

fastq_files = sorted(input_dir.glob(f"*.{ext}"))

if not fastq_files:
    sys.exit(f"Error: no *.{ext} files found in {input_dir}")

print(f"Found {len(fastq_files)} FASTQ files")
output_files = [process_fastq(file) for file in fastq_files]

# Separate single, R1, R2

single = []
r1 = []
r2 = []

for f in output_files:
    name = f.name

    if "_1." in name or ".1." in name:
        r1.append(f)
    elif "_2." in name or ".2." in name:
        r2.append(f)
    else:
        single.append(f)

# Build profiles

build_profile(single, out_prefix)
build_profile(r1, f"{out_prefix}R1")
build_profile(r2, f"{out_prefix}R2")