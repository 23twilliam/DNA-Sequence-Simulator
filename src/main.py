import argparse
import io_functions
from simulators import IlluminaSimulator, PacBioSimulator
import time

# Default profile paths per platform
DEFAULT_PROFILES = {
    "illumina": "simulators/illumina/profiles/ecoli_hiseqR1.txt",
    #"pacbio": "simulators/illumina/profiles/ecoli_hiseqR1.txt"
}

# Default read lengths per platform
DEFAULT_READ_LENGTHS = {
    "illumina": 150,
    "pacbio": 15000,
    "nanopore": 20000,
}
# PacBio length distribution coefficient of variation
# TODO: check literature
PACBIO_DEFAULT_CV = 0.3

def build_simulator(platform, read_length, coverage, profile_path):
    """Construct the appropriate simulator for the chosen platform."""
    if platform == "illumina":
        return IlluminaSimulator(
            read_length=read_length,
            coverage=coverage,
            profile_filename=profile_path
        )
    else:
        raise ValueError(f"Platform '{platform}' not yet implemented")

def run_simulation(platform, input_file, read_length, coverage, output_format, output_path, profile_path):
    simulator = build_simulator(platform, read_length, coverage, profile_path)

    start = time.time()
    reads = simulator.simulate_from_file(input_file, io_functions)
    elapsed = time.time() - start

    if output_format == "fastq":
        io_functions.write_fastq(reads, output_path)
    elif output_format == "fasta":
        io_functions.write_fasta(reads, output_path)
    else:
        print(f"Output format '{output_format}' not yet supported")
        return

    print(f"Generated {len(reads)} reads using {simulator.get_platform_name()} simulator")
    print(f"Read length: {read_length}, Coverage: {coverage}x")
    print(f"Simulation took {elapsed:.2f}s ({elapsed / len(reads) * 1000:.3f}ms per read)")

def main():
    # for all commands marked as required=False they must be changed to True
    parser = argparse.ArgumentParser(description='DNA Read Generator')

    parser.add_argument('command', choices=['run'],
                        help='Command to execute program')

    parser.add_argument('-i', '--input', required=False,
                        default="../Example_FASTA/GCF_000005845.2_ASM584v2_genomic.fna",
                        help='Input genome file (FASTA format)')

    parser.add_argument('-p', '--platform',
                        choices=['illumina', 'pacbio', 'nanopore'],
                        default='illumina',
                        help='Sequencing platform (default: illumina)')

    parser.add_argument('-l', '--read_length', type=int, default=None,
                        help='Defaults vary by platform: 150 for Illumina, 15000 for PacBio, 20000 for Nanopore.')

    parser.add_argument('-c', '--coverage', type=int, default=30,
                        help='Coverage depth (default: 30)')

    parser.add_argument('-f', '--format',
                        choices=['fasta', 'fastq', 'sam', 'bam', 'bed'],
                        default='fastq',
                        help='Output format (default: fastq)')

    parser.add_argument('-o', '--output',
                        default='simulated_reads',
                        required=False,
                        help='Output file or directory path')

    parser.add_argument('--profile',
                        default=None,
                        choices=['art_illumina_profile_hiseq'],
                        help='Illumina read model (default: art_illumina_profile_hiseq)')

    args = parser.parse_args()

    # We should warn that they want to put in a profile or accept the default chosen option
    # if args.platform == 'illumina' and args.model is NONE:
    #     parser.error('-p illumina requires -m option')

    # Apply platform-specific default read length if user did not specify one
    if args.read_length is None:
        args.read_length = DEFAULT_READ_LENGTHS[args.platform]

    # Apply platform-specific default profile if user did not specify one
    profile_path = args.profile if args.profile else DEFAULT_PROFILES.get(args.platform)
    if profile_path is None:
        parser.error(f"No default profile available for {args.platform}; please supply --profile")

    output_path = f"{args.output}.{args.format}"

    if args.command == 'run':
        run_simulation(
            platform=args.platform,
            input_file=args.input,
            read_length=args.read_length,
            coverage=args.coverage,
            output_format=args.format,
            output_path=output_path,
            profile_path=profile_path,
        )
        print("reads generated")


if __name__ == "__main__":
    main()
