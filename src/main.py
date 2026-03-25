import argparse
import io_functions
from simulators import RandomSimulator, IlluminaSimulator
import time

def run_example(read_length, coverage, output_format):
    """The original run_example() refactored to use the simulator class"""
    example_input = "../Example_FASTA/GCF_000005845.2_ASM584v2_genomic.fna"

    # Create simulator instance
    simulator = IlluminaSimulator(
        read_length=read_length,
        coverage=coverage,
        profile_filename="simulators/illumina/profiles/ecoli_hiseqR1.txt"
    )

    # Generate reads
    start = time.time()
    reads = simulator.simulate_from_file(example_input, io_functions)
    elapsed = time.time() - start

    # Write output
    if output_format == "fastq":
        io_functions.write_fastq(reads, f"simulated_reads.{output_format}")
    elif output_format == "fasta":
        io_functions.write_fasta(reads, f"simulated_reads.{output_format}")
    else:
        print("error")


    print(f"Generated {len(reads)} reads using {simulator.get_platform_name()} simulator")
    print(f"Read length: {read_length}, Coverage: {coverage}x")
    print(f"Simulation took {elapsed:.2f}s ({elapsed / len(reads) * 1000:.3f}ms per read)")

def main():
    # for all commands marked as required=False they must be changed to True
    parser = argparse.ArgumentParser(description='DNA Read Generator')

    parser.add_argument('command', choices=['run'],
                        help='Command to execute program')

    parser.add_argument('-i', '--input', required=False,
                        help='Input genome file (FASTA format)')

    parser.add_argument('-p', '--platform',
                        choices=['illumina', 'pacbio', 'nanopore'],
                        default='illumina',
                        help='Sequencing platform (default: illumina)')

    parser.add_argument('-l', '--read_length', type=int, default=100,
                        help='Read length in base pairs (default: 100)')

    parser.add_argument('-c', '--coverage', type=int, default=30,
                        help='Coverage depth (default: 30)')

    parser.add_argument('-f', '--format',
                        choices=['fasta', 'fastq', 'sam', 'bam', 'bed'],
                        default='fastq',
                        help='Output format (default: fastq)')

    parser.add_argument('-o', '--output', required=False,
                        help='Output file or directory path')

    parser.add_argument('-m', '--model',
                        choices=['art_illumina_profile_hiseq'],
                        help='Illumina read model (default: art_illumina_profile_hiseq)')

    args = parser.parse_args()

    # We should warn that they want to put in a profile or accept the default chosen option
    # if args.platform == 'illumina' and args.model is NONE:
    #     parser.error('-p illumina requires -m option')

    if args.command == 'run':
        if args.read_length == 150:  # default wasn't changed but nanopore and pacbio have longer reads
            if args.platform in ['pacbio', 'nanopore']:
                args.read_length = 10000

        run_example(
            coverage=args.coverage,
            read_length=args.read_length,
            output_format=args.format
        )
        print("reads generated")


if __name__ == "__main__":
    main()
