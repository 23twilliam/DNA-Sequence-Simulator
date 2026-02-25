import argparse
import io_functions
from simulators import RandomSimulator

def run_example(read_length, coverage, output_format):
    """The original run_example() refactored to use the simulator class"""
    example_input = "../Example_FASTA/gene.fna"

    # Create simulator instance
    simulator = RandomSimulator(
        read_length=read_length,
        coverage=coverage,
        error_rate=0.001
    )

    # Generate reads
    reads = simulator.simulate_from_file(example_input, io_functions)

    # Write output
    if output_format == "fastq":
        io_functions.write_fastq(reads, f"simulated_reads.{output_format}")
    elif output_format == "fasta":
        io_functions.write_fasta(reads, f"simulated_reads.{output_format}")
    else:
        print("error")


    print(f"Generated {len(reads)} reads using {simulator.get_platform_name()} simulator")
    print(f"Read length: {read_length}, Coverage: {coverage}x")

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

    parser.add_argument('-l', '--read_length', type=int, default=150,
                        help='Read length in base pairs (default: 150)')

    parser.add_argument('-c', '--coverage', type=int, default=30,
                        help='Coverage depth (default: 30)')

    parser.add_argument('-f', '--format',
                        choices=['fasta', 'fastq', 'sam', 'bam', 'bed'],
                        default='fastq',
                        help='Output format (default: fastq)')

    parser.add_argument('-o', '--output', required=False,
                        help='Output file or directory path')


    args = parser.parse_args()

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
