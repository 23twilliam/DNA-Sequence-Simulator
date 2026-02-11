import read_generation
import argparse

def main():
    parser = argparse.ArgumentParser(description='DNA Read Generator')
    parser.add_argument('command', choices=['run'],
                        help='Command to execute program')

    # Defaults move here from function
    parser.add_argument('-l', '--read_length', type=int, default=150,
                        help='Read length in base pairs (default: 150)')
    parser.add_argument('-c', '--coverage', type=int, default=30,
                        help='Coverage depth (default: 30)')
    args = parser.parse_args()

    if args.command == 'run':
        read_generation.run_example(coverage=args.coverage, read_length=args.read_length)
        print("reads generated")


if __name__ == "__main__":
    main()