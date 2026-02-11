import read_generation
import argparse


def main():
    parser = argparse.ArgumentParser(description='DNA Read Generator')
    parser.add_argument('command', choices=['run'],
                        help='Command to execute')

    args = parser.parse_args()

    if args.command == 'run':
        read_generation.run_example()
        print("program ran")


if __name__ == "__main__":
    main()