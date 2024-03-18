import os
import sys
import argparse
from pathlib import Path
from multiprocessing import cpu_count
from bin.fastq_analyzer import FastqAnalyzer
from bin.fasta_analyzer import FastaAnalyzer
from bin.utility import *


@measure_execution_time
def main():
    num_cpus = max(cpu_count() - 1, 1)
    parser = argparse.ArgumentParser(description="Bioinformatics Toolbox Demo")
    parser.add_argument("-c", "--cpu", type=int, default=num_cpus,
                        help="Number of CPU cores to use for analysis (default: max CPU - 1)")
    parser.add_argument("-p", "--path", help="Path to a directory containing FASTQ/FASTA files or a single FASTQ/FASTA"
                                             "FASTA file")
    parser.add_argument("-o", "--output", default=None, help="Path to output directory, if not provided, print output "
                                                             "to stdout")

    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # Subparser for the analyze command
    parser_analyze = subparsers.add_parser("analyze", help="Analyze Bioinformatics' files")
    parser_analyze.add_argument("path", help="Path to a directory containing FASTQ files or a single FASTQ file")

    # Subparser for the count_long_read command
    parser_count = subparsers.add_parser("count_reads_length_glt_threshold", help="Count reads with length greater then"
                                                                                  "the threshold in FASTQ files")
    parser_count.add_argument("-t", "--threshold", type=int, default=30,
                              help="Threshold for minimum nucleotide length (default: 30)")

    # Subparser for the find_most_frequent_sequences command
    parser_freq = subparsers.add_parser("find_most_frequent_sequences", help="Find most frequent sequences in FASTQ"
                                                                             " files")
    parser_freq.add_argument("-n", "--n_top", type=int, default=10, help="Number of most frequent sequences to return")
    parser_freq.add_argument("-k", "--chuck_size", type=int, default=1, help="Number of processes submitted to the "
                                                                              "pool each time")
    args = parser.parse_args()
    output = args.output

    # Check if no arguments were provided, then print usage
    if len(vars(args)) < 5:
        print("Error: Miss required arguments. Please check!!")
        parser.print_help()
        return

    # Check the path to output dir
    if output:
        try:
            assert os.path.isdir(output)
        except AssertionError:
            print("Error: input output path must be a directory!!")
            sys.exit(1)

    if args.command == "find_most_frequent_sequences":
        analyzer = FastaAnalyzer(args.path, num_cpus=num_cpus, method="find_most_frequent_sequences", n_top=args.n_top,
                                 chuck_size=args.chuck_size)
        results = analyzer.analyze()
        for seq_count in results:
            output_table = ["Rank\tSequence\tCount"]
            file_name = seq_count[0]
            output_table += [f"{ind + 1}\t{count[0]}\t{count[1]}" for ind, count in enumerate(seq_count[1])]
            output_table = "\n".join(output_table)
            if output:
                with open(os.path.join(output, Path(file_name).stem + f".frequency_count_top_{args.n_top}.tsv"),
                          "w+") as output_path:
                    print(output_table, file=output_path)
            else:
                print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
                print(f"Top {args.n_top} Most Frequent sequences")
                print(f"File: {file_name}:")
                print(output_table)
                print("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<")
    elif args.command == "count_reads_length_glt_threshold":
        analyzer = FastqAnalyzer(args.path, num_cpus=num_cpus, method="calculate_long_sequences_percentage",
                                 threshold=args.threshold)
        results = analyzer.analyze()
        output_table = ["File\tPercentage\tThreshold"]
        for file, _, percentage, threshold in results:
            output_table.append(f"{file}\t{percentage:.2f}\t{threshold}")
        output_table = "\n".join(output_table)
        if output:
            with open(os.path.join(output, f"fastq_reads_length_lgt_threshold.tsv"),
                      "w+") as output_path:
                print(output_table, file=output_path)
        else:
            print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
            print(f"The percent of sequences in that file that are greater than {args.threshold} nucleotides long:")
            print(output_table)
            print("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<")


if __name__ == "__main__":
    main()
