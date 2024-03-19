import os
import sys
import argparse
from pathlib import Path
from multiprocessing import cpu_count
from bin.fastq_analyzer import FastqAnalyzer
from bin.fasta_analyzer import FastaAnalyzer
from bin.annotator import GTFAnnotation
from bin.interval_analyzer import IntervalSummary
from bin.utility import *


@measure_execution_time
def main():
    num_cpus = max(cpu_count() - 1, 1)
    parser = argparse.ArgumentParser(description="Bioinformatics Toolbox Demo")
    parser.add_argument("-c", "--cpu", type=int, default=num_cpus,
                        help="Number of CPU cores to use for analysis (default: max CPU - 1)")
    parser.add_argument("-o", "--output", default=None, help="Specifies the path to the output directory. If not "
                                                             "provided, the output will be printed to the standard "
                                                             "output. When the command is 'annotation' and no "
                                                             "output directory is provided, the output file will be "
                                                             "placed in the working directory.")

    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # Subparser for the count_long_read command
    parser_count = subparsers.add_parser("count_reads_length_glt_threshold", help="Count reads with length greater then"
                                                                                  "the threshold in FASTQ files")
    parser_count.add_argument("-p", "--path", required=True, help="Path to a directory containing FASTQ files or a "
                                                                  "single FASTQ file")
    parser_count.add_argument("-t", "--threshold", type=int, default=30,
                              help="Threshold for minimum nucleotide length (default: 30)")

    # Subparser for the find_most_frequent_sequences command
    parser_freq = subparsers.add_parser("find_most_frequent_sequences", help="Find most frequent sequences in FASTQ"
                                                                             " files")
    parser_freq.add_argument("-p", "--path", required=True, help="Path to a directory containing FASTQ files or A "
                                                                 "single FASTQ file")
    parser_freq.add_argument("-n", "--n_top", type=int, default=10, help="Number of most frequent sequences to return")
    parser_freq.add_argument("-s", "--chuck_size", type=int, default=1, help="Number of processes submitted to the "
                                                                             "pool each time")

    # Subparser for the GTF annoatation command
    parser_anno = subparsers.add_parser("annotation", help="Annotate positions: [chr pos]")
    parser_anno.add_argument("-p", "--path", required=True, help="Path to a the file to annotate")
    parser_anno.add_argument("-k", "--knowledge_base", required=True, help="Path to the knowledge base providing the "
                                                                           "annotation")
    parser_anno.add_argument("-s", "--chuck_size", type=int, default=1, help="Number of processes submitted to the "
                                                                             "pool each time")
    # Subparser for the interval_summary command
    parser_interval = subparsers.add_parser("interval_summary", help="Summarize interval data by mean for each group "
                                                                     "of bins of specified column")
    parser_interval.add_argument("-p", "--path", required=True, help="Path to the interval data file")
    parser_interval.add_argument("-b", "--bins_num", type=int, default=10, help="Number of bins for %GC (default: 10)")
    parser_interval.add_argument("-g", "--group_by_key", default="%gc",
                                 help="Column name to group the data by (default: '%gc')")

    args = parser.parse_args()
    output = args.output

    # Check if no arguments were provided, then print usage
    if len(vars(args)) < 4:
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
        print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
        print(f"INFO: Find top {args.n_top}most frequent sequences in FASTA")
        analyzer = FastaAnalyzer(args.path, num_cpus=args.cpu, method="find_most_frequent_sequences", n_top=args.n_top,
                                 chuck_size=args.chuck_size)
        results = analyzer.analyze()
        output_file = None
        for seq_count in results:
            output_table = ["Rank\tSequence\tCount"]
            file_name = seq_count[0]
            output_table += [f"{ind + 1}\t{count[0]}\t{count[1]}" for ind, count in enumerate(seq_count[1])]
            output_table = "\n".join(output_table)
            if output:
                output_file = os.path.join(output, Path(file_name).stem + f".frequency_count_top_{args.n_top}.tsv")
                with open(output_file, "w+") as output_path:
                    print(output_table, file=output_path)
            else:
                print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
                print(f"Top {args.n_top} Most Frequent sequences")
                print(f"File: {file_name}:")
                print(output_table)
                print("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<")
        if output:
            print(f"INFO: output directory -> {output_file}")
            print("INFO: Analysis complete!")
            print("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<")
    elif args.command == "count_reads_length_glt_threshold":
        print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
        print(f"INFO: Count sequences with length larger than {args.threshold} in the input FASTQ")
        analyzer = FastqAnalyzer(args.path, num_cpus=args.cpu, method="calculate_long_sequences_percentage",
                                 threshold=args.threshold)
        results = analyzer.analyze()
        output_table = ["File\tPercentage\tThreshold"]
        for file, _, percentage, threshold in results:
            output_table.append(f"{file}\t{percentage:.2f}\t{threshold}")
        output_table = "\n".join(output_table)
        if output:
            output_file = os.path.join(output, f"fastq_reads_length_lgt_threshold.tsv")
            with open(output_file, "w+") as output_path:
                print(output_table, file=output_path)
        else:
            print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
            print(f"The percent of sequences in that file that are greater than {args.threshold} nucleotides long:")
            print(output_table)
            print("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<")
        if output:
            print(f"INFO: output directory -> {output}")
            print("INFO: Analysis complete!")
            print("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<")
    elif args.command == "annotation":
        if not output:
            output = "."  # default: output to work directory
        print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
        print(f"INFO: Annotate input file {args.path}")
        output_file = os.path.join(output, Path(args.path).stem + ".annotated.tsv")
        print(f"INFO: Output file -> {output_file}")
        annotator = GTFAnnotation(args.path, args.knowledge_base, output_file, num_cpus=args.cpu,
                                  chuck_size=args.chuck_size)
        annotator.annotate_all_positions()
        print("INFO: Annotation complete!")
        print("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<")
    elif args.command == "interval_summary":
        if not output:
            output = "."  # default: output to work directory
        print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
        print("INFO: Summarizing interval data")
        processor = IntervalSummary(args.path, args.group_by_key, args.bins_num)
        output_file = os.path.join(output, Path(args.path).stem + ".group_mean.tsv")
        processor.process_data(output_file)
        print(f"INFO: Mean for each {args.group_by_key} bin written to {output_file}")
        print("INFO: Interval summary complete!")
        print("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<")


if __name__ == "__main__":
    main()
