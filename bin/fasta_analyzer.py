from collections import Counter
from Bio import SeqIO
from bin.utility import *


class FastaAnalyzer:
    @log
    def __init__(self, directory, num_cpus, method, n_top=10, chuck_size=1):
        """
        Initializes FastaAnalyzer class.

        Args:
            directory (str): Directory containing FASTA files.
            num_cpus (int): Number of CPUs to utilize for processing.
            method (str): Method for analysis.
            batch_size (int): Size of each batch of sequences (default is subfile_row_count_per_file).
            n_top (int): Number of top frequent sequences to return (default is 10).
            turn_on_split_file (bool): Flag indicating whether to split the input files (default is False).
        """
        self.directory = directory
        self.num_cpus = num_cpus
        self.method = method
        self.n_top = n_top
        self.chuck_size = chuck_size
        self.fasta_files = find_input_files(self.method, self.directory)

    @log
    def find_most_frequent_sequences(self, fasta_file_path):
        """
        Finds the most frequent sequences among multiple sequence counters.

        Args:
            file_name (str): Name of the FASTA file being processed.
            sequence_counter_iter (iterable): Iterable of sequence counters.

        Returns:
            tuple: A tuple containing the file name and a list of tuples of the most frequent sequences and their counts.
            :param file_name:
            :param sequence_counter:
        """
        if isinstance(fasta_file_path, str):
            sequences = SeqIO.parse(fasta_file_path, "fasta")  # Parsing FASTA file
        else:
            raise ValueError("Input is not a single FASTQ file or a batch of SeqIO records!!")
        sequence_counter = Counter(str(record.seq) for record in sequences)  # Counting sequences
        most_common = sequence_counter.most_common(self.n_top)  # Getting most common sequences
        return fasta_file_path, most_common  # Returning list of most common sequences and their counts

    @log
    def analyze(self):
        """
        Analyzes FASTA files using multiprocessing.

        Returns:
            list: A list of results.
        """

        results = parallel_process(n_cpu=self.num_cpus, func=getattr(self, self.method), sources=self.fasta_files,
                                   lazy_map=True, chunk_size=self.chuck_size)
        return results
