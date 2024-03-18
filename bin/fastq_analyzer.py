from bin.utility import *
from Bio import SeqIO

logger = logging.getLogger("FastqAnalyzer")


class FastqAnalyzer(object):
    """
    Analyzes FASTQ files and performs various analyses.
    """
    @log
    def __init__(self, directory, num_cpus, method, threshold=None, n_top=None):
        """
        Initializes the FastqAnalyzer instance.

        Args:
            directory (str): The directory containing FASTQ files or a single FASTQ file path.
            num_cpus (int): Number of CPU cores to utilize for multiprocessing.
            method (str): Method for calculating long sequences percentage.
            threshold (int): Minimum length threshold for considering a sequence as long.
            n_top (int): Number of most frequent sequences to return.
        Raises:
            ValueError: If the provided path is invalid or the file is not a FASTQ file.
        """
        self.directory = directory
        self.fastq_files = []
        self.num_cpus = num_cpus
        self.method = method
        self.threshold = threshold
        self.n_top = n_top
        self.fastq_files = find_input_files(self.method, self.directory)

    @log
    def calculate_long_sequences_percentage(self, file_path):
        """
        Calculates the percentage of sequences longer than the threshold in a FASTQ file.

        Args:
            file_path (str): Path to the FASTQ file.

        Returns:
            tuple: A tuple containing file path, total sequence count, percentage of long sequences, and the threshold.
        """
        long_seq_count = 0
        total_seq_count = 0

        with open(file_path, "r") as handle:
            for record in SeqIO.parse(handle, "fastq"):
                total_seq_count += 1
                if len(record.seq) > self.threshold:
                    long_seq_count += 1
        if total_seq_count == 0:
            return file_path, total_seq_count, 0, self.threshold
        else:
            percentage = (long_seq_count / total_seq_count) * 100
            return file_path, total_seq_count, percentage, self.threshold

    @log
    def analyze(self):
        """
        Analyzes FASTQ files using multiprocessing.

        Returns:
            list: A list of results containing file path, total sequence count, percentage of long sequences, and threshold.
        """
        num_processes = min(self.num_cpus, len(self.fastq_files))
        return parallel_process(n_cpu=num_processes, func=getattr(self, self.method), sources=self.fastq_files)
