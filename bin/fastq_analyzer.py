import os
import logging
import functools
from collections import Counter
from multiprocessing import Pool
from Bio import SeqIO

filename_suffix = {
    "find_most_frequent_sequences": ".fasta",
    "calculate_long_sequences_percentage": ".fastq"
}

logger = logging.getLogger("FastqAnalyzer")
logger.setLevel(logging.DEBUG)
console_handler = logging.StreamHandler()
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)


def log(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        try:
            result = func(*args, **kwargs)
            return result
        except Exception as e:
            logger.exception(f">>>>>>>>>>>>>>>\n Exception raised in {func.__name__}: {str(e)}\n<<<<<<<<<<<<<<<<<<<<<<<"
                             f"<<<<<<<<<<<<<<<<<<<<<<<<<<")
            raise

    return wrapper


# Not use
# for potential improvement: split file for parallel
def batch_iterator(iterator, batch_size):
    """
    Returns lists of length batch_size.
    """
    batch = []
    for entry in iterator:
        batch.append(entry)
        if len(batch) == batch_size:
            yield batch
            batch = []
    if batch:
        yield batch


# Not use
# for potential improvement: split file for parallel
def file_reader(file_path, row_count_limit=10000):
    """
    Split file.
    :param file_path:
    :param row_count_limit:
    :return:
    """
    for i, batch in enumerate(batch_iterator(file_path, row_count_limit)):
        filename = "group_%i.fastq" % (i + 1)
        with open(filename, "w") as handle:
            count = SeqIO.write(batch, handle, "fastq")
        print("Wrote %i records to %s" % (count, filename))


class FastqAnalyzer(object):
    """
    Analyzes FASTQ files and performs various analyses.
    """

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
        self.fastq_files = []
        self.find_fastq_files()

    @log
    def find_fastq_files(self):
        """
        Finds all FASTQ files in the specified directory.

        Returns:
            list: A list of paths to FASTQ files.
        """
        suffix = filename_suffix[self.method]
        if os.path.isdir(self.directory):
            for root, dirs, files in os.walk(self.directory):
                for file in files:
                    if file.endswith(suffix):
                        self.fastq_files.append(os.path.join(root, file))
        elif os.path.isfile(self.directory) and self.directory.endswith(suffix):
            self.fastq_files.append(self.directory)
        if not self.fastq_files:
            raise ValueError(f"Invalid path provided or file is not a valid type for the chosen method, {suffix} files!"
                             f" are required by {self.method}.")
        for file in self.fastq_files:
            if os.path.getsize(file) == 0:
                logger.warning(f"File '{self.directory}' has zero file size.")

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
    def find_most_frequent_sequences(self, fasta_file_path):
        """
        Finds the most frequent sequences in a given FASTA file.

        Args:
            fasta_file_path (str): Path to the FASTA file.

        Returns:
            list: A list of tuples containing the most frequent sequences and their counts.
        """
        sequences = SeqIO.parse(fasta_file_path, "fasta")
        sequence_counter = Counter(str(record.seq) for record in sequences)
        most_common = sequence_counter.most_common(self.n_top)
        return fasta_file_path, most_common

    @log
    def analyze(self):
        """
        Analyzes FASTQ files using multiprocessing.

        Returns:
            list: A list of results containing file path, total sequence count, percentage of long sequences, and threshold.
        """
        num_processes = min(self.num_cpus, len(self.fastq_files))
        with Pool(num_processes) as pool:
            results = pool.map(getattr(self, self.method), self.fastq_files)
            return results
