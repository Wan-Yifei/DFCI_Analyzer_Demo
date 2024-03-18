import tempfile
from collections import Counter
from Bio import SeqIO
from bin.utility import *


@log
def fast_file_reader(file_path, batch_size=subfile_count_per_file, turn_on_split=False):
    """
    Split file.
    :param batch_size:
    :param turn_on:
    :param file_path:
    :param row_count_limit:
    :return: list of temp files.
    """
    temp_files = []
    if not turn_on_split:
        temp_files.append(file_path)
        return temp_files
    record_iter = SeqIO.parse(open(file_path), "fasta")
    temp_dir = tempfile.mkdtemp()
    print("Write temp FASTQ files:")
    for i, batch in enumerate(batch_iterator(record_iter, batch_size)):
        temp_file = os.path.join(temp_dir, "group_%i.fasta" % (i + 1))
        with open(temp_file, "w") as handle:
            count = SeqIO.write(batch, handle, "fasta")
        print("Wrote %i records to %s" % (count, temp_file))
        temp_files.append(temp_file)
    return temp_files


class FastaAnalyzer:
    def __init__(self, directory, num_cpus, method, n_top=10, turn_on_split_file=False):
        """
        Initializes FastaAnalyzer class.

        Args:
            directory (str): Directory containing FASTA files.
            num_cpus (int): Number of CPUs to utilize for processing.
            method (str): Method for analysis.
            n_top (int): Number of top frequent sequences to return (default is 10).
        """
        self.directory = directory
        self.num_cpus = num_cpus
        self.method = method
        self.n_top = n_top
        self.turn_on_split_file = turn_on_split_file
        self.fasta_files = find_input_files(self.method, self.directory)  # Finding input FASTA files

    @log
    def count_sequences(self, fasta_file_path):
        """
        Finds the most frequent sequences in a given FASTA file.

        Args:
            fasta_file_path (str): Path to the FASTA file.

        Returns:
            list: A list of tuples containing the most frequent sequences and their counts.
        """
        sequences = SeqIO.parse(fasta_file_path, "fasta")  # Parsing FASTA file
        sequence_counter = Counter(str(record.seq) for record in sequences)  # Counting sequences
        return sequence_counter

    @log
    def find_most_frequent_sequences(self, file_name, sequence_counter_iter):
        merged_counter = Counter()
        for seq_counter in sequence_counter_iter:
            merged_counter.update(dict(seq_counter))
        most_common = merged_counter.most_common(self.n_top)  # Getting most common sequences
        return file_name, most_common  # Returning list of most common sequences and their counts

    def downstream_trigger(self):
        downstream_map = {
            "count_sequences": "find_most_frequent_sequences"
        }
        return getattr(self, downstream_map[self.method])


    @log
    def analyze(self):
        """
        Analyzes FASTA files using multiprocessing.

        Returns:
            list: A list of resultsd.
        """
        summary = []
        for fasta in self.fasta_files:
            temp_fasta_files = fast_file_reader(fasta, turn_on_split=self.turn_on_split_file)
            num_processes = min(self.num_cpus, len(temp_fasta_files))
            results = parallel_process(n_cpu=num_processes, func=getattr(self, self.method), sources=temp_fasta_files)
            summary.append(self.downstream_trigger()(fasta, results))
        return summary
