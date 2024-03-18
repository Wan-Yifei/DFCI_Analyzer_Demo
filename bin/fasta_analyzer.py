import tempfile
from collections import Counter
from Bio import SeqIO
from bin.utility import *




class FastaAnalyzer:
    @log
    def __init__(self, directory, num_cpus, method, batch_size=subfile_row_count_per_file, n_top=10, turn_on_split_file=False):
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
        self.batch_size = batch_size
        self.turn_on_split_file = turn_on_split_file
        self.fasta_files = find_input_files(self.method, self.directory)  # Finding input FASTA files

    @log
    def fast_file_reader(self, file_path, turn_on_split=False):
        """
        Reads a FASTA file and splits it into smaller temporary files.

        Args:
            file_path (str): Path to the input FASTA file.
            turn_on_split (bool): Flag indicating whether to split the file (default is False).

        Returns:
            list: List of paths to temporary files containing batches of sequences.
        """
        temp_files = []
        if not turn_on_split:
            temp_files.append(file_path)
            return temp_files
        record_iter = SeqIO.parse(open(file_path), "fasta")
        temp_dir = tempfile.mkdtemp()
        print("Write temp FASTQ files:")
        for i, batch in enumerate(batch_iterator(record_iter, self.batch_size)):
            temp_file = os.path.join(temp_dir, "group_%i.fasta" % (i + 1))
            with open(temp_file, "w") as handle:
                count = SeqIO.write(batch, handle, "fasta")
            print("Wrote %i records to %s" % (count, temp_file))
            temp_files.append(temp_file)
        return temp_files

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
        """
        Finds the most frequent sequences among multiple sequence counters.

        Args:
            file_name (str): Name of the FASTA file being processed.
            sequence_counter_iter (iterable): Iterable of sequence counters.

        Returns:
            tuple: A tuple containing the file name and a list of tuples of the most frequent sequences and their counts.
        """
        merged_counter = Counter()
        for seq_counter in sequence_counter_iter:
            merged_counter.update(dict(seq_counter))
        most_common = merged_counter.most_common(self.n_top)  # Getting most common sequences
        return file_name, most_common  # Returning list of most common sequences and their counts

    @log
    def downstream_trigger(self):
        """
        Gets the downstream method based on the selected analysis method.

        Returns:
            method: The downstream method to execute.
        """
        downstream_map = {
            "count_sequences": "find_most_frequent_sequences"
        }
        return getattr(self, downstream_map[self.method])

    @log
    def analyze(self):
        """
        Analyzes FASTA files using multiprocessing.

        Returns:
            list: A list of results.
        """
        summary = []
        for fasta in self.fasta_files:
            temp_fasta_files = self.fast_file_reader(fasta, turn_on_split=self.turn_on_split_file)
            num_processes = min(self.num_cpus, len(temp_fasta_files))
            results = parallel_process(n_cpu=num_processes, func=getattr(self, self.method), sources=temp_fasta_files)
            summary.append(self.downstream_trigger()(fasta, results))
        return summary
