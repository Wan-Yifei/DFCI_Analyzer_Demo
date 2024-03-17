import os
import functools
import logging
from multiprocessing import Pool
from bin.constant import *

logger = logging.getLogger(__name__)
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
            logger.exception(f"Raised in {func.__name__}: {str(e)}\n"
                             f"------------------------------------------------")
            raise

    return wrapper


@log
def parallel_process(n_cpu, func, sources):
    """
    Execute a function concurrently using multiple processes.

    Parameters:
        n_cpu (int): Number of CPU cores to utilize.
        func (function): The function to be executed concurrently.
        sources (list): List of input data to be processed by the function.

    Returns:
        list: List containing the results of applying the function to each input in sources.
    """
    with Pool(n_cpu) as pool:
        # Map the function over the input data and execute it concurrently using multiple processes
        results = pool.map(func, sources)
        return results


@log
def find_input_files(method, input_path):
    """
    Finds all input files in the specified directory.

    Returns:
        list: A list of paths to input files.
    """
    suffix = filename_suffix[method]
    input_files = []
    if os.path.isdir(input_path):
        for root, dirs, files in os.walk(input_path):
            for file in files:
                if file.endswith(suffix):
                    input_files.append(os.path.join(root, file))
    elif os.path.isfile(input_path) and input_path.endswith(suffix):
        input_files.append(input_path)
    if not input_files:
        raise ValueError(f"Invalid path provided or file is not a valid type for the chosen method, '{suffix}' files!"
                         f" are required by {method}.")
    for file in input_files:
        if os.path.getsize(file) == 0:
            print(f"Warning: File '{input_path}' has zero file size.")
    return input_files


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


def fast_file_reader(file_path, row_count_limit=10000):
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