import os
import functools
import logging
import timeit
from tqdm import tqdm
from pathlib import Path
from multiprocessing import Pool, cpu_count
from bin.constant import *

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
console_handler = logging.StreamHandler()
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)


def measure_execution_time(func):
    """
    Decorator to measure the execution time of a function.

    Args:
        func (function): The function to be decorated.

    Returns:
        function: The wrapped function.
    """

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        """
        Wrapper function that measures the execution time of the decorated function.

        Args:
            *args: Variable length argument list.
            **kwargs: Arbitrary keyword arguments.

        Returns:
            Any: The result of the decorated function.
        """
        start_time = timeit.default_timer()
        result = func(*args, **kwargs)
        end_time = timeit.default_timer()
        execution_time = end_time - start_time
        print(f"------------------------------------------------")
        print(f"Execution time of {func.__name__}: {execution_time:.6f} seconds")
        print(f"------------------------------------------------")
        return result

    return wrapper


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
def parallel_process(n_cpu, func, sources, lazy_map=False, chunk_size=1):
    """
    Execute a function concurrently using multiple processes.

    Parameters:
        n_cpu (int): Number of CPU cores to utilize.
        func (function): The function to be executed concurrently.
        sources (list): List of input data to be processed by the function.
        lazy_map (bool): If True, uses lazy evaluation with imap. (default is False)
        chunk_size (int): Number of items to send to the worker process at a time. (default is 2)

    Returns:
        list: List containing the results of applying the function to each input in sources.
    """
    with Pool(n_cpu) as pool:
        if lazy_map:
            results = list(tqdm(pool.imap(func, sources, chunk_size), total=len(sources)))
        else:
            # Map the function over the input data and execute it concurrently using multiple processes
            results = pool.map(func, sources)
        return results


@log
def find_input_files(method, input_path, is_file_flag=False):
    """
    Finds all input files in the specified directory.

    Returns:
        list: A list of paths to input files.
    """
    suffix = filename_suffix[method]
    input_files = []
    assert os.path.exists(input_path), f"Input path {input_path} doesn't Exist!!"
    if os.path.isdir(input_path):
        assert not is_file_flag, f"Input path {input_path} is required as a single file, please check!!"
        for root, dirs, files in os.walk(input_path):
            for file in files:
                if file.endswith(suffix):
                    input_files.append(os.path.join(root, file))
    elif os.path.isfile(input_path) and input_path.endswith(suffix):
        input_files.append(input_path)
    if not input_files:
        raise FileNotFoundError(f"Invalid path f{input_path} provided or file is not a valid type for the chosen "
                                f"method, '{suffix}' files! are required by {method}.")
    for file in input_files:
        if os.path.getsize(file) == 0:
            print(f"Warning: File '{input_path}' has zero file size.")
    return input_files


@log
def check_path_is_file(file_path):
    """
    Check if the given path exists and is a file.

    Args:
    - path: A string representing the path to be checked.

    Returns:
    - The path if it exists and is a file, otherwise None.
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"The path '{file_path}' does not exist.")
    elif not os.path.isfile(file_path):
        raise NotADirectoryError(f"The path '{file_path}' exists but is not a file.")


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


def make_output_path(output_dir, input_file, suffix_extension):
    """
    Create an output path for a file based on the input file path and the specified output directory and suffix extension.

    Args:
        output_dir (str): The directory where the output file will be placed.
        input_file (str): The path to the input file.
        suffix_extension (str): The suffix and extension to be appended to the input file's stem.

    Returns:
        str: The output file path.

    Example:
        If output_dir='/output', input_file='/path/to/input.txt', and suffix_extension='.output.txt',
        then the output file path will be '/output/input.output.txt'.
    """
    return os.path.join(output_dir, Path(input_file).stem + suffix_extension)
