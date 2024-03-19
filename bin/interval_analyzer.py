import pandas as pd
import tempfile
from bin.utility import *


class IntervalSummary:
    """
    A class to summarize interval data based on specified parameters.

    Attributes:
    - file_path (str): The path to the input file containing interval data.
    - bins_num (int): The number of bins to use for grouping.
    - group_by_key (str): The column name to group the data by.
    """

    @log
    def __init__(self, input_path, group_by_key, bins_num=10):
        """
        Initialize the IntervalSummary object.

        Parameters:
        - input_path (str): The path to the input file containing interval data.
        - bins_num (int, optional): The number of bins to use for grouping. Default is 10.
        - group_by_key (str, optional): The column name to group the data by. e.g. '%gc'.
        """
        self.file_path = input_path
        self.bins_num = bins_num
        self.group_by_key = group_by_key
        check_path_is_file(self.file_path)

    @log
    def pre_fix_input(self):
        """
        Pre-parse the input file, removing empty columns and writing to a temporary CSV file.

        Returns:
        - str: The path to the temporary CSV file.
        """
        temp_csv_path = None
        with open(self.file_path, 'r') as input_file, \
                tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.csv') as temp_file:
            for line in input_file:
                # Split the line by tabs
                columns = line.strip().split('\t')

                # Remove empty columns
                columns = [col for col in columns if col.strip()]

                # Join the columns back with tabs
                pre_parsed_line = '\t'.join(columns)

                # Write the pre-parsed line to the temporary file
                temp_file.write(pre_parsed_line)
                temp_file.write('\n')  # Add a new line after each line

            temp_csv_path = temp_file.name

        return temp_csv_path

    def process_data(self, output_file):
        """
        Process the interval data to calculate mean coverage for each group and write the results to a CSV file.

        Parameters:
        - output_file (str): The path to the output CSV file.
        """
        temp_csv_path = self.pre_fix_input()

        # Read the pre-parsed data into a DataFrame
        data = pd.read_csv(temp_csv_path, delimiter="\t")

        # Define bins for grouping
        bins = [(i / self.bins_num) for i in range(self.bins_num + 1)]

        # Generate the column name for bins dynamically
        bins_column = f"{self.group_by_key}_bins"

        # Group the data by specified key and calculate mean coverage for each bin
        data[bins_column] = pd.cut(data[self.group_by_key], bins)
        grouped = data.groupby(bins_column)['mean_coverage'].mean().reset_index()
        grouped.rename(columns={'mean_coverage': 'mean_target_coverage', bins_column: self.group_by_key}, inplace=True)

        # Write the results to a CSV file
        grouped.to_csv(output_file, index=False, sep="\t", na_rep='NA')
