import pandas as pd
import numpy as np
from collections import defaultdict
from bin.constant import *
from bin.utility import *


class GTFAnnotation:
    """
    A class for annotating positions in a genomic dataset using GTF annotation data.

    Parameters:
        knowledge_file (str): Path to the GTF annotation file.
        input_file (str): Path to the input genomic data file.
        output_file (str): Path to the output file where annotated data will be written.
        num_cpus (int): Number of CPU cores to use for parallel processing.
        chuck_size (int, optional): Chunk size for parallel processing. Default is 1.

    Attributes:
        knowledge_df (DataFrame): DataFrame containing GTF annotation data.
        input_positions_df (DataFrame): DataFrame containing input genomic data.
        output_file (str): Path to the output file.
        num_cpus (int): Number of CPU cores to use.
        chuck_size (int): Chunk size for parallel processing.

    """
    @log
    def __init__(self, input_file, knowledge_file, output_file, num_cpus, chuck_size=1):
        self.input_knowledge_file = find_input_files("annotation", knowledge_file, is_file_flag=True)[0]
        self.knowledge_df = self.load_annotation(self.input_knowledge_file)
        self.input_positions_df = pd.read_csv(input_file, sep='\t', header=None, names=['chromosome', 'position'])
        self.output_file = output_file
        self.num_cpus = num_cpus
        self.chuck_size = chuck_size

    @log
    def load_annotation(self, annotation_file):
        """
        Load GTF annotation data into memory.

        Parameters:
            annotation_file (str): Path to the GTF annotation file.

        Returns:
            DataFrame: DataFrame containing GTF annotation data.

        """
        knowledge_df = pd.read_csv(annotation_file, comment='#', header=None, sep='\t', index_col=False,
                                   names=gtf_columns)
        return knowledge_df

    @log
    def extract_annotation(self, attributes_field):
        """
        Extract annotation from GTF attributes field.

        Parameters:
            attributes_field (str): Attributes field from GTF annotation.

        Returns:
            str: Extracted annotation.

        """
        fields = attributes_field.rstrip(";").split(';')
        annotation_dict = defaultdict(lambda: ".")
        for field in fields:
            annotation_dict[field.split('"')[0].strip()] = field.split('"')[1]
        annotation = "\t".join([annotation_dict[key] for key in annotation_fields])
        return annotation

    @log
    def annotate_position(self, chromosome, position):
        """
        Annotate a single genomic position.

        Parameters:
            chromosome (str): Chromosome identifier.
            position (int): Genomic position.

        Returns:
            str: Annotated position.

        """
        # Filter annotation DataFrame for the given chromosome
        overlapping_genes = self.knowledge_df[
            (self.knowledge_df['chr'] == chromosome) &
            (self.knowledge_df['start'] <= position) &
            (self.knowledge_df['end'] >= position)
            ]

        # Extract annotation
        annotation_rows = overlapping_genes['attributes'].apply(self.extract_annotation).tolist()
        if not annotation_rows:
            annotation_row = str(chromosome) + "\t" + str(position) + "\t" + "\t".join(["." for i in annotation_fields])
        else:
            annotation_rows = [str(chromosome) + "\t" + str(position) + "\t" + row for row in annotation_rows]
            annotation_row = "\n".join(annotation_rows)
        return annotation_row

    @log
    def execute_annotation(self, dataframe):
        """
        Perform annotation on a DataFrame containing genomic positions.

        Parameters:
            dataframe (DataFrame): DataFrame containing 'chromosome' and 'position' columns.

        Returns:
            DataFrame: DataFrame with 'annotation' column added.

        """
        dataframe["annotation"] = np.vectorize(self.annotate_position)(dataframe['chromosome'], dataframe['position'])
        return dataframe

    @log
    def annotate_all_positions(self):
        """
        Annotate all genomic positions and write to the output file.

        """
        num_cpus = max(self.num_cpus - 1, 1)
        sub_df = np.array_split(self.input_positions_df, num_cpus)
        results = parallel_process(n_cpu=num_cpus, func=self.execute_annotation, sources=sub_df,
                                   lazy_map=True, chunk_size=self.chuck_size)
        annotated_df = pd.concat(results)
        with open(self.output_file, 'w') as f:
            header = "\t".join(["chr", "pos"] + annotation_fields)
            f.write(header + '\n')
            for index, annotation_row in annotated_df.iterrows():
                for row in annotation_row['annotation'].strip().split("\n"):
                    f.write(f"{row}\n")
