import csv
import ast
import requests
from requests.exceptions import RequestException
from utility import *
from constant import *


class EnsemblVariantInfoRetriever:
    """A class to retrieve variant information from Ensembl databases."""

    @log
    def __init__(self, species, input_file, output_file, skip_bad_call=False):
        """
        Initialize the EnsemblVariantInfoRetriever object.

        Args:
            species (str): The species for which variant information is retrieved.
            input_file (str): Path to the input file containing variant IDs.
            output_file (str): Path to the output file to write the variant information.
            skip_bad_call (bool, optional): If True, skip variant IDs for which information
                cannot be retrieved. Defaults to False.
        """
        self.species = species
        self.input_file = input_file
        self.output_file = output_file
        self.variant_ids = self.read_variant_ids()
        self.skip_bad_call = skip_bad_call
        check_path_is_file(self.input_file)

    @log
    def read_variant_ids(self):
        """
        Read variant IDs from the input file.

        Returns:
            list: A list of variant IDs.
        """
        variant_ids = []
        with open(self.input_file, 'r') as file:
            reader = csv.reader(file)
            for row in reader:
                variant_ids.extend(row)
        return variant_ids

    @log
    def get_variant_info(self, variant_id):
        """
        Retrieve variant information from Ensembl Variant Effect Predictor (VEP) API.

        Args:
            variant_id (str): The variant ID to retrieve information for.

        Returns:
            dict: Variant information in JSON format.
        """
        try:
            url = VEP_API_URL.format(self.species, variant_id)
            response = requests.get(url, headers={"Content-Type": "application/json"})
            response.raise_for_status()  # Raise an exception for HTTP errors
            return response.json()
        except requests.exceptions.RequestException as e:
            if self.skip_bad_call:
                return None
            if "No variant found with ID" in response.text:
                print(f"\nWARNING: Variant Info API - {ast.literal_eval(response.text)['error']}")
                return None
            raise RequestException(f"Error fetching variant info for {variant_id}: {e}")

    @log
    def get_transcript_info(self, transcript_id):
        """
        Retrieve transcript information from Ensembl Transcript API.

        Args:
            transcript_id (str): The transcript ID to retrieve information for.

        Returns:
            dict: Transcript information in JSON format.
        """
        try:
            url = TRANSCRIPT_API_URL.format(transcript_id)
            response = requests.get(url, headers={"Content-Type": "application/json"})
            response.raise_for_status()  # Raise an exception for HTTP errors
            return response.json()
        except requests.exceptions.RequestException as e:
            if self.skip_bad_call:
                return None
            if "not found" in response.text:
                print(f"\nWARNING: Transcript Info API - {ast.literal_eval(response.text)['error']}")
                return None
            raise RequestException(f"Error fetching transcript info for {transcript_id}: {e}")


    @log
    def fetch_and_write_variant_info(self):
        """
        Fetch variant information for each variant ID and write to the output file.
        """
        with open(self.output_file, 'w', newline='') as csvfile:
            csv_writer = csv.writer(csvfile, delimiter='\t')
            csv_writer.writerow(["Variant ID", "Alleles", "Location", "Effects", "Gene"])

            for variant_id in tqdm(self.variant_ids):
                variant_info = self.get_variant_info(variant_id)
                if variant_info:
                    for colocated_variant in variant_info:
                        alleles = colocated_variant.get("allele_string")
                        location = colocated_variant.get("seq_region_name")
                        effects = colocated_variant.get("most_severe_consequence")
                        transcript_id = colocated_variant.get("transcript_consequences")[0].get("transcript_id")

                        transcript_info = self.get_transcript_info(transcript_id)
                        if transcript_info:
                            gene_name = transcript_info.get("display_name")
                            csv_writer.writerow([variant_id, alleles, location, effects, gene_name])
                        else:
                            csv_writer.writerow([variant_id, ".", ".", ".", "."])
                else:
                    csv_writer.writerow([variant_id, ".", ".", ".", "."])

        print(f"Data has been written to {self.output_file}")
