filename_suffix = {
    "find_most_frequent_sequences": ".fasta",
    "calculate_long_sequences_percentage": ".fastq",
    "annotation": ".gtf"
}

# Single input file only cmd

# Annotation
gtf_columns = ['chr', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
annotation_fields = ["gene_id", "transcript_id", "exon_number ", "exon_id", "gene_name"]

# REST API URLs
VEP_API_URL = "https://rest.ensembl.org/vep/{}/id/{}"
TRANSCRIPT_API_URL = "https://rest.ensembl.org/lookup/id/{}?expand=1"