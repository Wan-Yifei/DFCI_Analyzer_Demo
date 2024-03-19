# Bioinformatics Toolbox (Demo)

## Description
The Bioinformatics Toolbox is a collection of command-line tools designed to perform various bioinformatics analyses. It includes functionalities such as analyzing FASTQ files, annotating genomic positions, summarizing interval data, fetching variant information from Ensembl, and more.

****
## Features
- **count_reads_length_glt_threshold**: Count reads with a length greater than the specified threshold in FASTQ files.
- **find_most_frequent_sequences**: Find the most frequent sequences in FASTQ files.
- **annotation**: Annotate genomic positions using a knowledge base.
- **interval_summary**: Summarize interval data by calculating the mean for a specific group of bins.
- **variant_info**: Fetch variant information from Ensembl for provided variant IDs.
****
## Prerequisites
- Python 3.8 or higher
- pip 24.0 or higher
****
## Installation

```bash
# Clone this repository:
git clone https://github.com/Wan-Yifei/DFCI_Analyzer_Demo.git

# Navigate to the project directory:
cd DFCI_Analyzer_Demo

# Install the required dependencies and setup
pip install .

# Verify the installation
dfci-analysis -h
```
****
## Basic Usage

### Q1. 
**Recursively find all FASTQ files in a directory and report each file name and the percent of 
sequences in that file that are greater than 30 nucleotides long.**
  
- Source code:
The module used for counting sequences can be found in `/bin/fastq_analyzer.py`.
```bash
# The default threshold is set to 30 but can be specified using the -t argument.
# All FASTQ files in the input directory will be processed in parallel.
# If an output path is provided, the output file will be generated separately for each input file.
# Use the -c option to set the number of CPU cores to use for analysis (default: max(CPU - 1)).

# Print output to stdout if the output directory path isn't provided by using the -o option.
dfci-analysis count_reads_length_glt_threshold -p <path_to_fastq_directory>

# Generate output file(s)
dfci-analysis -o <path_to_output_directory> count_reads_length_glt_threshold -p <path_to_fastq_directory>
```

### Q2. 
**Given a FASTA file with DNA sequences, find 10 most frequent sequences and return the 
sequence and their counts in the file.**
  
- Source code:
The module used for finding the top N most frequent sequences can be found in `/bin/fasta_analyzer.py`.
```bash
# The default the top N is set to 10 but can be specified using the -n argument.
# Input path can be the directory of FASTA files or a single FASTA file.
# All FASTA files in the input directory will be processed in parallel.
# If an output path is provided, the output file will be generated separately for each input file.
# Use the -c option to set the number of CPU cores to use for analysis (default: max(CPU - 1)).

# Print output to stdout if the output directory path isn't provided by using the -o option.
dfci-analysis find_most_frequent_sequences -p <path_to_fasta_directory>

# Generate output file(s)
dfci-analysis -o <path_to_output_directory>  find_most_frequent_sequences -p <path_to_fasta_directory>
```

### Q3. 
**Given a chromosome and coordinates, write a program for looking up its annotation. Output Annotated file of gene name that 
input position overlaps.**

- Source code:  
The module used for annotating genomic positions can be found in `/bin/annotator.py`.
```bash
# The input file will be split into smaller chunks and annotated in parallel.
# Use the -c option to set the number of CPU cores to use for analysis (default: max(CPU - 1)).
# Generate output file to work directory if the output directory path isn't provided by using the -o option.
dfci-analysis -o <path_to_output_directory> annotation -p <path_to_coordinates_to_annotate (TXT file)> -k <path_to_knowledgebase (GTF file)>
```

### Q4.
**Parse the given Example.hs_intervals.txt file. The file contains information on covereage on 
exon level in a hybrid capture panel. The file is a tab-delimited text file. Report the mean 
target coverage for the intervals grouped by GC% bins. Bin in 10%GC intervals (e.g. >= 0 to < 
10; >= 10 to < 20; etc). Note that in the file, GC values range from 0 to 1 rather than 
percentage.**
  
- Note:
Relevant Columns in the file: %gc and mean_target_coverage

- **Hints**:
1. There are certain rows in the input file `Example.hs_intervals.txt` that contain additional tabs, which can cause parsing issues. These tabs will be temporarily removed during execution.
2. Rows have additional tabs: `line # 1277, line # 2558, line # 3833`.
  
- Source code:
The module used for computing mean coverage can be found in `/bin/interval_analyzer.py`.
```bash
# The default number of bins is set to 10 but can be specified using the -b argument.
# The default column used to calculate the group mean is set to '%gc' but can be specified using the -g argument.
# Specifying a non-numeric column using the -g argument will cause an issue.
# Generate output file to work directory if the output directory path isn't provided by using the -o option.
dfci-analysis -o <path_to_output_directory> interval_summary -p <path_to_Example.hs_intervals.txt>
```

### Q5.
**Given a list of variant IDs, using Ensembl API retrieve information about alleles, locations, 
effects of variants in transcripts, and genes containing the transcripts.**

- Note:
1. You can use either PERL API or REST API.
2.  Example variant id: rs56116432

- Source code:
The module used for calling REST APIs can be found in `bin/ensembl_variant_info_retriever.py`
  
- Test input:
`tests/variant_ids/variant_ids.tsv`
  
- **Warning**:
The entire REST API is subject to rate limiting, with a limit of 55,000 requests per hour. Exceeding this limit will result in denial of further requests.
  
```bash
# Generate output file to work directory if the output directory path isn't provided by using the -o option.
dfci-analysis -o <path_to_output_directory> variant_info -p <path_to_variant_ids>
```
****
## Essay questions

### Q6.
**How would you architect a framework for sharing large files (10Gb-25Gb) on the 
cloud with access controls at the file level? We want to share the same file with 
multiple users without making a copy. The users should be able to have access to 
the data on any cloud platform to run bioinformatics analysis pipelines. The 
users can run any cloud service, there is no restriction. The framework’s 
responsibility is only to make data accessible with access controls.**

A: A highly effective solution for sharing large files is leveraging cloud object storage services like AWS S3. This platform boasts exceptional accessibility. Moreover, AWS S3 provides robust access control mechanisms to ensure secure file management.

Access control within AWS S3 is facilitated through a multi-layered approach, primarily categorized into resource-based policies and user policies. Resource-based policies, including bucket policies and access point policies, are directly linked to resources like buckets and objects. For instance, access point policies enable us to restrict usage to specific VPC endpoints or whitelist IP addresses, enhancing security measures.

On the other hand, user policies are attached to individual users within your AWS account. By granting fundamental read and write permissions to other AWS users, we ensure controlled access. Preferably, it's best practice to grant permissions to IAM roles rather than individual users. IAM roles function similarly to IAM users but are intended to be assumed by anyone who requires them, even across accounts. This approach enhances cross-account AWS S3 access control.

For users without AWS accounts, a viable solution is to provide a presigned URL, granting temporary access to download the file securely. This method ensures seamless file sharing while maintaining security measures.

Based on the aforementioned access control measures, all authorized users, whether within our AWS account or granted access via IAM roles or presigned URLs, can securely download files from our bucket.

### Q7.
**Evaluate the benefits and limitations of using containerization and container 
orchestration technologies, such as Docker and Kubernetes, for deploying and 
managing bioinformatics HPC workloads in the cloud.**

A: Using containerization offers several advantages such as isolation for applications and portability across environments. Docker images facilitate easy version control and encapsulate the entire software stack for simplified management. Additionally, containers are more resource-efficient than VMs as they don't require their own operating system.

Container orchestration technologies like Kubernetes further enhance scalability by dynamically allocating resources for optimal performance, ensuring an elastic environment.

Nevertheless, containerization may lead to increased image storage requirements, and configuring communication between containers can be complex. Moreover, OS kernel sharing poses security risks, expanding the attack surface if not appropriately managed.

### Q8.
**For the following SQL statement, what is wrong with it and how would you fix it**
```sql
-- Question:
SELECT UserId, AVG(Total) AS AvgOrderTotal
FROM Invoices
HAVING COUNT(OrderId) >= 1
```

A: The `HAVING` clause in SQL is specifically designed to filter results based on aggregate functions, and it must come after the `GROUP BY` clause. This means that when using aggregate functions, we first need to specify how to group your data. In this case, we can group the data by the `UserId` column:
```sql
SELECT UserId, AVG(Total) AS AvgOrderTotal
FROM Invoices
GROUP BY UserId
HAVING COUNT(OrderId) >= 1;
```

## Contributor
- Yifei Wan (wanyifei0123@gmail.com)
