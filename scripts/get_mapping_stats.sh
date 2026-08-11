#!/bin/bash
# This script is used to get the mapping statistics of the input bam file (from stdout of bbmap).
# Usage: bash get_mapping_stats.sh <input_stdout_text_file>
# Example: bash get_mapping_stats.sh bbmap_stdout.txt


# Check if the input file is provided
if [ -z $1 ]; then
    echo "Please provide the input file."
    exit 1
fi

# Get the input file
input_file=$1

# print the input file
echo "Input file: $input_file"

# Get the mapping statistics using regular expressions. get: sample names from e.g.:
# "Starting mapping for sample: ID_007579_Fnv_PNA79_16_h_C.fq.gz"
# total input reads from:
# "Reads Used:           	8349297	(726179831 bases)"
# total mapped reads from (extract third column):
# mapped:          	 86.7418% 	 11112135 	 91.9811% 	  1091516226
grep "Starting mapping for sample:" $input_file | sed 's/Starting mapping for sample: //' > sample_names.txt
grep "Reads Used:" $input_file | sed 's/Reads Used:           	//' | awk '{print $1}' > total_input_reads.txt
grep "mapped:" $input_file | awk '{print $3}' > total_mapped_reads.txt

# Combine the statistics into a single file. add header to the file
echo -e "Sample\tTotal_input_reads\tTotal_mapped_reads" > mapping_stats.txt
paste sample_names.txt total_input_reads.txt total_mapped_reads.txt >> mapping_stats.txt

# Remove the intermediate files
rm sample_names.txt total_input_reads.txt total_mapped_reads.txt

# Print the mapping statistics
echo "Mapping statistics are saved in mapping_stats.txt file."

