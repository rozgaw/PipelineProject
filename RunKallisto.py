# file for running kallisto
# impoting os for using os.system
import os
# importing sys and argparse
import sys
import argparse

# function to parse command line arguments
def check_arg(args=None):
    parser = argparse.ArgumentParser(description='Wrapper python script for running the differential expression pipeline')
    return parser.parse_args(args)

# retrieve command line arguments
arguments = check_arg(sys.argv[1:])

# assigns parent_folder to be the folder with the data (fastq files)
parent_folder = "data"

# iterates over each sample folder within data
for sample_folder in os.listdir(parent_folder):
    # get the sample names from the folder name
    sample_name = os.path.basename(sample_folder)

    # gets the paths to the paired-end FASTQ files
    fastq1 = os.path.join(sample_folder, f"{sample_name}_1.fastq")
    fastq2 = os.path.join(sample_folder, f"{sample_name}_2.fastq")
    
    # runs kallisto quant for the current sample
    os.system(f"kallisto quant -i index.idx -o results_kal/{sample_name} -b 10 -t 2 {fastq1} {fastq2}")