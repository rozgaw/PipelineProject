'''
Script for retrieving genbank genome for Human cytomegalovirus (HCMV)
and returning fasta file which includes only the CDS regions from the 
HCMV genome in a file called HCMV_CDS.fasta
'''
import argparse

# function to parse command line arguments
def check_arg(args=None):
    parser = argparse.ArgumentParser(description='Wrapper python script for running the differential expression pipeline')
    parser.add_argument('-email', help='input your email', required=True)
    return parser.parse_args(args)

# Retrieve command line arguments
arguments = check_arg()

# assigns the email argument from argparse to email variable
email = arguments.email

# imports Entrez and SeqIO for reading in and parsing genome
from Bio import Entrez
from Bio import SeqIO

# using Entrez.efetch to retrieve genbank_id entry and read with SeqIO
Entrez.email = email
handle = Entrez.efetch(db = "nucleotide", id = "NC_006273.2", rettype = "gb", retmode="text")
# assigns record to SeqRecord of the genome of interest
record = SeqIO.read(handle, "genbank")

# initiates empty string for storing CDS regions in fasta format
fasta_out = ''
# iterates through the features in SeqRecord
for feature in record.features:
    # if the feature is a CDS
    if feature.type == "CDS":
        # assigns cds_seq to the seq of that record
        cds_seq = feature.extract(record.seq)
        # assigns protein_id the protein_id or Unkown if missing an id
        protein_id = feature.qualifiers.get("protein_id", ["Unknown"])[0]
        # write a > andt the protein id, followed by a new line, cds_seq, and newline to fasta_out string
        fasta_out += f">{protein_id}\n{cds_seq}\n"

outfile = "HCMV_CDS.fasta"

with open(outfile, 'w') as o: # opens the outfile for writing
    # writing the output string of fasta_out to the outfile, removing whitespace from end
    o.write(fasta_out.strip()) # automatically closes file after done writing