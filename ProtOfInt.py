# importing sys and argparse packages
import sys
import argparse

# function to parse command line arguments
def check_arg(args=None):
    parser = argparse.ArgumentParser(description='takes proteinID and email, returns fasta file of that protein')
    parser.add_argument("-e", "--email", help = "email", required = True)
    parser.add_argument("-p", "--proteinID", help = "NCBI accession for protein of interest", required = True)
    return parser.parse_args(args)

# retrieves command line arguments
arguments = check_arg()

# assigns arguments to correct variables
email = arguments.email
ProtID = arguments.proteinID

# imports Entrez and SeqIO for reading in and parsing protein
from Bio import Entrez

# using Entrez.efetch to retrieve genbank_id entry and read with SeqIO
Entrez.email = email
# makes handle for retrieving the protein fasta file
handle = Entrez.efetch(db="protein", id=ProtID, rettype="fasta", retmode="text")
# saves the protein sequence to a FASTA file
with open(f"MDEprotein.fasta", "w") as fasta_file:
    fasta_file.write(handle.read())