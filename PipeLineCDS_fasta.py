# importing sys and argparse packages
import sys
import argparse

# function to parse command line arguments
def check_arg(args=None):
    parser = argparse.ArgumentParser(description='takes text file input with NCBI genbank id, returns fasta file of all CDS regions')
    parser.add_argument("-i", "--input", help = "input file", required = True)
    parser.add_argument("-o", "--output", help = "output file", required = True)
    return parser.parse_args(args)

# retrieve command line arguments
arguments = check_arg(sys.argv[1:])
infile = arguments.input
outfile = arguments.output

# infile will contain the GenBank ID for the genome of interest and email on newline
# opens infile for reading
with open(infile, "r") as f:
    # reads in the text file, assigning the first line to genbank_id
    genbank_id = f.readline().strip()
    # reads in the second line and assigns to email
    email = f.readline()

# imports Entrez and SeqIO for reading in and parsing genome
from Bio import Entrez
from Bio import SeqIO

# using Entrez.efetch to retrieve genbank_id entry and read with SeqIO
Entrez.email = email
handle = Entrez.efetch(db = "nucleotide", id = genbank_id, rettype = "gb", retmode="text")
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

with open(outfile, 'w') as o: # opens the outfile for writing
    # writing the output string of fasta_out to the outfile, removing whitespace from end
    o.write(fasta_out.strip()) # automatically closes file after done writing