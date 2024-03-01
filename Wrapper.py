# importing sys, argparse, os, and pandas packages
import sys
import argparse
import os
import pandas as pd

# input arguments
parser = argparse.ArgumentParser(description='Script for differential expression '
    + 'of donors at dif timepoints post infection.')
parser.add_argument('-email', help='input your email', required=True)
args = parser.parse_args()

# function to parse command line arguments
def check_arg(args=None):
    parser = argparse.ArgumentParser(description='Wrapper python script for running the differential expression pipeline')
    parser.add_argument('-email', help='input your email', required=True)
    return parser.parse_args(args)

# Retrieve command line arguments
arguments = check_arg()

# assigns the email argument from argparse to email variable
email = arguments.email

# create and open log file for storing outputs
# using a for appending, so content not overwritten each time
log_file = open("PipelineProject_Roza_Gawin/PipelineProject.log", "a")

# command for running Genome_to_CDS python code from terminal
GenomeToCDS = "python Genome_to_CDS.py -email {email}"
os.system(GenomeToCDS)

# counting CDS regions in HCMV genome
countCDS = 'grep -c ">" HCMV_CDS.fasta'
# saves the output number as numCDS variable
numCDS = os.popen(countCDS)
# writes to the log file
log_file.write("The HCMV genome (NC_006273.2) has {numCDS} CDS.\n")

# running kallisto
# creating kallisto index
kallistoIndex = "kallisto index -i index.idx HCMV_CDS.fasta"
sys.os(kallistoIndex)
# running kallisto on the samples by calling the RunKallisto python script
kallistoRun = "python RunKallisto.py"
os.system(kallistoRun)

# writes headers for TPM info to log file
log_file.write("sample\tcondition\tmin_tpm\tmed_tpm\tmean_tpm\tmax_tpm\n")

# dictionary of samples which stores the condition and kalisto result file path
samp_dict = {"SRR5660030":["2dpi", "/PipelineProject/results_kal/SRR5660030"], "SRR5660044":["2dpi", "/PipelineProject/results_kal/SRR5660044"], "SRR5660033":["6dpi", "/PipelineProject/results_kal/SRR5660033"], "SRR5660045":["6dpi", "/PipelineProject/results_kal/SRR5660033"]}

# function for finding TPM values for each sample using a pandas df for easier calculations
def TPMcalcs(file_path):
    '''makes pd df from abundance tsv file, returning min, med, mean, and max tpm values'''
    abundanceDF = pd.read_table(file_path, delimiter = '\t')
    minTPM = abundanceDF["tpm"].min()
    medTPM = abundanceDF["tpm"].median()
    meanTPM = abundanceDF["tpm"].mean()
    maxTPM = abundanceDF["tpm"].max()
    return minTPM, medTPM, meanTPM, maxTPM

# iterates throught the sample dictionary, writing the TPM calculations to the log file
# by calling the TPM calcs function on each dictionary key (each sample)
for key, values in samp_dict.items():
    log_file.write(key +"\t")
    log_file.write(values[0] + "\t")
    file_path = values[1]+"abundances.tsv"
    minTPM, medTPM, meanTPM, maxTPM = TPMcalcs(file_path)
    log_file.write(minTPM + "\t" + medTPM + "\t" + meanTPM + "\t" + maxTPM + "\n")

# create and open sleuth input
SluIn = open("SleuthInput.txt", "w+")
SluIn.write("sample\tcondition\tpath\n")
# iterates through sample dictionary to make appropriate input file for runnign sleuth
for key, values in samp_dict.items():
    SluIn.write(key +"\t")
    for value in values:
        SluIn.write(value +"\t")
    SluIn.write('\n')
SluIn.close()

# runs sleuth by calling R script for sleuth
# set path to sleuth rscript
sleuth_path = "/PipelineProject/Sleuth.R"
# run sleuth, writting output to the log file (adding on to the contents, and not overwritting)
os.system(f"Rscript {sleuth_path} >> {log_file}")

# creates new file to store Sleuth output
SluOut = open("SluOut.txt", "a")
# running again to get output in another file
os.system(f"Rscript {sleuth_path} >> {SluOut}")
# find most differencially expressed CDS from reference genome

with open("SluOut.txt", "r") as file:
    # read the first line
    MDEprotein = file.readline().strip()  # remove any leading/trailing whitespace
    # splits the line into a list on tab characters
    MDEproteinList = MDEprotein.split("\t")
# saves the first list item (the protein ID as PooID)
ProID = MDEproteinList[0]

# fetching fasta file of protein of interest
proteinFASTA = "python ProtOfInt.py --email {email} --proteinID {ProID}"
os.system(proteinFASTA)

# downloading virus genome for blast database
downloadDB = "datasets download virus genome taxon Betaherpesvirinae --refseq --include genome"
os.system(downloadDB)
# unzipping downloaded data
unzipDB = "unzip ncbi_dataset.zip"
os.system(unzipDB)
# creating blast data base using the downloaded data
blastDB = "makeblastdb -in ncbi_dataset/data/genomic.fna -out blastDB -title blastDB -dbtype nucl"
os.system(blastDB)

# now actually running blast
blastCMD = 'tblastn -query MDEprotein.fasta -db blastDB -out tblastn_results.tsv -outfmt "6 sacc pident length qstart qend sstart send bitscore evalue stitle"' 
os.system(blastCMD)

# writing the top 10 hits from the blast to the log file using counter to track number of lines
with open("tblastn.tsv") as f:
    count_lines = 0
    for line in f:
        log_file.write(f"{line}\n")
        count_lines += 1
        if count_lines == 10:
            break

# close log file
log_file.close()