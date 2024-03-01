# PipelineProject
Differential expression pipeline to compare HCMV transcriptomes of two individuals 2- and 6-days post-infection (2dpi, 6dpi).

### Requirements:
- Python: BioPython SeqIO, BioPython Entrez; sys; argparse; pandas
- R: sleuth, dplyr
- local machine/server: kallisto, blast, sra-tools 

### Before Running Pipeline
Before running this pipeline, download the transcriptomes from two patient donors from SRA and convert to paired-end fastq files. This can be done by running the following command in your terminal (where SRA_File_Links is a text file that contains each of the SRA links for the patients on a new line.

Download the SRA files from terminal:

wget -i SRA_File_Links  

For this example, SRA_File_Links.txt contains the following with a newline after each url (see SRA_File_Links for correct formatting):  

https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660030/SRR5660030  
https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660033/SRR5660033  
https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660044/SRR5660044  
https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660045/SRR5660045  

Next, convert the downloaded files into paired-end fastq files. This can be done by running the following code in the terminal, calling each SRA file by name.  

for file in SRR5660030 SRR5660033 SRR5660044 SRR5660045; do   
    fastq-dump -I --split-files "$file"   
done   

For different analyses, make sure to update the sample names to match the downloaded files.   

Next, rearrange the file organization. I downloaded the SRA files into a folder called "data" and then grouped each pair of paired-end reads into a folder of their sample name within the "data" folder.   
For example, here are the file paths for the paired-end reads from sample SRR5660030 where ~ signifies the working directory:
   
~/data/SRR5660030/SRR5660030_1.fastq    
~/data/SRR5660030/SRR5660030_2.fastq   

### Running the Pipeline
Finally, the wrapper script in this repository will do the following:
- Build transcriptome index for HCMV (NCBI accession NC_006273.2) using BioPython and Kallisto and return the number of CDS regions in the HCMV genome.
- Quantify the TPM of each CDS in each transcriptome using Kallisto returning the sample, condition, min_tpm, med_tpm, mean_tpm, and max_tpm for each sample.
- Use Sleuth to find differentially expressed genes between the two conditions returning the target_id, test_stat, pval, and qval of significant transcripts (FDR < 0.05)
- Retrieve a protein fasta file of the most differentially expressed CDS in the reference genome, running a blast of it against the Betaherpesvirinae subfamily, returning the subject accession, percent identity, alignment length, start of alignment in query, end of alignment in query, start of alignment in subject, end of alignment in subject, bit score, e-value, and subject title for the top ten hits.

For this pipeline, sample data is provided. This data is a fragment of the original data modified by running the code below in the terminal. Update data and sampledata to the correct input file name and output.     
    
head -n 40000 data.fastq > sampledata.fastq     
The sample data is stored in the data folder.    

### Notes
Notes:
- Make sure to set your working directory with a similar file path organization for the wrapper to run correctly.
- PipelineProject folder has data folder and all the required .py and r scripts.
- Data folder has the paired end fastq files for each sample in a folder named the sample name.
- Make sure to update the sample dictionary (samp_dict)in the wrapper python script to reflect your samples, conditions, and file paths to the kallisto output.

### Running wrapper.py
Finally to run the wrapper python script, run the following command in your terminal:     
"python Wrapper.py -email <your_email>" replacing the "<your_email>" with your actual email adress
