# PipelineProject
Differential expression pipeline to compare HCMV transcriptomes of two individuals 2- and 6-days post-infection.

Before running this pipeline, download the transcriptomes from two patient donors from SRA and convert to paired-end fastq files. This can be done by running the following command in your terminal (where SRA_File_Links is a text file that contains each of the SRA links for the patients on a new line.

Download the SRA files from terminal:

wget -i SRA_File_Links

For this example, SRA_File_Links.txt contains the following:
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

Finally, the wrapper script in this repository will do the following:
- Build transcriptome index for HCMV (NCBI accession NC_006273.2) using BioPython and Kallisto and return the number of CDS regions in the HCMV genome.
- Quantify the TPM of each CDS in each transcriptome using Kallisto returning the sample, condition, min_tpm, med_tpm, mean_tpm, and max_tpm for each sample.
- Use Sleuth to find differentially expressed genes between the two conditions returning the target_id, test_stat, pval, and qval of significant transcripts (FDR < 0.05)
- Retrieve a protein fasta file of the most differentially expressed CDS in the reference genome, running a blast of it against the Betaherpesvirinae subfamily, returning the subject accession, percent identity, alignment length, start of alignment in query, end of alignment in query, start of alignment in subject, end of alignment in subject, bit score, e-value, and subject title for the top ten hits.
