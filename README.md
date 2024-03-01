# PipelineProject
Differential expression pipeline to compare HCMV transcriptomes of two individuals 2- and 6-days post-infection.

Before running this pipeline, dowload the transcriptomes from two patient donors from SRA and convert to paired-end fastq files. This can be done by running the following commend in your terminal (where SRA_File_Links is a txt file that contains each of the SRA links for the patients on a new line:

wget -i SRA_File_Links

For this example, SRA_File_Links.txt contains the following:

