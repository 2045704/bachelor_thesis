This small repository contains the main codes developed by me during my internship in my third year. 
Data files were obtained from:
Patient data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE124814
Cell line data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE171117
The description of the scripts goes as follows:
*PYTHON SCRIPTS*
- R4: This script contains the class which sets up and runs the command line prompts required for usage of some of the bioinformatic tools
- PROCESS1: This script calls on the functions of the class contained in R4 to process the RNA-Sequencing data, processing as follows
   - Download of the samples with SRA toolkit
   - Quality control of downloaded samples with FastQC
   - Trimming of fastq files with Cutadapt tool
   - Quality control after trimming wit FastQC
   - Alignment with Salmon
   - Removal of downloaded fastq file for space-sake
- Patient_Analysis: Jupyter notebook in which the patient-derived data was analysed, then crossed with the sequencing data from the cell line
*R SCRIPTS*
- BM: Biomart script to obtain gene symbols
- EdgeR: Differential expression pipeline ran after obtaining the aligned reads from the cell line data.
