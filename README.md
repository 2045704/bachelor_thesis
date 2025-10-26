# Internship Bioinformatics Scripts

This repository contains the main scripts I developed during my internship in my third year.

## Data Sources

- **Patient data:** [GSE124814](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE124814)  
- **Cell line data:** [GSE171117](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE171117)

## Description of Scripts

### Python Scripts

**R4**  
This script defines a class that sets up and runs the command-line prompts required for several bioinformatic tools.

**PROCESS1**  
This script calls the functions of the class defined in *R4* to process RNA-sequencing data.  
The pipeline includes:  
- Downloading samples using the SRA Toolkit  
- Performing quality control on downloaded samples with FastQC  
- Trimming FASTQ files using Cutadapt  
- Running a second round of quality control with FastQC  
- Aligning reads with Salmon  
- Removing downloaded FASTQ files to save disk space  

**Patient_Analysis**  
A Jupyter Notebook in which the patient-derived data were analyzed and compared with the sequencing data from the cell line.


### R Scripts

**BM**  
Biomart script used to obtain gene symbols.

**EdgeR**  
Differential expression analysis pipeline performed after obtaining the aligned reads from the cell line data.
