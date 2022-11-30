#!/bin/bash
# Slurm options
#SBATCH --mail-user=*your e-mail address*
#SBATCH --mail-type=begin,fail,end
#SBATCH --job-name="fastqc"
#SBATCH --workdir=/path to your folder/
#SBATCH --time=06:00:00 
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=16G

# This code uses the following software/algorithms: 
# Trimmomatic: Bolger, A.M., Lohse, M., and Usadel, B. (2014). Trimmomatic: a flexible trimmer for Illumina sequence data. Bioinformatics 30, 2114–2120. 10.1093/bioinformatics/btu170.

# Put your code below this line

module load vital-it
module load UHTS/Quality_control/fastqc/0.11.7
module load UHTS/Analysis/trimmomatic/0.36
module load UHTS/Analysis/MultiQC/1.8

#Line 3
#controls (homozygous axN): 11, 12, 15, 17
#CRISPR mutants (homozygous myb-fl mutants): 6, 9, 20, 24

#Line 10
#controls (homozygous axN): 3, 6, 8, 12
#CRISPR mutants (homozygous myb-fl mutants): 2, 9, 10, 16

### running trimmomatic to remove any remaining adapters, and also quality trimming

trimmomatic SE -threads 4 fastq.gz_file_name ../path_to_file/ ILLUMINACLIP:indexadaptersfile:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36


#fastqc to see if this improved quality and removed all adapters
fastqc ../path_to_file/ -t 4


#multiqc to make a unified report of the fastqc runs
multiqc ../path_to_folder/
