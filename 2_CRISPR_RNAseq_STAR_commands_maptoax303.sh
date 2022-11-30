#!/bin/bash
# Slurm options
#SBATCH --mail-user=*your e-mail address*
#SBATCH --mail-type=begin,fail,end
#SBATCH --job-name="STAR mapping for CRISPR MYB-FL mutants"
#SBATCH --workdir=/path to your folder
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=12500M
#SBATCH --output=STAR_%A.out
#SBATCH --error=STAR_%A.err

# This code uses the following software/algorithms: 
# STAR: Dobin, A., Davis, C.A., Schlesinger, F., Drenkow, J., Zaleski, C., Jha, S., Batut, P., Chaisson, M., and Gingeras, T.R. (2013). STAR: ultrafast universal RNA-seq aligner. Bioinformatics 29, 15–21. 10.1093/bioinformatics/bts635.
# and is based on http://homer.ucsd.edu/homer/basicTutorial/mapping.html

# Put your code below this line

module load vital-it/7 
module load UHTS/Aligner/STAR/2.6.0c 

#align RNA-Seq data with STAR to a reference genome
#step 1: build a genome index

STAR  --runMode genomeGenerate --runThreadN 16 --genomeDir /path_to_P_axillaris_genome_folder/ --genomeFastaFiles /path_to_genome_file/ --sjdbGTFfile /path to gff file/ --sjdbGTFtagExonParentTranscript Parent --genomeChrBinNbits 16 --limitGenomeGenerateRAM 200000000000
# 16 CPUs detailed above, so 16 threads possible here

#Running mapping jobs
#need --genomeLoad LoadAndKeep 
#using 2pass method for downstream SNP analysis, and outputting as coord-sorted bams to save some steps, default usually is 1pass, 2pass goes and remaps again, more accurate

STAR --genomeDir /path_to_P_axillaris_genome/ --runThreadN 16 --readFilesIn /path_to_fastq.gz_file/ --sjdbOverhang 100 --outFilterType BySJout --outFilterMultimapNmax 20 --outSAMtype BAM SortedByCoordinate --twopassMode Basic --outFileNamePrefix axN_3_11_1_ --readFilesCommand zcat --genomeLoad NoSharedMemory


# output BAM files 
