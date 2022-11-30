#!/bin/bash
# Slurm options
#SBATCH --mail-user=*your e-mail address*
#SBATCH --mail-type=begin,fail,end
#SBATCH --job-name="featureCounts for CRISPR MYB-FL"
#SBATCH --workdir=/path to your folder
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=12500M

# This code uses the following software/algorithms:
# featureCounts:Liao, Y., Smyth, G.K., and Shi, W. (2014). featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. Bioinformatics 30, 923–930. 10.1093/bioinformatics/btt656.

# Put your code below this line
###16 cpu-per-task (implicitly 1 task, this is default and left out) and 16000M for mem-per-cpu. 16x12500M => 200G.
# this is one node with at least 16 cpu cores and at least 200G of RAM

module load vital-it/7
module load UHTS/Analysis/subread/1.6.0 # featureCounts was downloaded from UBELIX server, University of Bern


### counting genes from "genes" - this includes any intronic reads (generally a difference of less than 10% from our experiments, but closer to HTSeq counts)
#count by gene including exons, introns, 3UTR, and other features nested under gene
# -s = strandedness
# -t has exon as default -> we want to count by gene, not exon
# -g ID determines what the gene ID is -> in annotation file

featureCounts -T 16 -t gene -s 2 -g ID \
	-a /path_to_gff_file/ \
	-o ./counts/Crispr_mybfl_gene_counts_pax303_s2.txt \
	axN_3_11_1_Aligned.sortedByCoord.out.bam \
	axN_3_12_1_Aligned.sortedByCoord.out.bam \
	axN_3_15_3_Aligned.sortedByCoord.out.bam \
	axN_3_17_1_Aligned.sortedByCoord.out.bam \
	mut_3_06_1_Aligned.sortedByCoord.out.bam \
	mut_3_09_1_Aligned.sortedByCoord.out.bam \
	mut_3_20_2_Aligned.sortedByCoord.out.bam \
	mut_3_24_1_Aligned.sortedByCoord.out.bam \
	axN_10_03_1_Aligned.sortedByCoord.out.bam \
	axN_10_06_2_Aligned.sortedByCoord.out.bam \
	axN_10_08_1_Aligned.sortedByCoord.out.bam \
	axN_10_12_1_Aligned.sortedByCoord.out.bam \
	mut_10_02_1_Aligned.sortedByCoord.out.bam \
	mut_10_09_1_Aligned.sortedByCoord.out.bam \
	mut_10_16_1_Aligned.sortedByCoord.out.bam \
	mut_10_20_1_Aligned.sortedByCoord.out.bam 
	
# ./counts tells bash that counts is in the current directory
