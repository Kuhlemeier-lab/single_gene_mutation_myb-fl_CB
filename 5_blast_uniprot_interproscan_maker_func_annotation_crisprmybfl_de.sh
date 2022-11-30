#!/bin/bash
#SBATCH --mail-user=*your e-mail address*
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=32GB
#SBATCH --workdir=/path to your folder/
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -J "blastp and interproscan manual annotation"

# This code uses the following software/algorithms:
# Maker v2.31.9: University of Utah, Department of Human Genetics

module load vital-it/7
module load Blast/ncbi-blast/latest
module load SequenceAnalysis/SequenceAlignment/interproscan/5.33.72.0
module load SequenceAnalysis/GenePrediction/maker/2.31.9
module load UHTS/Analysis/seqtk/1.2


#make list of genes you want re-annotated, copy from excel and make into list crispr_de_list.txt
#you need two versions of this file - the first needs to have the ending .1.mrna1 and the second needs to have .1.path1 to match different files for subsetting
#the pep file has genes listed as Peaxi162Scf00000g00013.1.mrna1 and the gff has genes listed as both but use .path1
	#cp crispr_de_list.txt crispr_de_list_for_grep.txt
	#sed -i 's/path/mrna/g' crispr_de_list.txt #find and replace the word path with mrna
	sed -i 's/.1.path1//g' crispr_de_list_for_grep.txt #find and replace so that only the gene name itself is present, to get all parts of the gff

#get protein sequences from list - use pep1.fasta because asterisks at the ends of the sequences have been removed
seqtk subseq ~/path_to_protein_sequence_file/ crispr_de_list.txt > Pax_crisprmybfl_de_pep.fasta

### blastp functional annotation

blastp -query Pax_crisprmybfl_de_pep.fasta -db ~/path_to_uniprot_fasta_file/ -evalue 1e-6 -max_hsps 1 -max_target_seqs 1 -outfmt 6 -out file_name_output_file_blastp


### interproscan functional annotation

interproscan.sh -appl pfam -dp -f TSV -goterms -iprlookup -pa -t p -i Pax_crisprmybfl_de_pep.fasta -o file_name_output_file_interproscan


#Rename gene models and other data.
bedtools sort -i crisprmybfl_de_genes.gff > crisprmybfl_de_genes.sorted.gff

maker_map_ids --prefix Peaxi --justify 8 crisprmybfl_de_genes.sorted.gff > map
map_gff_ids map crisprmybfl_de_genes.sorted.gff
map_fasta_ids map Pax_crisprmybfl_de_pep.fasta
map_data_ids map Pax_crisprmybfl_de_pep.fasta.blastp
map_data_ids map Pax_crisprmybfl_de_pep.fasta.iprscan


###integrate results into gff
ipr_update_gff crisprmybfl_de_genes.sorted.gff  Pax_crisprmybfl_de_pep.fasta.iprscan > crisprmybfl_de_genes.sorted.all.gff
iprscan2gff3 Pax_crisprmybfl_de_pep.fasta.iprscan crisprmybfl_de_genes.sorted.gff >> crisprmybfl_de_genes.sorted.all.gff 

maker_functional_gff /$PATH/uniprot/uniprot_sprot.fasta Pax_crisprmybfl_de_pep.fasta.blastp crisprmybfl_de_genes.sorted.gff > crisprmybfl_de_genes.sorted.all.gff

