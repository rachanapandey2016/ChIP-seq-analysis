#!/bin/bash
#
#SBATCH --time=96:00:00
#SBATCH --ntasks=20
#SBATCH --mem=50g
#SBATCH --tmp=120g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pande250@umn.edu 


#Indexing the reference genome
chromap -i -r /scratch.global/pande250/GCD-8141/Project-2/human_genome/GRCh38.p14.genome.fa.gz -o /scratch.global/pande250/GCD-8141/Project-2/indexed_human_genome/GRCh38_chromap_index
