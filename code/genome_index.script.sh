#!/bin/bash
#
#SBATCH --time=8:00:00
#SBATCH --ntasks=20
#SBATCH --mem=50g
#SBATCH --tmp=120g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=uour_emailid 


# Set this to your project directory (or export PROJECT_DIR)
wkDir="${PROJECT_DIR:-/path/to/your/project}"


#Indexing the reference genome
chromap -i -r $wkDir/human_genome/GRCh38.p14.genome.fa.gz -o $wkDir/indexed_human_genome/GRCh38_chromap_index

