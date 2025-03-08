#!/bin/bash -l

#SBATCH --job-name=chip
#SBATCH --cpus-per-task=2
#SBATCH -n 1
#SBATCH --mem=60gb
#SBATCH -t 8:00:0
#SBATCH -o bedtobam.%j.o
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pande250@umn.edu

# Define necessary directories
wkDir="/scratch.global/pande250/GCD-8141/Project-2"
bedDir="$wkDir/chromap_output_bed"  # Contains BED files
bamDir="$wkDir/chromap_bed_into_bam"  # Output directory for BAM files
bigWigDir="$wkDir/bigWig"           # Output directory for BigWig files
indexDir="$wkDir/samtools_indx_genome"

# Reference Genome Path
ref="/scratch.global/pande250/GCD-8141/Project-2/human_genome/GRCh38.p14.genome.fa"

# Create output directories if not present
mkdir -p "$bamDir"
mkdir -p "$bigWigDir"
mkdir -p "$indexDir"  # Create directory for samtools-indexed genome


# Activate Conda environment for bedtools, samtools, and deepTools
source ~/miniconda3/etc/profile.d/conda.sh
conda activate chipseq

#Loading  the necessary module
module load bedtools2/2.31.0-gcc-8.2.0-7j35k74
module load samtools/1.16.1-gcc-8.2.0-egljrr3

# Navigate to reference genome directory
cd "$indexDir"

# Index the genome using samtools (This will generate a .fai file), but before that our reference gneome must be gunzipped/decompressed if it is in .gz format
echo "Indexing the reference genome with samtools..."
samtools faidx "$ref" --fai-idx GRCh38_samtools_index.fai
cut -f1,2 GRCh38_samtools_index.fai > GRCh38sizes.genome
echo "Genome indexing completed. Genome sizes file (GRCh38sizes.genome) created successfully."

# Convert Chromap BED to BAM, then sort, index, and generate BigWig
cd "$bedDir"

for bed in *.chromap.bed; do
    samp="${bed/.chromap.bed/}"
    echo "Processing sample: $samp"

    # Convert BED to BAM using the newly indexed genome
    bedtools bedtobam -i "$bedDir/$bed" -g "$indexDir/GRCh38sizes.genome" > "$bamDir/$samp.bam"

    # Sort BAM
    samtools sort -o "$bamDir/$samp.sorted.bam" "$bamDir/$samp.bam"

    # Index BAM
    samtools index "$bamDir/$samp.sorted.bam"

    # Generate BigWig
    bamCoverage -p 4 -b "$bamDir/$samp.sorted.bam" --normalizeUsing RPKM -o "$bigWigDir/$samp.norm.bw"

    echo "Finished processing $samp"
done

echo "All samples processed successfully!"
