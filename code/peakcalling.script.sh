#!/bin/bash -l

#SBATCH --job-name=peakcalling
#SBATCH --cpus-per-task=2
#SBATCH -n 1
#SBATCH --mem=60gb
#SBATCH -t 8:00:0
#SBATCH -o peakcalling.%j.o
#SBATCH --mail-user=pande250@umn.edu

### Define directories
wkDir="/scratch.global/pande250/GCD-8141/Project-2"
fastqDir="$wkDir/input_fastq"
bamDir="$wkDir/chromap_bed_into_bam" 
indexDir="$wkDir/indexed_human_genome"
peakDir="$wkDir/peak_calling"
ref="/scratch.global/pande250/GCD-8141/Project-2/human_genome/GRCh38.p14.genome.fa"

##get access to needed software 
conda activate chipseq
module load fastqc/0.12.1
module load samtools/1.16.1-gcc-8.2.0-egljrr3
module load macs/2.1.1
module load bedtools2/2.31.0-gcc-8.2.0-7j35k74

# Define Control and Treatment Samples
control="$bamDir/Control_shRNA.sorted.bam"
treatments=(
    "GATA3_shRNA"
    "JUN_oe_dox"
    "GATA3_shRNA_JUN_oe"
)

# Remove old FRiP score file
rm -f $peakDir/FRiP_scores.tsv

# Loop through each treatment sample
for sample in "${treatments[@]}"; do
    macs2 callpeak \
        -t "$bamDir/${sample}.sorted.bam" \
        -c "$control" \
        -f BAM \
        -g hs \
        --nomodel --shift -100 \
        --extsize 200 \
        -n "$sample" \
        --outdir "$peakDir"

    # Compute FRiP Score
    readsInSample=$(samtools view -c "$bamDir/${sample}.sorted.bam")
    readsInPeaks=$(bedtools intersect -u -a "$bamDir/${sample}.sorted.bam" -b "$peakDir/${sample}_peaks.narrowPeak" -ubam | samtools view -c)
    frip_score=$(echo "scale=6; ${readsInPeaks} / ${readsInSample}" | bc)

    # Save FRiP score to file
    echo -e "${sample}\t${frip_score}" >> $peakDir/FRiP_scores.tsv

    echo "Finished processing $sample."
done

echo "All samples processed successfully!"
