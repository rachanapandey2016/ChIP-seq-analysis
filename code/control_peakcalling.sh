#!/bin/bash -l

#SBATCH --job-name=peakcallingcontrol
#SBATCH --cpus-per-task=2
#SBATCH -n 1
#SBATCH --mem=60gb
#SBATCH -t 2:00:0
#SBATCH -o controlpeakcalling.%j.o
#SBATCH --mail-user=pande250@umn.edu

### Define directories
wkDir="/scratch.global/pande250/GCD-8141/Project-2"
fastqDir="$wkDir/input_fastq"
bamDir="$wkDir/chromap_bed_into_bam"
indexDir="$wkDir/indexed_human_genome"
peakDir="$wkDir/peak_calling"
ref="/scratch.global/pande250/GCD-8141/Project-2/human_genome/GRCh38.p14.genome.fa"

## Load necessary software
source ~/miniconda3/etc/profile.d/conda.sh
conda activate chipseq
module load fastqc/0.12.1
module load samtools/1.16.1-gcc-8.2.0-egljrr3
module load macs/2.1.1
module load bedtools2/2.31.0-gcc-8.2.0-7j35k74

# Define control sample
control="$bamDir/Control_shRNA.sorted.bam"
control_peak="$peakDir/Control_shRNA_peaks.narrowPeak"

# Step 1: Call peaks using MACS2
macs2 callpeak \
    -t "$control" \
    -f BAM \
    -g hs \
    --nomodel --shift -100 \
    --extsize 200 \
    -n "Control_shRNA" \
    --outdir "$peakDir"

echo "Control sample peak calling completed."

# Step 2: Compute FRiP Score for Control Sample
readsInControl=$(samtools view -c "$control")
readsInControlPeaks=$(bedtools intersect -u -a "$control" -b "$control_peak" -ubam | samtools view -c)
frip_control_score=$(echo "scale=6; ${readsInControlPeaks} / ${readsInControl}" | bc)

# Save FRiP score to file
echo -e "Control_shRNA\t${frip_control_score}" >> $peakDir/FRiP_scores.tsv

echo "FRiP score for control sample: ${frip_control_score}"
echo "Control sample processing completed successfully!"


