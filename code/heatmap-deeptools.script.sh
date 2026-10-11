#!/bin/bash -l

#SBATCH --job-name=chip_heatmap
#SBATCH --partition=amdsmall
#SBATCH --cpus-per-task=2
#SBATCH -n 1
#SBATCH --mem=60gb
#SBATCH -t 8:00:0
#SBATCH -o heatmap.%j.o
#SBATCH --mail-user=your_emailid

# Set this to your project directory
wkDir="${PROJECT_DIR:-/path/to/your/project}"

# Load necessary modules
source ~/miniconda3/etc/profile.d/conda.sh
conda activate chipseq

# Change to working directory
cd $wkDir/heatmap/heatmap-2

# Define input and output files
BW_FILES="Control_shRNA.bw GATA3_shRNA.bw GATA3_shRNA_JUN_OE.bw JUN_OE_DOX.bw"
PEAK_FILE="Control_shRNA_summits.bed GATA3_shRNA_summits.bed GATA3_shRNA_JUN_OE_summits.bed JUN_OE_DOX_summits.bed"

# Step 1: Compute matrix for heatmap
computeMatrix reference-point \
    --referencePoint center \
    -S ${BW_FILES} \
    -R ${PEAK_FILE} \
    -b 3000 -a 3000 \
    --skipZeros \
    -p 2 \
    -o peaks_matrix.gz
# Step 2: Generate heatmap
plotHeatmap \
    -m peaks_matrix.gz \
    -out chipseq_heatmap.pdf \
    --dpi 300 \
    --colorMap Reds \
    --whatToShow 'heatmap and colorbar'


