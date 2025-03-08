#!/bin/bash -l

#SBATCH --job-name=deeptools_correlation
#SBATCH --cpus-per-task=10
#SBATCH --mem=20GB
#SBATCH --time=2:00:00
#SBATCH -o deeptools_correlation.%j.o
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pande250@umn.edu

# Activate Conda
source ~/miniconda3/etc/profile.d/conda.sh
conda activate chipseq

# Change to working directory
cd /scratch.global/pande250/GCD-8141/Project-2/bigWig

multiBigwigSummary bins -b *.bw --binSize 500 -o correlation_matrix.npz --outRawCounts raw_counts.txt && \

plotCorrelation -in correlation_matrix.npz --corMethod spearman --whatToPlot heatmap --removeOutliers --plotTitle "ChIP-seq Sample Correlation" -o correlation_heatmap.png && \

plotCorrelation -in correlation_matrix.npz --corMethod spearman --whatToPlot scatterplot --plotTitle "ChIP-seq Sample Correlation" -o correlation_scatterplot.png


echo "DeepTools correlation analysis completed!"


