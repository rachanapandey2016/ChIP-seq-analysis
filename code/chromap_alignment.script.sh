#!/bin/bash -l

#SBATCH --job-name=chip
#SBATCH --cpus-per-task=2
#SBATCH -n 1
#SBATCH --mem=60gb
#SBATCH -t 8:00:0
#SBATCH -o chip.%j.o
#SBATCH --mail-user=pande250@umn.edu

### Define directories
wkDir="/scratch.global/pande250/GCD-8141/Project-2"
fastqDir="$wkDir/input_fastq"
bamDir="$wkDir/chromap_output_bed"
indexDir="$wkDir/indexed_human_genome"

# Reference Genome Path
ref="/scratch.global/pande250/GCD-8141/Project-2/human_genome/GRCh38.p14.genome.fa"

# Create necessary directories
mkdir -p "$bamDir"
mkdir -p "$indexDir"

# Activate Conda and Ensure Chromap is Available
source ~/miniconda3/etc/profile.d/conda.sh
conda activate chipseq

# Debugging: Check if Chromap is accessible in Slurm
echo "Checking Chromap version inside Slurm job:"
which chromap
chromap --version

# Ensure Chromap index exists
cd "$indexDir"
if [ ! -f "GRCh38_chromap_index" ]; then

    echo "Creating Chromap index for the full genome..."

    chromap -i -r "$ref" -o GRCh38_chromap_index

else

    echo "Chromap index already exists."

fi

# Perform alignments with Chromap (FULL GENOME)
cd "$fastqDir"
for f in *.fastq; do
    samp=$(basename "$f" .fastq)

    echo "Processing file: $f"

    # Run Chromap alignment
    chromap --preset chip \
        -x "$indexDir/GRCh38_chromap_index" \
        -r "$ref" \
        -q 20 \
        --min-read-length 15 \
        -1 "$fastqDir/$f" \
        -o "$bamDir/$samp.chromap.bed"

    # Verify if the output file was created
    if [ -f "$bamDir/$samp.chromap.bed" ]; then
        echo "Chromap alignment successful for $samp. Converting to UNIX format..."
        dos2unix "$bamDir/$samp.chromap.bed"
    else
        echo "ERROR: Chromap failed for $samp. No output file found!"
    fi
done
 
