**ChIP-seq Analysis**
--- 
**Project Title:  Enhancer Reprogramming in MCF7 Breast Cancer: The Impact of GATA3 Knockdown and JUN Overexpression**  
This repository contains the ChIP-seq analysis pipeline for investigating enhancer reprogramming in MCF7 breast cancer cells under different experimental conditions. 
I performed ChIP-seq analysis as a class project for the Computational Genomics course. In this project, I analyzed publicly available ChIP-seq data to understand enhancer activity changes due to GATA3 knockdown and JUN overexpression in MCF7 breast cancer cells.  
**Reference Paper**  
Enhancer reprogramming driven by high-order assemblies of transcription factors promotes phenotypic plasticity and breast cancer endocrine resistance.(https://doi.org/10.1038/s41556-020-0514-z)  
GEO page where data is deposited-(https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA527779&o=acc_s%3Aa&s=SRR8742725,SRR8742736,SRR8742726,SRR8742716)   
GEO acession number- GSE128460  

Experimental Conditions  
---
I analyzed H3K27ac ChIP-seq data from four conditions:
1. Control_shRNA - Normal MCF7 cells
2. GATA3_shRNA - Knockdown of GATA3
3. JUN_OE_DOX - Overexpression of JUN using doxycycline
4. GATA3_shRNA_JUN_OE_DOX - Combined GATA3 knockdown and JUN overexpression with doxycycline

Objectives  
---
1. To determine the impact of GATA3 knockdown in chromatin accessibility (H3K27ac) in MCF7 cells.
2. To analyze the impact of JUN overexpression on enhancer remodeling and transcriptional regulation.
3. To evaluate whether combining GATA3 knockdown and JUN overexpression leads to a 
synergistic effect on enhancer reprogramming in MCF cells.

Tools and Packages Required  
---
1. miniconda- https://www.anaconda.com/docs/getting-started/miniconda/install
2. Sra-tools- https://github.com/ncbi/sra-tools
3. FastQC- https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
4. MultiQC- https://github.com/MultiQC/MultiQC
6. Samtools- https://www.htslib.org/
7. bedtools- https://bedtools.readthedocs.io/en/latest/
8. MACS2- https://pypi.org/project/MACS2/
9. Deep Tools- https://deeptools.readthedocs.io/en/develop/
10. R packages- ggplot2

Execution Modes  
---
This analysis uses both interactive sessions and batch processing for optimal performance. Interactive sessions are used for exploratory data analysis and small-scale processing, while batch jobs are submitted for large-scale computations, such as alignment and peak calling. The pipeline provides scripts that can be adapted for either execution mode based on system requirements.  

Pipeline Overview  
---  
**1. Data Preprocessing**   
- Download FASTQ files from GEO  
- Quality control using **FastQC** & **MultiQC**

**3. Alignment & Peak Calling**  
- Indexing the reference genome for Chromap
- Align reads to hg38 genome using Chromap  
- Convert BED files to BAM using samtools  
- Sort and index BAM files with samtools  
- Generate BigWig tracks for visualization  
- Peak calling using **MACS2**

**3. Quality Control**  
- Compute FRiP score (Fraction of Reads in Peaks)  

**4. Visualization**  
- Plot FRiP scores using ggplot2  
- Generate heatmaps, correlation matrices & scatterplot using deepTools  
- Identify unique & shared peaks using bedtools  
- Generate Venn diagrams to compare peak overlaps

Repository Structure
---

```
ChIP-seq-analysis/
├── code/
│   ├── chromap_bed_to_bam_script.sh         # for converting the output from chromap alignemnt which is in bed format to the bam fromat
│   ├── chromap_alignment.script.sh          # for mapping the fastq files to the reference genome using Chromap
│   ├── correlation_matrix-deeptools.script.sh       
│   ├── frip-score-visualization.R           # R scripts for visualizing the frip score in barplots using ggplot2
|   ├── genome_index.script.sh       #for indexing the reference genome for Chromap alignment
|   ├── heatmap-deeptools.script.sh   #For making heatmap using deeptools
|   ├── peakcalling.script.sh            
├── results/
│   ├── FRiP_scores.png  
│   ├── Frip_scores.xlsx
│   ├── correlation_heatmap.png
│   ├── correlation_scatterplot.png
│   ├── heatmap.png
│   ├── multiqc_report.html
│   ├── peak number details.png
├── .gitignore
├── README.md              
```
**0. Interactive Session and Miniconda Setup**  
- Before starting the analysis, it is recommended to initialize an interactive session and set up the Miniconda environment for installing required tools.
- This setup is optional but recommended to streamline the analysis. It helps ensure a controlled environment, but you can modify or skip this step based on your system setup.
- Starting an interactive session`srun --nodes=1 --ntasks-per-node=1 --mem=20g --time=4:00:00 -p agsmall --pty bash -i`
- Once the session starts, download Miniconda (Linux x86) from the official website provided above and install it in your home directory.
- To keep dependencies organized, create a Conda environment specifically for ChIP-seq analysis`conda create -n chipseq`
- Activate the chipseq conda environment`conda activate chipseq`
- Using a Conda environment allows us to install all necessary tools for this project while keeping our system clean. However, you can install the tools globally if preferred.

**1. Install Required Software for the Analysis not in the modules**  
- Install deepTools: activate conda environment if not activated`conda activate chipseq`.Then install deepTools~conda install -c conda-forge -c bioconda deeptools`.
- Install MultiQC`conda install -c bioconda multiqc`
- Install chromap`conda install -c bioconda chromap`
- Install Samtools`conda install -c bioconda samtools`

**2. Downloading FASTQ Files from GEO/SRA**
- Navigate to the GEO page and use the SRA Run Selector to generate an Accession List for the FASTQ files you need and download the lsit or for this project you can directly download the`SRR_Acc_List.txt` file.
- Load the sratoolkit module to access SRA tools:`module load sratoolkit/3.0.0`
- Download FASTQ files using prefetch from SRA accession list(SRR_Acc_List.txt)`cat SRR_Acc_List.txt | xargs prefetch`.
- Convert downloaded SRA files to FASTQ format`cat SRR_Acc_List.txt | xargs fasterq-dump`. we have 4 fastq files for this project

**3. Quality Control Using FasttQC and MultiQC**  
- Once we have all the fastq files( which is 4 for this project), we will run FastQC to the the quality of the sequencing read. We will run this interactively as we only have 4 files, but you can submit it as batch job as well.
- Make sure to activate the chipseq conda environment `conda activate chipseq` and load the fastqc module `module load fastqc/0.12.1` or any verison of fastqc. the version 0.12.1 is just an example.
- Now run FastQC on all the input file `fastqc -o ../fastqc_output/ input_fastq/*.fq`.
- Fastqc will output an HTML report and a zipped archive containing the report's data and plots for each of the input fastq files.
- Next we will run MultiQC on .zip files from fastqc output `multiqc . -o ../multiqc_output/`

**4. Alignment/Mapping**  
**4.1 Reference Genome**  
- I used GRCh38.p14 as a reference genome for this project from GENCODE. You can use the same or use the reference genome from other website like ensembl or NCBI.
- To download the reference genome from GENCODE use
```
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/GRCh38.primary_assembly.genome.fa.gz
```  

 **4.2 Indexing the Reference Genome**  
 - The reference genome must be indexed before running the chromap alignment. The indexed genome from other alignment tools like STAR, Bowtie2 or BWA cannot be used for the Chromap. Each of these aligners should have their own indexed genome. We will submit job for creating the indexed genome for chromap as it is computationally intensive. The script for indexing is at genome_index.script.sh
```
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

#Indexing the reference genome
chromap -i -r /scratch.global/pande250/GCD-8141/Project-2/human_genome/GRCh38.p14.genome.fa.gz -o /scratch.global/pande250/GCD-8141/Project-2/indexed_human_genome/GRCh38_chromap_index
```
**4.2 Genome Mapping using Chromap**  
- Chromap (https://github.com/haowenz/chromap) besides aligning the reads to the reference genome it by default comes with all preprocessing steps like trimming the low quality reads and adapters and removing the duplicated reads.
- The full script for alignment can be found at chromap_alignment.script.sh.
- Basically for alignmetn this is the core script that we use:
```
chromap --preset chip \
        -x "$indexDir/GRCh38_chromap_index" \
        -r "$ref" \
        -q 20 \
        --min-read-length 15 \
        -1 "$fastqDir/$f" \
        -o "$bamDir/$samp.chromap.bed"
```
**4.3 Post mapping- Processing Chromap Outputs**  
- Chromap outputs in the bed format, we need to convert this into bam format using samtools and then perform sorting and indexing, indexing the referencing genome using samtools, and finally generating the bigwig files. script name for this process is chrom_bed_to_bam_script.sh.

```
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
```  




