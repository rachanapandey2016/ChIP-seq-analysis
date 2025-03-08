# Chip-seq-analysis  
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
- Convert BED files to BAM using bedtools  
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
ChIP-seq_GATA3_JUN_Enhancer/
├── code/
│   ├── chipseq_analysis.sh   # Full pipeline script
│   ├── preprocessing.sh      # FASTQ to BAM conversion
│   ├── peak_calling.sh       # MACS2 peak calling
│   ├── visualization.R       # R scripts for plots
├── data/
│   ├── README.md             # Dataset description (No raw data)
├── results/
│   ├── figures/              # Plots (heatmaps, Venn diagrams)
│   ├── summary_tables/       # FRiP scores, peak stats
├── README.md                 # Main project documentation
├── LICENSE                   # Open-source license (optional)
├── .gitignore                # Ignore large files (e.g., BAM, FASTQ)
```
**0) Interactive Session and Miniconda Setup**  
- Before starting the analysis, it is recommended to initialize an interactive session and set up the Miniconda environment for installing required tools.
- This setup is optional but recommended to streamline the analysis. It helps ensure a controlled environment, but you can modify or skip this step based on your system setup.
- Starting an interactive session`srun --nodes=1 --ntasks-per-node=1 --mem=20g --time=4:00:00 -p agsmall --pty bash -i`
- Once the session starts, download Miniconda (Linux x86) from the official website provided above and install it in your home directory.
- To keep dependencies organized, create a Conda environment specifically for ChIP-seq analysis`conda create -n chipseq`
- Activate the chipseq conda environment`conda activate chipseq`
- Using a Conda environment allows us to install all necessary tools for this project while keeping our system clean. However, you can install the tools globally if preferred.

**1) Install Required Software for the Analysis**  
- Install deepTools




