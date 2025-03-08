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
ChIP-seq-analysis/
├── code/
│   ├── chromap_bed_to_bam_script.sh         # for converting the output from chromap alignemnt which is in bed format to the bam fromat
│   ├── chromap_alignment.script.sh          # for mapping the fastq files to the reference genome using Chromap
│   ├── correlation_matrix-deeptools.script.sh       
│   ├── frip-score-visualization.R           # R scripts for visualizing the frip score in barplots using ggplot2
|   ├── genome_index.script.sh       #for indexing the reference genome for Chromap alignment
|   ├── heatmap-deeptools.script.sh   #For making heatmap using deeptools
|   ├── peakcalling.script.sh
├── README.md             
├── results/
│   ├── FRiP_scores.png  
│   ├── Frip_scores.xlsx
│   ├── correlation_heatmap.png
│   ├── correlation_scatterplot.png
│   ├── heatmap.png
│   ├── multiqc_report.html
│   ├── peak number details.png
├── .gitignore                
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

**3. Quality Control Using FastTQC and MultiQC**  
-  

- 



