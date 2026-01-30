# EasyAmplicon

**Absolute Quantification Pipeline for Amplicon Sequencing Using Accu16S
/ AccuITS**

## Overview

**EasyAmplicon** is a reproducible Linux-based pipeline for **absolute
quantification of amplicon sequencing data** using synthetic **spike-in
standards (Accu16S / AccuITS)**.\
The workflow integrates spike-in read detection, ASV inference, and
standard-curve-based absolute copy number estimation at **DNA level**
and **sample level**.

The pipeline is designed for **paired-end Illumina amplicon data**
(e.g. 16S rRNA gene, ITS), and supports downstream analysis using
**QIIME2**, **BLAT**, **BLAST**, and **R**.

------------------------------------------------------------------------

## System Requirements

-   **Operating system**
    -   Linux (Ubuntu ≥ 22.04, CentOS ≥ 7.7)
    -   Windows 11 via **WSL (Windows Subsystem for Linux)** is
        supported
-   **Hardware**
    -   Multi-core CPU recommended (≥ 8 threads)
    -   ≥ 16 GB RAM recommended for large datasets

------------------------------------------------------------------------

## Software Dependencies

### Conda environment (basic tools)

``` bash
conda create -n meta_env -y
conda activate meta_env
conda install -c bioconda blat fastx_toolkit seqtk parallel -y
conda install -c conda-forge coreutils gzip awk gawk -y
```

### QIIME2 (Amplicon analysis)

``` bash
wget http://www.imeta.science/db/qiime2/qiime2-amplicon-ubuntu-latest-conda.yml
conda env create \
  --name qiime2-amplicon-2025.4 \
  --file qiime2-amplicon-ubuntu-latest-conda.yml

conda activate qiime2-amplicon-2025.4
qiime --help
```

### BLAST+

``` bash
wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.17.0+-x64-linux.tar.gz
tar -zxvf ncbi-blast-2.17.0+-x64-linux.tar.gz
export PATH=$PATH:/path/to/ncbi-blast-2.17.0+/bin
```

------------------------------------------------------------------------

## Directory Structure

``` text
EasyAmplicon/
├── seq/
├── db/
├── script/
├── temp/
├── result/
└── README.md
```

------------------------------------------------------------------------

## Pipeline Overview

### Step 0. Preparation

``` bash
sed -i 's/\r//;s/ //g' result/metadata.txt
```

### Step 1. Spike-in Read Identification and Separation

-   FASTQ → FASTA (seqtk)
-   BLAT alignment to spike-in reference
-   Alignment filtering
-   Spike-in read identification
-   FASTQ splitting (client vs spike-in reads)
-   Spike-in read counting and zero-detection

### Step 2. ASV Inference (QIIME2 + DADA2)

-   Spike-in reads ASV inference
-   Client reads ASV inference
-   Export ASV tables and representative sequences

### Step 3. Absolute Quantification

-   Spike-in IR copy extraction
-   BLAST-based spike-in ASV matching
-   Standard curve fitting
-   Absolute copy number calculation

``` bash
Rscript script/quantitation_calculate_absolute_copies.R \
  -o feature-table.txt \
  -s spike_otu.txt \
  -i spike_IR_copies.txt \
  -d DNA_quantity.txt \
  --output result/qiime2/6.final
```

------------------------------------------------------------------------

## Output Files

All final results are located in:

``` text
result/6.final/
```

  File                         Description
  ---------------------------- ------------------------------
  feature-table2.txt           ASV abundance table
  dna-sequences.fasta          ASV representative sequences
  standard_curve_formula.txt   Standard curve equations
  absolute_otu.txt             Absolute ASV copy numbers
  unit_dna_otu_copies.txt      DNA-normalized copies
  unit_sample_otu_copies.txt   Sample-level copies
  spike_otu.txt                Spike-in ASV table
  spike_reads_percent.txt      Spike-in recovery

------------------------------------------------------------------------

## Applications

-   Absolute quantification of microbial communities
-   Spike-in--based correction of sequencing bias
-   16S rRNA gene / ITS amplicon studies

------------------------------------------------------------------------

## Citation

If you use **EasyAmplicon**, please cite or acknowledge the pipeline
appropriately.
