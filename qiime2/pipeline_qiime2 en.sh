[TOC]

# QIIME2 2023.2 Analysis Process 

## 0. Software Installation (Appendix 1)

    # Only compatible with Linux/Mac systems，
    # For Windows users, use the built-in Linux Subsystem (supports right-click paste) or install and run in a Linux server environment
    # Refer to the official website for detailed instructions,https://docs.qiime2.org/2023.2/，Or the QIIME2 tutorial collection https://mp.weixin.qq.com/s/farGisfX3fVL_5WgXS8lIg
    # Install Windows Subsystem for Linux，https://mp.weixin.qq.com/s/0PfA0bqdvrEbo62zPVq4kQ

## 1. Preparatory Work

    # Set up the working directory, for example, on the server it would be ~/amplicon/qiime2, and for the Windows Subsystem：
    wd=/mnt/d/Easyamplicon/qiime2/
    #  Navigate to the Working Directory
    mkdir -p ${wd}
    cd ${wd}
    # Activate the QIIME2 working environment. For older versions of conda, use 'source' instead of 'conda activate
    conda activate qiime2-2023.2
    
    # Prepare sample metadata metadata.txt、Raw data seq/*.fq.gz
    
    ## The code provided here is based on the metadata.txt file to download sequencing data from public sources. It downloads data according to GSA's CRA (batch) and CRR (sample) identifiers
    mkdir -p seq
    # # downloadbload from the public Network
    # cd seq
    # awk '{system("wget -c ftp://download.big.ac.cn/gsa/"$5"/"$6"/"$6"_f1.fq.gz -O "$1"_1.fq.gz")}' \
    #     <(tail -n+2 ../metadata.txt)
    #  awk '{system("wget -c ftp://download.big.ac.cn/gsa/"$5"/"$6"/"$6"_r2.fq.gz -O "$1"_2.fq.gz")}' \
    #     <(tail -n+2 ../metadata.txt)
    # cd .. && ls -lsh seq
    # Link from another location (without taking up additional space)
    ln /mnt/c/EasyAmplicon/seq/* seq/
    ln /mnt/c/EasyAmplicon/result/metadata.txt ./
    # Generate a manifest file based on the metadata
    awk 'NR==1{print "sample-id\tforward-absolute-filepath\treverse-absolute-filepath"} \
      NR>1{print $1"\t$PWD/seq/"$1"_1.fq.gz\t$PWD/seq/"$1"_2.fq.gz"}' \
      metadata.txt > manifest
    head -n3 manifest
    
    # Import data into QIIME2 in paired-end 33 format
    qiime tools import \
      --type 'SampleData[PairedEndSequencesWithQuality]' \
      --input-path manifest \
      --output-path demux.qza \
      --input-format PairedEndFastqManifestPhred33V2
    # The fq file takes 7M for 1G, the fq.gz compression format is only 34s, and the test data is 9s


## 2. Generate feature tables and representative sequences
### Method 1: DADA2 (Slow, prone to errors due to multiple dependencies on R and Python packages)

    # Support multi-thread acceleration, 900,000 PE250 data，0/96p, 34m；24p, 44m；8p, 77m；1p, 462m
    # 270,000, 8p, 9m; 4p, 11m;
    time qiime dada2 denoise-paired \
      --i-demultiplexed-seqs demux.qza \
      --p-n-threads 4 \
      --p-trim-left-f 29 --p-trim-left-r 18 \
      --p-trunc-len-f 0 --p-trunc-len-r 0 \
      --o-table dada2-table.qza \
      --o-representative-sequences dada2-rep-seqs.qza \
      --o-denoising-stats denoising-stats.qza
    # Identify the results using DADA2 and import them into the main workflow
    cp dada2-table.qza table.qza
    cp dada2-rep-seqs.qza rep-seqs.qza

### Method2. External import of feature table and representative sequences (commonly used)

    # Upload the generated OTU table (otutab.txt) and representative sequences (otus.fa)
    # Convert the text file to Biom 1.0 format, keeping in mind that biom --version 2.1.5/8 works, but version 2.1.7 causes errors
    biom convert -i otutab.txt -o otutab.biom \
      --table-type="OTU table" --to-json
    # Importing the feature table took 9 seconds
    qiime tools import --input-path otutab.biom \
      --type 'FeatureTable[Frequency]' --input-format BIOMV100Format \
      --output-path table.qza
    # Importing the representative sequences took 8 seconds
    qiime tools import --input-path otus.fa \
      --type 'FeatureData[Sequence]' \
      --output-path rep-seqs.qza

### Statistics for the feature table and representative sequences

    qiime feature-table summarize \
      --i-table table.qza \
      --o-visualization table.qzv \
      --m-sample-metadata-file metadata.txt
    qiime feature-table tabulate-seqs \
      --i-data rep-seqs.qza \
      --o-visualization rep-seqs.qzv
    # Download the QZV file for online viewing. DADA2 resulted in just over 1,000 ASVs


## 3. Alpha and beta diversity analysis

### Build a phylogenetic tree for diversity analysis, took 53s

    qiime phylogeny align-to-tree-mafft-fasttree \
      --i-sequences rep-seqs.qza \
      --o-alignment aligned-rep-seqs.qza \
      --o-masked-alignment masked-aligned-rep-seqs.qza \
      --o-tree unrooted-tree.qza \
      --o-rooted-tree rooted-tree.qza

### Calculate core diversity 

    # 13s，Sampling depth is typically chosen as the minimum value, as indicated by the "table.qzv" file.
    qiime diversity core-metrics-phylogenetic \
      --i-phylogeny rooted-tree.qza \
      --i-table table.qza \
      --p-sampling-depth 7439 \
      --m-metadata-file metadata.txt \
      --output-dir core-metrics-results

### Significance analysis and visualization of inter-group alpha diversity
    # 7s, Optional alpha indices include faith_pd、shannon、observed_features、evenness
    index=observed_features
    qiime diversity alpha-group-significance \
      --i-alpha-diversity core-metrics-results/${index}_vector.qza \
      --m-metadata-file metadata.txt \
      --o-visualization core-metrics-results/${index}-group-significance.qzv

### Alpha diversity rarefaction curve

    # 25s, Choose the maximum value for max-depth, derived from table.qzv
    qiime diversity alpha-rarefaction \
      --i-table table.qza \
      --i-phylogeny rooted-tree.qza \
      --p-max-depth 10298 \
      --m-metadata-file metadata.txt \
      --o-visualization alpha-rarefaction.qzv
    # The available indices for the results are observed_otus, shannon, and faith_pd

### Significance of analysis and visualization of inter-group beta diversity

    # Optional beta diversity indices include unweighted_unifrac、bray_curtis、weighted_unifrac and jaccard
    # 7s, Specifying groups reduces computational load; permutation tests are more time-consuming
    distance=weighted_unifrac
    column=Group
    qiime diversity beta-group-significance \
      --i-distance-matrix core-metrics-results/${distance}_distance_matrix.qza \
      --m-metadata-file metadata.txt \
      --m-metadata-column ${column} \
      --o-visualization core-metrics-results/${distance}-${column}-significance.qzv \
      --p-pairwise


## 4. Species composition analysis

    # Species annotation using databases as listed in the appendix. You can refer to these databases,silva-138-99-nb-classifier.qza or 2022.10.backbone.full-length.nb.qza
    # 1m Optional specific primer training sets, such as classifier_gg_13_8_99_V5-V7.qza, are trained with V5-V7 regions. Refer to the appendix or official tutorial for more details
    time qiime feature-classifier classify-sklearn \
      --i-classifier 2022.10.backbone.full-length.nb.qza \
      --i-reads rep-seqs.qza \
      --o-classification taxonomy.qza
    # Visualize species annotations
    qiime metadata tabulate \
      --m-input-file taxonomy.qza \
      --o-visualization taxonomy.qzv
    # Display data with stacked bar charts
    qiime taxa barplot \
      --i-table table.qza \
      --i-taxonomy taxonomy.qza \
      --m-metadata-file metadata.txt \
      --o-visualization taxa-bar-plots.qzv


## 5. Differential analysis using ancom

    # Format the feature table and add pseudo-counts，4s
    qiime composition add-pseudocount \
      --i-table table.qza \
      --o-composition-table comp-table.qza
    
    # Calculate differential features, specifying the type of grouping for comparison，1m
    column=Group
    time qiime composition ancom \
      --i-table comp-table.qza \
      --m-metadata-file metadata.txt \
      --m-metadata-column ${column} \
      --o-visualization ancom-${column}.qzv
    
    # Merge by genus level and perform statistics
    ## Merge at the genus level，6s
    qiime taxa collapse \
      --i-table table.qza \
      --i-taxonomy taxonomy.qza \
      --p-level 6 \
      --o-collapsed-table table-l6.qza
    # Format the feature table and add pseudo-counts，6s
    qiime composition add-pseudocount \
      --i-table table-l6.qza \
      --o-composition-table comp-table-l6.qza
    # Calculate differential genera, specifying the type of grouping for comparison，16s
    qiime composition ancom \
      --i-table comp-table-l6.qza \
      --m-metadata-file metadata.txt \
      --m-metadata-column ${column} \
      --o-visualization ancom-l6-${column}.qzv


# Appendix

## 1. To install qiime2 2023.2 

### Install Conda

    # Download,Install, and  Launch conda
    wget -c https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh -b -f
    ~/miniconda3/condabin/conda init
    # Close the terminal and reopen it.

### Method 1. Online installation of QIIME using Conda.

    # Attached is the online installation of the software and packaging code  
    n=qiime2-2023.2
    # Download the list of Software
    wget -c https://data.qiime2.org/distro/core/${n}-py38-linux-conda.yml
    # Alternative Link
    wget -c http://www.imeta.science/db/conda/${n}-py38-linux-conda.yml
    # Install in a new environment, and once successfully installed on different computer servers, package and distribute.
    conda env create -n ${n} --file ${n}-py38-linux-conda.yml
    # Environment packaging(Optional，1.2G)
    conda pack -n ${n} -o ${n}.tar.gz

### Method 2. Local Installation of QIIME

    n=qiime2-2023.2
    # Installation package download link. 
    wget -c ftp://download.nmdc.cn/tools/conda/${n}.tar.gz
    # Installation in a new environment 
    mkdir -p ~/miniconda3/envs/${n}
    tar -xzf ${n}.tar.gz -C ~/miniconda3/envs/${n}
    # Activate and initialize the environment
    conda activate ${n}
    conda unpack

## 2. Species annotation data training Set

### Silva 138 99% OTUs full-length sequences

    # Download from the Official Website
    wget -c https://data.qiime2.org/2023.2/common/silva-138-99-nb-classifier.qza
    # Alternative Link
    wget -c ftp://download.nmdc.cn/tools/amplicon/silva/silva-138-99-nb-classifier.qza

### Greengenes2 2022.10 full length sequences

    # Download from the Official Website
    wget -c http://ftp.microbio.me/greengenes_release/2022.10/2022.10.backbone.full-length.nb.qza
    # Alternative Link
    wget -c ftp://download.nmdc.cn/tools/amplicon/GreenGenes/2022.10.backbone.full-length.nb.qza
    

## 3. Species Annotation data training set

    wd=/mnt/c/amplicon/qiime2
    mkdir -p $wd
    cd $wd
    # Download the Database file(greengenes, 320M)
    # wget -c ftp://greengenes.microbio.me/greengenes_release/gg_13_5/gg_13_8_otus.tar.gz
    # Domestic Alternate Links
    # wget -c ftp://download.nmdc.cn/tools/amplicon/GreenGenes/gg_13_8_otus.tar.gz
    # Backup of Core 99 Database within the country (60MB).
    wget -c ftp://download.nmdc.cn/tools/amplicon/GreenGenes/gg_13_8_otus_99.tar.gz
    mv gg_13_8_otus_99.tar.gz gg_13_8_otus.tar.gz
    # unpack
    tar -zxvf gg_13_8_otus.tar.gz
    
    # Use the 99_otus.fasta data in the rep_set file and the 99_OTU_taxonomy.txt data in taxonomy as reference species annotations
    # Import reference sequence, 50s
    qiime tools import \
      --type 'FeatureData[Sequence]' \
      --input-path gg_13_8_otus/rep_set/99_otus.fasta \
      --output-path 99_otus.qza
    # Fontconfig error: Cannot load default config file Does not affect the results
    
    # Importing species taxonomy informationtook 6s
    qiime tools import \
      --type 'FeatureData[Taxonomy]' \
      --input-format HeaderlessTSVTaxonomyFormat \
      --input-path gg_13_8_otus/taxonomy/99_otu_taxonomy.txt \
      --output-path ref-taxonomy.qza
    
    # Train the classifier——full-length，took 30m
    time qiime feature-classifier fit-classifier-naive-bayes \
      --i-reference-reads 99_otus.qza \
      --i-reference-taxonomy ref-taxonomy.qza \
      --o-classifier classifier_gg_13_8_99.qza

    # Primer extraction of amplified segments of the reference sequence Extract reference reads
    # Commonly use Greengenes 13_8 99% OTUs from 341F CCTACGGGNGGCWGCAG/805R GACTACHVGGGTATCTAATCC region of sequences（分类器描述），提供测序的引物序列，截取对应的区域进行比对，达到分类的目的。
    # This time used Primer 799F-1193R, Please replace according to the actual situation, 8m
    time qiime feature-classifier extract-reads \
      --i-sequences 99_otus.qza \
      --p-f-primer AACMGGATTAGATACCCKG \
      --p-r-primer ACGTCATCCCCACCTTCC \
      --o-reads ref-seqs.qza
    # Train the classifier
    # Based on the filtered specified segments, an experiment-specific classifier is generated，took 7 min
    time qiime feature-classifier fit-classifier-naive-bayes \
      --i-reference-reads ref-seqs.qza \
      --i-reference-taxonomy ref-taxonomy.qza \
      --o-classifier classifier_gg_13_8_99_V5-V7.qza
    
    # Common issue1：The scikit-learn version is not compatible, simply refit the classifier using Naive Bayes to build the training set will resolve the issue.
    Plugin error from feature-classifier:
      The scikit-learn version (0.21.2) used to generate this artifact does not match the current version of scikit-learn installed (0.22.1). Please retrain your classifier for your current deployment to prevent data-corruption errors.
    Debug info has been saved to /tmp/qiime2-q2cli-err-5ngzk2hm.log
    
    # Common issue2:  Error while running dada2 due to incomplete environment configuration. Run 'conda unpack' for initialization.
    Plugin error from dada2:
    An error was encountered while running DADA2 in R (return code 255), please inspect stdout and stderr to learn more.
    Debug info has been saved to /tmp/qiime2-q2cli-err-utwt1cmu.log
    
## Citation

    If used this script, please cited:
    
    Evan Bolyen, Jai Ram Rideout, Matthew R. Dillon, Nicholas A. Bokulich, Christian C. Abnet, Gabriel A. Al-Ghalith, Harriet Alexander, Eric J. Alm, Manimozhiyan Arumugam, Francesco Asnicar, Yang Bai, Jordan E. Bisanz, Kyle Bittinger, Asker Brejnrod, Colin J. Brislawn, C. Titus Brown, Benjamin J. Callahan, Andrés Mauricio Caraballo-Rodríguez, John Chase, Emily K. Cope, Ricardo Da Silva, Christian Diener, Pieter C. Dorrestein, Gavin M. Douglas, Daniel M. Durall, Claire Duvallet, Christian F. Edwardson, Madeleine Ernst, Mehrbod Estaki, Jennifer Fouquier, Julia M. Gauglitz, Sean M. Gibbons, Deanna L. Gibson, Antonio Gonzalez, Kestrel Gorlick, Jiarong Guo, Benjamin Hillmann, Susan Holmes, Hannes Holste, Curtis Huttenhower, Gavin A. Huttley, Stefan Janssen, Alan K. Jarmusch, Lingjing Jiang, Benjamin D. Kaehler, Kyo Bin Kang, Christopher R. Keefe, Paul Keim, Scott T. Kelley, Dan Knights, Irina Koester, Tomasz Kosciolek, Jorden Kreps, Morgan G. I. Langille, Joslynn Lee, Ruth Ley, **Yong-Xin Liu**, Erikka Loftfield, Catherine Lozupone, Massoud Maher, Clarisse Marotz, Bryan D. Martin, Daniel McDonald, Lauren J. McIver, Alexey V. Melnik, Jessica L. Metcalf, Sydney C. Morgan, Jamie T. Morton, Ahmad Turan Naimey, Jose A. Navas-Molina, Louis Felix Nothias, Stephanie B. Orchanian, Talima Pearson, Samuel L. Peoples, Daniel Petras, Mary Lai Preuss, Elmar Pruesse, Lasse Buur Rasmussen, Adam Rivers, Michael S. Robeson, Patrick Rosenthal, Nicola Segata, Michael Shaffer, Arron Shiffer, Rashmi Sinha, Se Jin Song, John R. Spear, Austin D. Swafford, Luke R. Thompson, Pedro J. Torres, Pauline Trinh, Anupriya Tripathi, Peter J. Turnbaugh, Sabah Ul-Hasan, Justin J. J. van der Hooft, Fernando Vargas, Yoshiki Vázquez-Baeza, Emily Vogtmann, Max von Hippel, William Walters, Yunhu Wan, Mingxun Wang, Jonathan Warren, Kyle C. Weber, Charles H. D. Williamson, Amy D. Willis, Zhenjiang Zech Xu, Jesse R. Zaneveld, Yilong Zhang, Qiyun Zhu, Rob Knight, J. Gregory Caporaso. 2019. Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2. **Nature Biotechnology** 37: 852-857. https://doi.org/10.1038/s41587-019-0209-9
    
    Copyright 2016-2023 Yong-Xin Liu <liuyongxin@caas.cn> 测试整理
