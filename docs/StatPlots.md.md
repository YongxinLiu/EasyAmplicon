
# 易扩增子 (EasyAmplicon 2) 流程教程-StatPlots

**易扩增子 (EasyAmplicon 2)** 是一个跨平台、开源的宏基因组扩增子序列分析流程，旨在简化生物群落测序数据的分析过程。该管道集成了30多个常用模块，涵盖从原始数据质控、序列合并，到去重、聚类/去噪、嵌合体检测、特征表生成、分类学多样性与组成分析，以及生物标志物发现和高质量可视化等几乎所有分析环节。版本2.0 相较于1.0版，扩展支持第三代全长扩增子测序数据（如 PacBio、Nanopore、Qitan 等平台），并集成了DADA2、Emu等主流工具，实现了从多种测序平台产生的原始序列到可直接发表的可视化图形的一站式工作流。易扩增子采用R语言和脚本化操作，用户可通过RStudio调用各类分析脚本，同时也提供了命令行工具，方便批量化处理。
(EasyAmplicon 2 is a cross-platform, open-source metagenomic amplicon sequence analysis pipeline designed to simplify the analysis of biome sequencing data. The pipeline integrates over 30 commonly used modules, covering nearly every analytical step, from raw data quality control and sequence merging to duplicate removal, clustering/denoising, chimera detection, feature table generation, taxonomic diversity and composition analysis, biomarker discovery, and high-quality visualization. Compared to version 1.0, version 2.0 expands support for third-generation, full-length amplicon sequencing data (e.g., from PacBio, Nanopore, and Qitan platforms) and integrates mainstream tools such as DADA2 and Emu, enabling a one-stop workflow from raw sequences generated from various sequencing platforms to publication-ready visualizations. EasyAmplicon utilizes the R programming language and scripting capabilities, allowing users to call various analysis scripts through RStudio. Command-line tools are also provided for batch processing.)

管道可部署于Windows、Mac OS和Linux系统，安装使用简单迅速。以下以易扩增子2.0版为例，介绍其主要功能模块与使用方法（所有命令均以Rscript或bash脚本形式运行）。
(The pipeline can be deployed on Windows, Mac OS, and Linux systems, making installation and use quick and easy. Below, using EasyAmplicon version 2.0 as an example, we'll introduce its key functional modules and usage instructions (all commands are run as Rscript or bash scripts).)

## 安装与环境配置 Installation and environment configuration

在使用易扩增子前，需要设置工作目录和软件库路径，并将脚本路径添加到环境变量中。以下示例假设将软件安装于`EasyMicrobiome`文件夹中，工作目录为`/c/amplicon2`（Windows）；Linux/Mac用户路径可相应调整。每次在RStudio或终端中运行前，执行以下命令：
(Before using EasyMicrobiome, you need to set your working directory and library path, and add the script path to your environment variables. The following examples assume the software is installed in the `EasyMicrobiome` folder, and the working directory is `/c/amplicon2` (Windows); Linux/Mac users can adjust the path accordingly. Before running each script in RStudio or Terminal, execute the following command:)
```bash
wd=/c/amplicon2
=db=/c/EasyMicrobiome
export PATH=$PATH:${db}/win
cd ${wd}
```

完成环境配置后，可根据数据类型（Illumina、PacBio或Nanopore）进入相应子目录，执行后续分析命令。
(After completing the environment configuration, you can enter the corresponding subdirectory according to the data type (Illumina, PacBio, or Nanopore) to execute subsequent analysis commands.)

## 多样性分析 (Diversity Analysis)

易扩增子提供丰富的微生物多样性分析工具，支持**Alpha多样性**和**Beta多样性**可视化。
(EasyAmplicon provides a wealth of microbial diversity analysis tools and supports visualization of **Alpha diversity** and **Beta diversity**.)

### Alpha多样性 Alpha diversity

#### 1. 箱线图 (Boxplot)
```bash
Rscript ${db}/script/alpha_div_box.R --input result/alpha/vegan.txt --metadata result/metadata.txt --alpha_index richness,chao1,ACE,shannon,simpson,invsimpson --group Group --out_prefix result/alpha/vegan
```

#### 2. 稀释曲线 (Rarefaction Curve)
```bash
Rscript ${db}/script/alpha_rare_curve.R --input result/alpha/alpha_rare.txt --design result/metadata.txt --group Group --output result/alpha/ --width 120 --height 78
```

#### 3. Venn图 Venn plot
```bash
Rscript ${db}/script/venn.R --input result/alpha/otu_group_exist.txt --groups All,KO,OE,WT --output result/alpha/venn.pdf
```

### Beta多样性 Beta diversity

#### 1. 距离矩阵热图 (pheatmap)
```bash
bash ${db}/script/sp_pheatmap.sh -f result/beta/bray_curtis.txt -H 'TRUE' -u 6 -v 5
```

#### 2. 主坐标分析 (PCoA)
```bash
Rscript ${db}/script/beta_PCoA.R --input result/otutab.txt --metadata result/metadata.txt --group Group --output result/beta/PCoa2.pdf
```

#### 3. 限制性主坐标分析 (CPCoA)
```bash
Rscript ${db}/script/beta_cpcoa.R --input result/beta/bray_curtis.txt --design result/metadata.txt --group Group --output result/beta/bray_curtis.cpcoa.pdf --width 89 --height 59
```

## 物种组成 (Taxonomy) 分析 Species composition (Taxonomy) analysis

### 堆叠柱状图 (Stacked Bar Plot)
```bash
Rscript ${db}/script/tax_stackplot.R --input result/tax/sum_p.txt --design result/metadata.txt --group Group -t 10 --color manual1 --legend 7 --width 89 --height 59 --output result/tax/sum_p.stackplot
```

### 圆弦图 (Circos/Circlize Plot)
```bash
i=c
Rscript ${db}/script/tax_circlize.R --input result/tax/sum_${i}.txt --design result/metadata.txt --group Group --legend 5
mv circlize.pdf result/tax/sum_${i}.circlize.pdf
mv circlize_legend.pdf result/tax/sum_${i}.circlize_legend.pdf
```

### 气泡图 (Bubble Plot)
```bash
i=g
Rscript ${db}/script/tax_bubble.R -i result/tax/sum_${i}.txt -g result/metadata.txt -c Group -w 7 -e 4 -s 15 -n 15 -o result/tax/sum_g.bubble3.pdf
```

### 甜甜圈图 (Donut) & 雷达图 (Radar Plot)
```bash
Rscript ${db}/script/Donut_plot.R --otu_table result/otutab.txt --metadata result/metadata.txt --taxonomy result/taxonomy.txt --output_dir result/tax/
```

### 核心微生物 (Core Microbiome)
```bash
Rscript ${db}/script/core_ASVS_and_other.R --otu_table result/otutab.txt --metadata result/metadata.txt --taxonomy result/taxonomy.txt --output_dir result/tax/
Rscript ${db}/script/Core_Abundance_ScatterPlot.R --otu_table result/otutab.txt --metadata result/metadata.txt --output result/tax/
```

## 差异分析 (Differential Analysis)

### 火山图 (Volcano Plot)
```bash
Rscript ${db}/script/volcano2.R --input result/compare/KO-OE2.txt --group KO-OE --output result/compare/
```

### 热图
```bash
bash ${db}/script/compare_heatmap.sh -i result/compare/KO-OE.txt -l 7 -d result/metadata.txt -A Group -t result/taxonomy.txt -w 12 -h 20 -s 14 -o result/compare/KO-OE
```

### 三元图 (Ternary Plot)
```bash
Rscript ${db}/script/Ternary_plot.R --input result/tax/sum_p2.txt --metadata result/metadata2.txt --group Group --taxlevel Phylum --output result/compare/ternary_p.pdf --topn 10
```

### STAMP差异分析
```bash
Rscript ${db}/script/compare_stamp.R --input result/tax/sum_g.txt --metadata result/metadata.txt --group Group --compare KO-OE --threshold 0.1 --method "t.test" --pvalue 0.2 --fdr "none" --width 80 --height 30 --output result/tax/stamp_KO-OE
```

## 网络分析 (Network Analysis)
```bash
Rscript ${db}/script/Spearman_network01.R --input result/tax/otutab_amplicon.txt --group result/tax/sample_amplicon.txt --tax result/tax/taxonomy_amplicon.txt --output result/tax/
```

## 随机森林模型 (Random Forest)
```bash
Rscript ${db}/script/random_forest01.R --input result/tax/PacBio_data.txt --group result/tax/PacBio_metadata.txt --output result/tax/RF_model/
Rscript ${db}/script/RF_plot01.R --input result/tax/RF_model/Species_imp_rf21.txt --tax result/tax/taxonomy_amplicon.txt --optimal 4 --output result/tax/RF_model/
```

## 功能差异分析 (Functional Analysis)
```bash
Rscript ${db}/script/function_data_process.R --input result/tax/KEGG.PathwayL2.raw.txt --group result/tax/metadata_amplicon.txt --output result/tax/pathway/
Rscript ${db}/script/function_diff_plot01.R --input result/tax/pathway/Difference_pathway21.txt --pathway result/tax/pathway/pathway_count_data.txt --statistic result/tax/pathway/pathway_percent_abundance2.txt --output result/tax/pathway/
```

References:

使用此脚本，请引用下文：

If used this script, please cited:

Yong-Xin Liu, Lei Chen, Tengfei Ma, Xiaofang Li, Maosheng Zheng, Xin Zhou, Liang Chen, Xubo Qian, Jiao Xi, Hongye Lu, Huiluo Cao, Xiaoya Ma, Bian Bian, Pengfan Zhang, Jiqiu Wu, Ren-You Gan, Baolei Jia, Linyang Sun, Zhicheng Ju, Yunyun Gao, Tao Wen, Tong Chen. 2023. EasyAmplicon: An easy-to-use, open-source, reproducible, and community-based pipeline for amplicon data analysis in microbiome research. iMeta 2: e83. https://doi.org/10.1002/imt2.83