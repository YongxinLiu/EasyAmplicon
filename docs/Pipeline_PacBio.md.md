
# 易扩增子 (EasyAmplicon 2) 流程教程-PacBio Pipeline
EasyAmplicon 2 Tutorial - PacBio Pipeline

本工作流处理全长扩增子测序数据（如 PacBio 的 16S），生成 OTU 或 ASV 表格、分类注释以及多样性分析结果。使用 Bash/Conda 环境中的 USEARCH、VSEARCH、Cutadapt、seqkit 等命令行工具，并辅以 R 脚本进行元数据处理和过滤。主要步骤如下。

This workflow processes full-length amplicon sequencing data (such as PacBio's 16S) to generate OTU or ASV tables, taxonomic annotations, and diversity analysis results. It uses command-line tools such as USEARCH, VSEARCH, Cutadapt, and seqkit in a Bash/Conda environment, supplemented by R scripts for metadata processing and filtering. The main steps are as follows.

## 1. 设置与输入数据 / Setup and Input Data

- **工作目录与环境 / Working directory and environment:**
```bash
wd=/c/EasyAmplicon2
db=/c/EasyMicrobiome
PATH=$PATH:${db}/win
cd ${wd}
```
激活 Conda 环境 / Activate the Conda environment:
```bash
conda activate easyAmplicon2
```

- **元数据 / Metadata:** 准备 `result/metadata.txt`，至少包含 `SampleID` 与 `Description`。如需去掉 Windows 回车符：
- **Metadata:** Prepare `result/metadata.txt`, which contains at least `SampleID` and `Description`. To remove Windows carriage returns:
```bash
sed -i 's/\r//' result/metadata.txt
```

- **测序数据 / Sequencing reads:** 将原始 PacBio reads 放入 `seq/`，统计文件：
- **Sequencing reads:** Put raw PacBio reads into `seq/`, statistics files:
```bash
seqkit stat seq/*.fastq > result/seqkit.txt
```

- **参考数据库 / Reference databases:** 确保 16S 数据库（如 SILVA）已解压可用。
- **Reference databases:** Make sure the 16S database (such as SILVA) is unpacked and available.

## 2. Reads 合并与重命名 / Read Merging and Renaming
```bash
for i in $(cut -f1 result/metadata.txt | tail -n+2); do
  usearch -fastx_relabel seq/${i}.fastq -fastqout temp/${i}.merged.fastq -prefix ${i}.
done
cat temp/*.merged.fastq > temp/all.fq
```

## 3. 引物去除与质量过滤 / Primer Trimming and Quality Filtering

- **Cutadapt 去除引物 / Primer removal:**
```bash
cutadapt -g "AGRGTTYGATYMTGGCTCAG...RGYTACCTTGTTACGACTT" \
         --error-rate=0.1 -j 10 --discard-untrimmed \
         -o temp/allfilter.fastq temp/all.fq
```

- **VSEARCH 长度过滤 / Length filtering:**
```bash
vsearch --fastx_filter temp/allfilter.fastq \
        --fastq_minlen 1200 --fastq_maxlen 1800 \
        --fastaout temp/filtered.fa --fastq_qmax 93
```

## 4. 去重与聚类/去噪 / Dereplication and Clustering/Denoising

- **去重 / Dereplication:**
```bash
vsearch --derep_fulllength temp/filtered.fa --fasta_width 0 \
        --sizeout --relabel Uni_ --minuniquesize 2 \
        --threads 8 --output temp/uniques.fa
```

- **OTU 聚类 (97%) 或 ASV 去噪 (UNOISE3) / OTU clustering (97%) or ASV denoising (UNOISE3):**
```bash
usearch -cluster_otus temp/uniques.fa -minsize 2 -otus temp/otus.fa -relabel OTU_
# 或 / or
usearch -unoise3 temp/uniques.fa -minsize 2 -zotus temp/zotus.fa
sed 's/Zotu/ASV_/g' temp/zotus.fa > temp/otus.fa
```

## 5. Feature 表生成与过滤 / Feature Table Generation and Filtering
```bash
vsearch --usearch_global temp/filtered.fa --db result/raw/otus.fa \
        --id 0.97 --threads 4 --otutabout result/raw/otutab.txt
```
- 可选：去除叶绿体/线粒体 / Remove chloroplasts/mitochondria (optional) 使用 SINTAX 分类器及 `otutab_filter_nonBac.R`

## 6. Alpha 与 Beta 多样性 / Alpha and Beta Diversity

- **稀释 / Rarefaction:**
```bash
Rscript ${db}/script/otutab_rare.R --input result/otutab.txt --depth 327 --seed 1 \
  --normalize result/otutab_rare.txt --output result/alpha/vegan.txt
```

- **Alpha 多样性 / Alpha diversity:**
```bash
usearch -alpha_div result/otutab_rare.txt -output result/alpha/alpha.txt
```

- **Beta 多样性 / Beta diversity:**
```bash
usearch -cluster_agg result/otus.fa -treeout result/otus.tree
usearch -beta_div result/otutab_rare.txt -tree result/otus.tree -filename_prefix result/beta/
```

## 7. 分类汇总 / Taxonomic Summaries

- **OTU-to-taxon 映射 / OTU-to-taxon map:**
```bash
cut -f1,4 result/otus.sintax | sed 's/\td/\tk/;s/:/__ /g;s/,/;/g;s/"//g' > result/taxonomy2.txt
```

- **分类计数 / Taxonomic counts:**
```bash
usearch -sintax_summary result/otus.sintax -otutabin result/otutab_rare.txt -rank p -output result/tax/sum_p.txt
```

## 8. 参考比对定量 (可选, PICRUSt2) / Reference-Based Quantification (Optional, PICRUSt2)
```bash
mkdir -p result/gg
usearch -otutab temp/filtered.fa -otus ${db}/gg/97_otus.fa -otutabout result/gg/otutab.txt -threads 4
```
- **PICRUSt2 分析 / PICRUSt2 pipeline:**
```bash
picrust2_pipeline.py -s result/otus.fa -i result/otutab.txt -o result/picrust2/out -p 8
```

## 9. 选定特征的系统发育树 / Phylogenetic Tree of Selected Features

1. 提取选定 OTU / Extract selected OTUs:
```bash
usearch -fastx_getseqs result/otus.fa -labels selected_otus.txt -fastaout selected_otus.fa
```

2. MUSCLE 比对 / Align sequences with MUSCLE:
```bash
muscle -in selected_otus.fa -out otus_aligned.fas
```

3. 构建树 / Build tree (FastTree 或 IQ-TREE):
```bash
fasttree -gtr -nt otus_aligned.fas > otus.nwk
# or
iqtree -s otus_aligned.fas -bb 1000 -alrt 1000 -nt AUTO -pre iqtree/otus
```

4. 可视化 / Visualize in iTOL or ETE.

---

本工作流将原始 PacBio 扩增子数据转化为计数表、多样性指标、系统发育树及预测功能谱，便于后续分析。
This workflow converts raw PacBio amplicon data into count tables, diversity metrics, phylogenetic trees, and predicted functional profiles for subsequent analysis.

**References:**

使用此脚本，请引用下文：

If used this script, please cited:

Yong-Xin Liu, Lei Chen, Tengfei Ma, Xiaofang Li, Maosheng Zheng, Xin Zhou, Liang Chen, Xubo Qian, Jiao Xi, Hongye Lu, Huiluo Cao, Xiaoya Ma, Bian Bian, Pengfan Zhang, Jiqiu Wu, Ren-You Gan, Baolei Jia, Linyang Sun, Zhicheng Ju, Yunyun Gao, Tao Wen, Tong Chen. 2023. EasyAmplicon: An easy-to-use, open-source, reproducible, and community-based pipeline for amplicon data analysis in microbiome research. iMeta 2: e83. https://doi.org/10.1002/imt2.83

