# 易扩增子 (EasyAmplicon 2) 流程教程-二代扩增子分析
EasyAmplicon Tutorial - Second-Generation Amplicon Analysis

**简介:** EasyAmplicon 是一个跨平台且社区维护的开源扩增子序列分析流程，集成了数据质量控制、双端序列合并、去冗余、聚类/去噪、嵌合体检测、特征表构建、物种注释、多样性分析和可视化等模块 。该管线可在 RStudio 中一站式运行，显著降低了学习成本。  
*EasyAmplicon is an open-source, cross-platform pipeline for amplicon sequencing data analysis. It includes modules for quality control, paired-end read merging, dereplication, clustering/denoising, chimera detection, feature table construction, taxonomic analysis, and visualization .*

## 1. 系统要求 / System Requirements
- **操作系统:** Windows 10 或更新版本；macOS 10.12 或以上；Ubuntu 20.04 或更高。  
  **OS:** Windows 10+; MacOS 10.12+; or Ubuntu 20.04+.
- **软件依赖:** R (建议 RStudio 环境)、vsearch、usearch、QIIME 2、muscle、iqtree/fasttree 等。  
  **Software:** R (with RStudio), **vsearch**, **usearch**, **QIIME 2**, **MUSCLE**, **IQ-TREE**/**FastTree**, etc.
- **数据库:** 内置了常用 16S/18S/ITS 序列数据库（RDP、SILVA、UNITE）和 Greengenes 97% OTUs 等。首次使用时需要解压并准备这些参考数据库。  
  **Databases:** Pre-packaged reference databases (e.g., RDP/SILVA/UNITE for taxonomy, Greengenes 97% OTUs for PICRUSt) should be downloaded and uncompressed prior to use.

## 2. 输入数据准备 / Input Data Preparation
- **分析脚本:** `pipeline.sh` (或相应的 RStudio 脚本) 是主要的流程控制脚本。  
  **Pipeline script:** `pipeline.sh` (or R scripts) controls the overall workflow.
- **样本元数据 (`metadata.txt`):** 保存每个样本的设计信息。至少需要包含 SampleID 列（与序列文件名对应）和分组信息（例如 Description 或 Group 列）。文件以制表符分隔。  
  **Metadata:** Prepare a sample metadata file (`metadata.txt`) with at least columns `SampleID` (matching sequence file names) and group/description info. Use tab-delimited format.
- **测序数据:** 将配对或单端的 FASTQ 文件放入 `seq/` 目录。文件名格式建议为 `SampleID_1.fq.gz`（前向）和 `SampleID_2.fq.gz`（反向）。务必确保文件名与元数据中的 SampleID 一致。  
  **Sequencing data:** Place paired-end or single-end FASTQ files in `seq/`. Name files like `SampleID_1.fq.gz` and `SampleID_2.fq.gz`. Ensure sample IDs match metadata.
- **工作目录结构:** 创建下列文件夹：`seq/`（原始数据）、`result/`（输出结果）、`temp/`（临时文件）。  
  **Directories:** Create directories `seq/` (input), `result/` (results), and `temp/` (temporary files).

## 3. 序列合并与重命名 / Merge Paired-End Reads
- **合并双端序列:** 使用 **vsearch** 的 `--fastq_mergepairs` 命令将每个样本的前后向序列合并，并用样本 ID 重命名序列名称。
```bash
vsearch --fastq_mergepairs seq/Sample1_1.fq.gz --reverse seq/Sample1_2.fq.gz \
  --fastqout temp/Sample1.merged.fq --relabel Sample1.
```
- **批量处理:** 对所有样本进行批量处理，可用 `for` 循环遍历元数据中的 SampleID。完成后检查 `.merged.fq` 文件的序列名。
- **单端数据:** 使用 `vsearch --fastq_convert` 或 `usearch -fastx_relabel` 给序列添加样本前缀。
- **整合所有样本:** 将所有样本合并到一个文件：
```bash
cat temp/*.merged.fq > temp/all.fq
```

## 4. 引物切除与质控 / Primer Trimming & Quality Filtering
- **切除引物:** 使用 `vsearch --fastx_filter` 并指定 `--fastq_stripleft` 和 `--fastq_stripright`。
```bash
vsearch --fastx_filter temp/all.fq \
  --fastq_stripleft 29 --fastq_stripright 18 \
  --fastq_maxee_rate 0.01 \
  --fastaout temp/filtered.fa
```
- **质量过滤:** `--fastq_maxee_rate 0.01` 移除高错误率序列。
- **检查结果:** 使用 `head` 查看 `filtered.fa` 确认格式正确。

## 5. 去冗余与OTU/ASV 构建 / Dereplication & OTU/ASV Clustering
- **序列去冗余:**
```bash
vsearch --derep_fulllength temp/filtered.fa \
  --minuniquesize 10 --sizeout --relabel Uni_ \
  --output temp/uniques.fa
```
- **ASV 去噪:** UNOISE3：
```bash
usearch -unoise3 temp/uniques.fa -minsize 10 -zotus temp/zotus.fa
sed 's/Zotu/ASV_/g' temp/zotus.fa > temp/otus.fa
```
- **OTU 聚类（可选）:**
```bash
usearch -cluster_otus temp/uniques.fa -minsize 10 -otus temp/otus.fa -relabel OTU_
```
- **vsearch 替代:**
```bash
vsearch --cluster_size temp/uniques.fa --id 0.97 --centroids temp/otus_raw.fa
vsearch --uchime_denovo temp/otus_raw.fa --nonchimeras temp/otus.fa
```

## 6. 去嵌合 / Chimera Removal
- **参考去嵌合:**
```bash
vsearch --uchime_ref temp/otus.fa --db ${db}/usearch/rdp_16s_v18.fa --nonchimeras result/raw/otus.fa
```
- **UNOISE3 已去嵌合:**
```bash
cp temp/otus.fa result/raw/otus.fa
```

## 7. 特征表构建与过滤 / Feature Table Construction & Filtering
- **构建特征表:**
```bash
vsearch --usearch_global temp/filtered.fa --db result/raw/otus.fa \
  --id 0.97 --otutabout result/raw/otutab.txt
```
- **去除非细菌:**
```bash
vsearch --sintax result/raw/otus.fa --db ${db}/usearch/rdp_16s_v18.fa \
       --sintax_cutoff 0.1 --tabbedout result/raw/otus.sintax
Rscript ${db}/script/otutab_filter_nonBac.R \
  --input result/raw/otutab.txt --taxonomy result/raw/otus.sintax \
  --output result/otutab.txt --discard result/raw/otus.sintax.discard
```
- **格式转换:**
```bash
cut -f1 result/otutab.txt | tail -n+2 > result/otutab.id
usearch -fastx_getseqs result/raw/otus.fa -labels result/otutab.id -fastaout result/otus.fa
```
- **抽样归一化:**
```bash
Rscript ${db}/script/otutab_rare.R --input result/otutab.txt \
  --depth 10000 --seed 1 --normalize result/otutab_rare.txt \
  --output result/alpha/vegan.txt
```

## 8. 多样性分析 / Diversity Analysis
- **α 多样性:**
```bash
usearch -alpha_div result/otutab_rare.txt -output result/alpha/alpha.txt
usearch -alpha_div_rare result/otutab_rare.txt -output result/alpha/alpha_rare.txt -method without_replacement
Rscript ${db}/script/alpha_boxplot.R --alpha_index richness --input result/alpha/vegan.txt --design result/metadata.txt --group Group --output result/alpha/ --width 89 --height 59
Rscript ${db}/script/alpha_rare_curve.R --input result/alpha/alpha_rare.txt --design result/metadata.txt --group Group --output result/alpha/ --width 120 --height 59
```
- **β 多样性:**
```bash
usearch -cluster_agg result/otus.fa -treeout result/otus.tree
usearch -beta_div result/otutab_rare.txt -tree result/otus.tree -filename_prefix result/beta/
Rscript ${db}/script/beta_pcoa.R --input result/beta/bray_curtis.txt --design result/metadata.txt --group Group --output result/beta/bray_curtis.pcoa.pdf
```

## 9. 物种注释与组成 / Taxonomy Assignment & Composition
- **生成注释表:** `result/taxonomy.txt`
- **各级丰度统计:**
```bash
for rank in p c o f g; do
  usearch -sintax_summary result/otus.sintax -otutabin result/otutab_rare.txt \
    -rank ${rank} -output result/tax/sum_${rank}.txt
done
```
- **可视化:**
```bash
Rscript ${db}/script/tax_stackplot.R --input result/tax/sum_p.txt --design result/metadata.txt --group Group --output result/tax/ --width 89 --height 59
Rscript ${db}/script/tax_circlize.R --input result/tax/sum_c.txt --design result/metadata.txt --group Group --legend 5
```

## 10. 功能预测 / Functional Prediction
- **PICRUSt2:**
```bash
conda activate picrust2
picrust2_pipeline.py -s result/otus.fa -i result/otutab.txt -o result/picrust2_out -p 8
```
- **FAPROTAX:**
```bash
biom convert -i result/otutab_rare.txt -o otutab.biom --table-type "OTU table" --to-json
biom add-metadata -i otutab.biom --observation-metadata-fp result/taxonomy.txt -o otutab_tax.biom --sc-separated taxonomy
python $FAPROTAX/collapse_table.py -i otutab_tax.biom -g $FAPROTAX/FAPROTAX.txt --collapse_by_metadata taxonomy -o faprotax.txt -r faprotax_report.txt
```
- **BugBase:**
```bash
Rscript ${db}/script/BugBase/bin/run.bugbase.r -L ${db}/script/BugBase -i result/gg/otutab.txt -m result/metadata.txt -c Group -o result/bugbase/
```

## 11. 差异丰度分析 / Differential Abundance Analysis
```bash
Rscript ${db}/script/compare.R --input result/otutab.txt --design result/metadata.txt --group Group --compare KO-WT --method wilcox --pvalue 0.05 --fdr 0.2 --output result/compare/
Rscript ${db}/script/compare_volcano.R --input result/compare/KO-WT.txt --output result/compare/KO-WT.volcano.pdf
bash ${db}/script/compare_heatmap.sh -i result/compare/KO-WT.txt -l 7 -d result/metadata.txt -A Group -t result/taxonomy.txt -w 8 -h 5 -s 7 -o result/compare/KO-WT
bash ${db}/script/compare_manhattan.sh -i result/compare/KO-WT.txt -t result/taxonomy.txt -p result/tax/sum_c.txt -w 183 -v 59 -s 7 -l 10 -L Class -o result/compare/KO-WT.manhattan.c.pdf
Rscript ${db}/script/alpha_boxplot.R --alpha_index ASV_2 --input result/otutab.txt --design result/metadata.txt --transpose TRUE --scale TRUE --group Group --output result/compare/feature_ASV2
```

## 12. 机器学习分析 / Machine Learning
```bash
conda activate qiime2-amplicon-2024.2
cd ${db}/script/slime2
./slime2.py otutab.txt design.txt --normalize --tag rf_e4 rf -n 10000
./slime2.py otutab.txt design.txt --normalize --tag ab_e4 ab -n 10000
```

## 13. 系统发育树构建与注释 / Phylogenetic Tree Construction & Annotation
```bash
usearch -otutab_trim result/otutab_rare.txt -min_otu_freq 0.005 -output tree/otutab.txt
cut -f1 tree/otutab.txt > tree/otutab_high.id
usearch -fastx_getseqs result/otus.fa -labels tree/otutab_high.id -fastaout tree/otus.fa
muscle -in tree/otus.fa -out tree/otus_aligned.fas
iqtree -s tree/otus_aligned.fas -bb 1000 -alrt 1000 -nt AUTO -pre tree/iqtree/otus
fasttree -gtr -nt tree/otus_aligned.fas > tree/otus.nwk
```
*iTOL 注释使用 `table2itol.R` 生成。*

## 14. 附加工具与分析资源
- QIIME 2: `qiime2/pipeline_qiime2.sh`
- 视频教程: Meta-genome 微信公众号及课程页
- 常见问题: FAQ 提供 Fastq 质量编码、重命名、多线程、Windows 换行符等解决方案

**参考文献:**

使用此脚本，请引用下文：

If used this script, please cited:

Yong-Xin Liu, Lei Chen, Tengfei Ma, Xiaofang Li, Maosheng Zheng, Xin Zhou, Liang Chen, Xubo Qian, Jiao Xi, Hongye Lu, Huiluo Cao, Xiaoya Ma, Bian Bian, Pengfan Zhang, Jiqiu Wu, Ren-You Gan, Baolei Jia, Linyang Sun, Zhicheng Ju, Yunyun Gao, Tao Wen, Tong Chen. 2023. EasyAmplicon: An easy-to-use, open-source, reproducible, and community-based pipeline for amplicon data analysis in microbiome research. iMeta 2: e83. https://doi.org/10.1002/imt2.83
