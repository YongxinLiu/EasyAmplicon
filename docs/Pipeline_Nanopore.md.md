
# 易扩增子 (EasyAmplicon 2) 流程教程-Nanopore Pipeline

**EasyAmplicon 2** 是一个用于扩增子测序数据的全流程分析管线，支持从原始数据到物种注释、丰度分析和多样性分析等多种常见步骤 。该流程集成了质控、去冗余、去嵌合、物种注释、多样性计算和可视化等模块，并且所有步骤都可在 RStudio 中运行  。下文按步骤详细介绍针对第三代全长扩增子（如 Nanopore）测序数据的分析流程。

**EasyAmplicon 2** is a comprehensive analysis pipeline for amplicon sequencing data, supporting common steps from raw data to species annotation, abundance analysis, and diversity analysis. The pipeline integrates modules for quality control, de-redundancy, de-chimerism, species annotation, diversity calculation, and visualization, and all steps can be run in RStudio. The following steps provide a detailed, step-by-step analysis pipeline for third-generation, full-length amplicon sequencing data (e.g., Nanopore sequencing).


## 1. 起始文件 Startup file

1. **工作目录**：在分析前需设置好工作目录和软件/数据库目录（`wd` 和 `db`），并将这些路径加入环境变量。脚本在 Windows Git Bash 下运行，部分命令需切换到 Linux Bash。
(**Working Directory**: Before analysis, set the working directory and software/database directories (`wd` and `db`) and add these paths to the environment variables. This script is run in Windows Git Bash; some commands require switching to Linux Bash.) 

2. **元数据（metadata）**：在 `result/metadata.txt` 中准备样本信息表。表格至少包含三列：第一列为样本 ID（SampleID），最后一列为样本说明（Description）。可以用 `csvtk -t stat result/metadata.txt` 检查行数和列数；并用 `sed` 命令移除 Windows 换行符 `^M` 。   
(**Metadata**: Prepare a table of sample information in `result/metadata.txt`. The table should contain at least three columns: the first column is the sample ID (SampleID), and the last column is the sample description (Description). You can use `csvtk -t stat result/metadata.txt` to check the number of rows and columns; and use the `sed` command to remove the Windows line break characters `^M`.)

3. **测序数据**：将原始 Nanopore FASTQ 数据放到 `seq/` 目录，每个样本对应一个 `.fq` 或 `.fastq` 文件。通常数据以 `.gz` 压缩包形式提供，可用 `gunzip seq/*.gz` 解压（或使用支持直接读取压缩文件的软件）。可用 `seqkit stat seq/*.fastq` 统计每个文件的序列数量和长度等基本信息。
(**Sequencing Data**: Place the raw Nanopore FASTQ data in the `seq/` directory, with each sample corresponding to a `.fq` or `.fastq` file. Data are typically provided as `.gz` compressed files, which can be decompressed using `gunzip seq/*.gz` (or software that supports direct reading of compressed files). `seqkit stat seq/*.fastq` can be used to count basic information such as the number and length of sequences in each file.)

   - 示例：`ls -sh seq/` 查看文件大小；`head seq/E1.fastq | cut -c1-60` 查看序列开头长度。  
   - **可选**：如果要评估数据质量分布，可使用 NanoPlot 等工具绘制读长、质量分布图（该步骤耗时较长）。  
   - Example: `ls -sh seq/` to check the file size; `head seq/E1.fastq | cut -c1-60` to check the length of the beginning of the sequence.
   - **Optional**: To assess data quality distribution, use a tool such as NanoPlot to plot the read length and quality distribution (this step may take some time).

4. **数据库准备**：首次运行时需解压下载好的参考数据库，如 SILVA、GTDB 和 Greengenes。将需要的数据库（例如 `SILVA_modified.fasta`）放到 `${db}/usearch/` 目录，并用 `seqkit stat` 检查序列数量。Greengenes OTU 库也需解压为 `97_otus.fa` 。
(**Database Preparation**: For the first run, you will need to unzip the downloaded reference databases, such as SILVA, GTDB, and Greengenes. Place the desired database (e.g., `SILVA_modified.fasta`) in the `${db}/usearch/` directory and check the sequence counts using `seqkit stat`. The Greengenes OTU library also needs to be unzipped to `97_otus.fa`.)

## 2. 序列合并与重命名 Sequence merging and renaming

对于 Nanopore 全长扩增子数据，每个样本已有一条长读，一般无需合并双端。该步骤主要对序列进行“重命名”和整合：
(For Nanopore full-length amplicon data, each sample already has a long read, so there is generally no need to merge paired ends. This step mainly "renaming" and integrating the sequences:)

- **重命名（Prefix）**：为了追踪来源，每条序列添加样本前缀。例如，使用 `usearch -fastx_relabel seq/E1.fastq -fastqout temp/E1.merged.fastq -prefix E1.` 命令将文件中所有序列 ID 都加上 `E1.` 前缀。可以利用循环批量处理 `metadata.txt` 中列出的所有样本 。  
(- **Rename (Prefix): To track the source, add a sample prefix to each sequence. For example, use `usearch -fastx_relabel seq/E1.fastq -fastqout temp/E1.merged.fastq -prefix E1.` to prefix all sequence IDs in the file with `E1.`. This can be batch processed using a loop for all samples listed in `metadata.txt`.)

- **整合所有样本**：将重命名后的所有样本序列合并到一个文件，`cat temp/*merged.fastq > temp/all.fq`。请注意样本名前缀不能包含“.”，否则后续生成的特征表列名会异常过多，影响分析效率。
(**Merge all samples**: Merge all renamed sample sequences into a single file using `cat temp/*merged.fastq > temp/all.fq`. Note that the sample name prefix cannot contain a dot (.). Otherwise, the generated feature table will have too many columns, affecting analysis efficiency.)  

整合后可用 `head` 检查序列名，确保前缀正确无误。这样每条序列的 ID 形式类似 `E1.read123…`，方便后续区分样本。
(After integration, you can use `head` to check the sequence names to ensure the prefix is ​​correct. This will give each sequence an ID similar to `E1.read123…`, making it easier to distinguish samples later.)

## 3. 切除引物与质量控制 Remove primers and quality control

切除扩增子的引物序列并进行质量过滤，以保证后续分析的准确性。一般使用 `cutadapt` 或 `NanoFilt` 等工具：
(Remove the primer sequence of the amplicon and perform quality filtering to ensure the accuracy of subsequent analysis. Generally, tools such as `cutadapt` or `NanoFilt` are used:)

- **引物去除**：使用 `cutadapt` 命令同时指定正向和反向引物序列（例如 `-g "AGAGTTTGATCCTGGCTCAG...AAGTCSTAACAAGGTADCCSTA"`），`--discard-untrimmed` 丢弃未找到引物的序列，输出到 `temp/allfilter.fq`。运行时可加多线程，如 `-j 20` 以提高速度。  
(- **Primer removal**: Use `cutadapt` to specify both forward and reverse primer sequences (e.g. `-g "AGAGTTTGATCCTGGCTCAG...AAGTCSTAACAAGGTADCCSTA"`), and `--discard-untrimmed` to discard sequences where no primers were found. Output is to `temp/allfilter.fq`. You can use multithreading, such as `-j 20`, to increase speed.)

- **长度和质量过滤**：第三代测序读长较长，常先根据长度范围（例如 1200–1800 bp）进行筛选。方法1：使用 `vsearch --fastx_filter` 设置 `--fastq_minlen 1200 --fastq_maxlen 1800`，同时设定 `--fastq_maxee`（最大期望错误数）和 `--fastq_qmax`（质量编码上限）参数完成质量过滤，输出为 fasta 格式。示例：`vsearch --fastx_filter temp/allfilter.fq --fastq_minlen 1200 --fastq_maxlen 1800 --fastaout temp/filtered.fa --fastq_maxee 30 --fastq_qmax 93`。这里采用较高的最大期望错误（maxEE=30，对应Q≈17）以保留较多有效序列。  
(**Length and quality filtering**: Third-generation sequencing reads are long, so they are often filtered based on a length range (e.g., 1200–1800 bp). Method 1: Use `vsearch --fastx_filter` with `--fastq_minlen 1200 --fastq_maxlen 1800`, along with `--fastq_maxee` (maximum expected error) and `--fastq_qmax` (quality encoding limit) parameters to perform quality filtering and output in fasta format. Example: `vsearch --fastx_filter temp/allfilter.fq --fastq_minlen 1200 --fastq_maxlen 1800 --fastaout temp/filtered.fa --fastq_maxee 30 --fastq_qmax 93`. A higher maximum expected error (maxEE=30, corresponding to Q≈17) is used here to retain more valid sequences.)

- **可选方法2**：使用 `NanoFilt`，例如 `cat temp/allfilter.fq | NanoFilt -l 1200 -q 20 --headcrop 10 --maxlength 1800 > temp/filtered.fastq`，先删去前 10 bp，然后保留长度在 [1200,1800] 且质量≥20 的序列。  
(- **Optional Method 2**: Use `NanoFilt`, for example, `cat temp/allfilter.fq | NanoFilt -l 1200 -q 20 --headcrop 10 --maxlength 1800 > temp/filtered.fastq`, first delete the first 10 bp, and then retain sequences with a length in [1200,1800] and a quality ≥ 20.)

过滤后得到的 `temp/filtered.fa`（或 `.fastq`）为通过质量筛选的序列集，可用 `head` 等命令查看格式是否正确。
(The `temp/filtered.fa` (or `.fastq`) obtained after filtering is the sequence set that passed the quality screening. You can use commands such as `head` to check whether the format is correct.)

## 4. 去冗余、聚类/去噪 De-redundancy, clustering/denoising

### 4.1 序列去冗余（dereplicate）

使用 **VSEARCH** 等工具将完全相同的序列合并为一个，以减少冗余计算：例如 `vsearch --derep_fulllength temp/filtered.fa --fasta_width 0 --sizeout --relabel Uni_ --output temp/uniques.fa --minuniquesize 1`。此命令输出 `uniques.fa`，其中每条序列名类似 `Uni_1;size=5` 表示该序列在所有样本中出现 5 次。第三代长读往往无冗余（size=1），可根据需要调整 `--minuniquesize` 过滤极低丰度序列。
(Use tools like **VSEARCH** to merge identical sequences into a single one to reduce redundant calculations: for example, `vsearch --derep_fulllength temp/filtered.fa --fasta_width 0 --sizeout --relabel Uni_ --output temp/uniques.fa --minuniquesize 1`. This command outputs `uniques.fa`, where each sequence is named like `Uni_1;size=5`, indicating that the sequence occurs five times across all samples. Third-generation long reads are often non-redundant (size=1), so `--minuniquesize` can be adjusted as needed to filter out extremely low-abundance sequences.)

### 4.2 OTU 聚类或 ASV 去噪 OTU clustering or ASV denoising

有两种常见策略：**一）97% OTU 聚类** 或 **二）ASV 去噪**（如 UNOISE3）：
(There are two common strategies: **i) 97% OTU clustering** or **ii) ASV denoising** (such as UNOISE3):)

- **方法1：97% OTU 聚类**（适用于数据量较大或需要传统 OTU 分析时）。可用 USEARCH 的 `-cluster_otus`：
(- **Method 1: 97% OTU clustering** (suitable for large data sets or when traditional OTU analysis is required). Use `-cluster_otus` of USEARCH:)
  ```bash
  usearch -cluster_otus temp/uniques.fa -minsize 1 -otus temp/otus.fa -relabel OTU_
  ```  
  该命令会输出非嵌合体 OTU 表型序列 `otus.fa`，并去除二次及以上冗余低丰度嵌合体。运行后会提示产生的 OTU 数量及去除的嵌合体数。
  (This command will output the non-chimeric OTU phenotype sequence `otus.fa` and remove secondary or higher redundant low-abundance chimeras. After running, the number of OTUs generated and the number of chimeras removed will be displayed.)

- **方法2：UNOISE3 ASV 去噪**（可获得单碱基分辨率的序列变体）。示例：
(- **Method 2: UNOISE3 ASV denoising** (can obtain sequence variants with single-base resolution). Example:)
  ```bash
  usearch -unoise3 temp/uniques.fa -minsize 1 -zotus temp/zotus.fa
  sed 's/Zotu/ASV_/g' temp/zotus.fa > temp/otus.fa
  ```  
  这样生成的 `otus.fa` 中每条序列为一个 ASV（用 `ASV_` 前缀标识）。注意针对第三代数据，高错误率可适当放宽 `-minsize` 要求。因数据量大，UNOISE3 可能耗时较长。
  (Each sequence in the generated `otus.fa` is an ASV (identified by the `ASV_` prefix). Note that for third-generation data, the `-minsize` requirement can be relaxed appropriately if the error rate is high. UNOISE3 may take a long time due to large data volumes.)

- **方法3（可选）**：若数据量过大导致 USEARCH 无法处理，可参考 VSEARCH 的替代方法（脚本 FAQ 部分提供了 vsearch 聚类和去嵌合的示例）。
(- **Method 3 (Optional)**: If the data volume is too large for USEARCH to handle, refer to the alternative method of VSEARCH (the script FAQ section provides examples of vsearch clustering and demosaicing).)

### 4.3 引物和质量过滤的替代：有时可直接使用参考比对去纠错和打磨 Alternatives to primer and quality filtering: Sometimes reference alignment can be used directly for error correction and polishing

**可选：Minimap2 + Racon/Medaka** —— 对 ASV/OTU 进行校正以提高序列准确度。该流程首先用 `minimap2` 将原始读段比对到草图 OTU 上，再用 `racon` 和 `medaka` 对 OTU 序列进行纠错打磨。由于步骤较多、计算量大，需在高性能服务器上进行，最终生成打磨后的 `otus.fa`。此步骤对于大规模研究可选 。
(**Optional: Minimap2 + Racon/Medaka** — Correct ASV/OTU sequences to improve sequence accuracy. This workflow first uses `minimap2` to align raw reads to the draft OTUs, then uses `racon` and `medaka` to correct and polish the OTU sequences. Due to the numerous steps and high computational effort, this process should be performed on a high-performance server. The final result is `otus.fa`. This step is optional for large-scale studies.)

### 4.4 参考库去嵌合（可选）Reference library dechimerization (optional)

默认采用 **de novo** 去嵌合（如上 `cluster_otus` 自带 de novo 去嵌合）。也可进行基于参考数据库的去嵌合，例如用 VSEARCH 的 `--uchime_ref` 对 `otus.fa` 中序列基于 SILVA 库进行查重，输出非嵌合体序列。示例命令：
(By default, de novo dechimerization is used (as described above with `cluster_otus`). Reference-based dechimerization can also be performed, for example, by using `--uchime_ref` with VSEARCH to check the sequences in `otus.fa` against the SILVA library and output non-chimera sequences. Example command:)
```bash
vsearch --uchime_ref temp/otus.fa --db ${db}/usearch/SILVA_modified.fasta --nonchimeras result/raw/otus.fa
```
注：参考库去嵌合可能造成假阴性（高丰度非嵌合序列被保留），仅在需要时使用。否则可直接复制 `temp/otus.fa` 到结果目录。
(Note: Reference library dechimerization may result in false negatives (highly abundant non-chimeric sequences are retained), so use this only when necessary. Otherwise, copy `temp/otus.fa` directly to the result directory.)

## 5. 特征表构建和筛选 Feature table construction and screening

OTU/ASV 序列被视为特征（Feature），需要构建每个样本中这些特征的丰度表。
(OTU/ASV sequences are considered as features, and the abundance table of these features in each sample needs to be constructed.)

### 5.1 生成特征表 Generate feature table

使用 **USEARCH** 或 **VSEARCH** 将每条质控后序列比对到 OTU/ASV 序列，统计每个 OTU 在每个样本中的读数。方法示例：
(Use **USEARCH** or **VSEARCH** to align each quality control sequence to the OTU/ASV sequence and count the number of reads for each OTU in each sample. Example method:)
- **VSEARCH method**：
  ```bash
  vsearch --usearch_global temp/filtered.fa --db result/raw/otus.fa --id 0.97 --otutabout result/raw/otutab.txt --threads 16
  ```  
  该命令以 97% 相似度比对所有序列到 OTU，生成 OTU 样本矩阵 `otutab.txt`。执行后用 `csvtk -t stat` 检查特征表行列数，列数应与样本数一致。
  (This command aligns all sequences to OTUs at 97% similarity, generating the OTU sample matrix `otutab.txt`. After executing this command, use `csvtk -t stat` to check the number of rows and columns in the feature table; the number of columns should match the number of samples.)  

- （可选）**USEARCH 方法**：速度稍快但多线程效率有限，适合样本数较少时使用。
(- (Optional) **USEARCH method**: Slightly faster but with limited multi-threading efficiency, suitable for use when the number of samples is small.)

需要注意 Windows 下的换行符问题，可用 `sed -i 's/\r//g'` 将生成的结果转为标准 Unix 格式。然后通过 `head` 命令预览表格前几行，确认结果正确。
(Be careful about line breaks on Windows. Use `sed -i 's/\r//g'` to convert the output to standard Unix format. Then use the `head` command to preview the first few rows of the table to confirm the output is correct.)

### 5.2 物种注释与非细菌/真核过滤 Species annotation and non-bacterial/eukaryotic filtering

- **物种注释**：用 `vsearch --sintax` 对最终的 OTU/ASV 序列文件进行物种分类注释。常用 SILVA 数据库，如 `SILVA_modified.fasta`（已包含物种级别标签）。例如：
(- **Species Annotation**: Use `vsearch --sintax` to annotate the final OTU/ASV sequence files with species taxonomy. The SILVA database is commonly used, such as `SILVA_modified.fasta` (which already contains species-level labels). For example:)
  ```bash
  vsearch --sintax result/raw/otus.fa --db ${db}/usearch/SILVA_modified.fasta --sintax_cutoff 0.1 --tabbedout result/raw/otus.sintax
  ```  
  结果 `otus.sintax` 包含每个 OTU 的分类注释和置信度，可用 `sed` 处理移除多余符号。
  (The resulting `otus.sintax` file contains the taxonomic annotation and confidence for each OTU, and can be processed with `sed` to remove redundant symbols.)

- **去除叶绿体、线粒体和非细菌**：通常使用 R 脚本（如 `otutab_filter_nonBac.R`）将特征表中过滤掉真核、叶绿体、线粒体和古菌等分类（保留细菌/古菌），并输出新的特征表 `result/otutab.txt` 和统计文件 `result/raw/otutab_nonBac.stat`。如果认为不需要过滤，也可跳过此步直接复制 `result/raw/otutab.txt` 作为输出结果。 
(**Remove Chloroplasts, Mitochondria, and Non-Bacteria**: Typically, an R script (e.g., `otutab_filter_nonBac.R`) is used to filter out Eukarya, Chloroplasts, Mitochondria, and Archaea from the feature table (retaining Bacteria/Archaea). This script then outputs a new feature table `result/otutab.txt` and a statistics file `result/raw/otutab_nonBac.stat`. If filtering is not necessary, skip this step and copy `result/raw/otutab.txt` directly as the output.)
 
- **物种注释整理**：将 sintax 注释提取为两列（OTU ID 和分类），或展开为 8 列（域、门、纲、目、科、属、种）输出为 `result/taxonomy.txt`。未注释的层级用 “Unassigned” 填充。
(- **Species annotation cleanup**: Extract sintax annotations into two columns (OTU ID and taxonomy), or expand them into eight columns (domain, phylum, class, order, family, genus, species) and output them as `result/taxonomy.txt`. Unannotated levels are filled with "Unassigned".)

- **统计特征表信息**：可用 `usearch -otutab_stats result/otutab.txt` 简要统计每个样本的序列数、均值、最小最大值等，检查测序深度分布，为下一步抽样提供参考。
(- **Statistical feature table information**: You can use `usearch -otutab_stats result/otutab.txt` to briefly count the number of sequences, mean, minimum and maximum values ​​of each sample, check the sequencing depth distribution, and provide a reference for the next sampling step.)

### 5.3 等量抽样标准化 Equal sampling standardization

对样本间测序深度差异进行标准化。使用 R 的 vegan 包脚本对 `result/otutab.txt` 进行重抽样稀释，使所有样本统一抽样深度（例如 139 reads），输出稀释后的特征表 `result/otutab_rare.txt` 和多样性指标 `result/alpha/vegan.txt`。然后再次用 `usearch -otutab_stats` 检查稀释后表的统计信息。
(Normalize for differences in sequencing depth between samples. Use the R package vegan script to resample `result/otutab.txt` to a uniform sampling depth for all samples (e.g., 139 reads). Output the resampled feature table `result/otutab_rare.txt` and the diversity metric `result/alpha/vegan.txt`. Then, use `usearch -otutab_stats` again to examine the statistics of the resampled table.)

## 6. α 多样性 α diversity

计算 α 多样性指数：使用 USEARCH 的 -alpha_div 命令计算诸如 Shannon、Simpson、ACE 等 14 种多样性指数（但 Chao1 在此工具中计算有误，建议不用）。例：
(Calculate alpha diversity index: Use the -alpha_div command of USEARCH to calculate 14 diversity indices such as Shannon, Simpson, ACE, etc. (However, Chao1 is calculated incorrectly in this tool and is not recommended). Example:)

```bash
usearch -alpha_div result/otutab_rare.txt -output result/alpha/alpha.txt
```

输出文件包含各样本的多样性指数值。
(The output file contains the diversity index value of each sample.)

稀释曲线：使用 usearch -alpha_div_rare result/otutab_rare.txt -output result/alpha/alpha_rare.txt 生成从 1% 到 100% 的采样深度对应的 OTU 丰富度。可处理输出中的 - 异常值为 0。
(Rarefaction curve: Generate OTU richness corresponding to sampling depth from 1% to 100% using usearch -alpha_div_rare result/otutab_rare.txt -output result/alpha/alpha_rare.txt. Outliers in the output are handled as 0.)

高丰度特征筛选：可根据分组计算特征在各组中的平均丰度，然后筛选出高丰度特征（例如平均丰度 > 0.1%），生成各组共有和特有 OTU 列表。也可使用在线工具绘制 Venn 图或热图展示分组特征。
(High-abundance feature screening: Calculate the average abundance of features within each group, then screen for high-abundance features (e.g., average abundance > 0.1%) to generate lists of shared and unique OTUs for each group. Alternatively, use online tools to create Venn diagrams or heat maps to visualize grouped features.)

## 7. β 多样性 β diversity

计算样本间差异和可视化所需的距离矩阵及树构建：
(Calculate the distance matrix and tree construction required for sample differences and visualization:)

系统发育树：使用 usearch -cluster_agg 对 result/otus.fa 序列进行对齐聚类，生成简易进化树 result/otus.tree（用于基于进化距离的多样性分析）。
(Phylogenetic tree: Use usearch -cluster_agg to align and cluster the result/otus.fa sequences to generate a simple phylogenetic tree result/otus.tree (for diversity analysis based on evolutionary distance).)

距离矩阵：使用 usearch -beta_div 基于稀释后特征表 (otutab_rare.txt) 生成多种样本间距离矩阵（如 Bray-Curtis、Jaccard、Unifrac 等），文件前缀为 result/beta/。这些可用于后续的PCoA/NMDS分析。
(Distance matrices: Use usearch -beta_div to generate various inter-sample distance matrices (e.g., Bray-Curtis, Jaccard, Unifrac, etc.) based on the rarefaction feature table (otutab_rare.txt). The files are prefixed with result/beta/. These can be used for subsequent PCoA/NMDS analysis.)

## 8. 注释汇总 Annotation summary

输出 OTU 与物种的对应关系：
(Output the correspondence between OTU and species:)

将 otus.sintax 中的分类注释整理为两列（OTUID 和 Species），或展开为多列（门、纲、目、科、属、种），并输出为 result/taxonomy.txt 便于查看。
(Arrange the taxonomic annotations in otus.sintax into two columns (OTUID and Species), or expand them into multiple columns (phylum, class, order, family, genus, species), and output them as result/taxonomy.txt for easy viewing.)

用 usearch -sintax_summary 对每个分类等级（phylum、class、order、family、genus、species）统计各样本的丰度分布，输出汇总表 result/tax/sum_*。
(Use usearch -sintax_summary to calculate the abundance distribution of each sample at each taxonomic level (phylum, class, order, family, genus, species) and output the summary table result/tax/sum_*.)

这些结果可用于绘制分类丰度图（柱状图、饼图等）。
(These results can be used to plot taxonomic abundances (histograms, pie charts, etc.).)

## 9. 参考数据库定量表 Reference database quantitative table

为功能预测准备参考OTU表：将序列比对到 Greengenes 97% OTUs 库，生成用于 PICRUSt 等工具的 OTU 表。示例命令：
(Prepare a reference OTU table for function prediction: Align the sequence to the Greengenes 97% OTUs library to generate an OTU table for tools such as PICRUSt. Example command:)

```bash
vsearch --usearch_global temp/filtered.fa --db ${db}/gg/97_otus.fa --otutabout result/gg/otutab.txt --id 0.97 --threads 12
```

该表描述了输入序列与 Greengenes OTU 的映射关系，可用 usearch -otutab_stats 查看统计信息。
(This table describes the mapping relationship between the input sequence and the Greengenes OTU. You can use usearch -otutab_stats to view the statistical information.)

## 10. 清理和提交 Clean and commit

删除中间大文件（如 temp/*.fq）释放磁盘空间。
(Delete large intermediate files (such as temp/*.fq) to free up disk space.)

计算原始序列文件的 MD5 值（md5sum seq/*.fastq > result/md5sum.txt）以便提交时校验数据完整性。
(Calculate the MD5 value of the original sequence file (md5sum seq/*.fastq > result/md5sum.txt) to verify data integrity when submitting.)

### Emu 软件流程（可选）Emu software pipeline (optional)

除了传统的 OTU/ASV 分析，EasyAmplicon2 还支持 Emu 软件的快速物种注释流程。Emu 通过将全长16S序列直接比对到参考数据库（如 SILVA、GTDB）来生成物种丰度谱，适用于需要快速物种组成分析的场景pubmed.ncbi.nlm.nih.gov。主要步骤包括：
(In addition to traditional OTU/ASV analysis, EasyAmplicon2 also supports the rapid species annotation workflow of Emu software. Emu generates species abundance profiles by directly aligning full-length 16S sequences to reference databases (such as SILVA and GTDB), making it suitable for scenarios requiring rapid species composition analysis (pubmed.ncbi.nlm.nih.gov). The main steps include:)

Reads 重标记：为每个样本的原始 FASTQ 文件加样本前缀（同样使用 usearch -fastx_relabel），以便混合后仍可区分来源。
(Read relabeling: Add a sample prefix to the original FASTQ file of each sample (also using usearch -fastx_relabel) so that the source can still be distinguished after mixing.)

单样本质控：分别对每个样本 FASTQ 文件进行引物修剪和质量过滤（使用 cutadapt 和 NanoFilt），输出每个样本独立的过滤后文件。
(Single-sample quality control: Perform primer trimming and quality filtering (using cutadapt and NanoFilt) on each sample FASTQ file separately, and output a separate filtered file for each sample.)

运行 Emu 注释：对每个样本运行 emu abundance，分别使用 SILVA 和 GTDB 数据库（需要预先构建好的 Emu 数据库）。该命令会输出每个样本的多级分类丰度表（species、genus、family 等）。
(Run Emu annotation: Run emu abundance on each sample, using the SILVA and GTDB databases (pre-built Emu databases are required). This command will output a multi-level taxonomic abundance table (species, genus, family, etc.) for each sample.)

结果合并：运行 emu combine-outputs 命令在不同分类等级上合并所有样本结果，生成跨样本的丰度矩阵，方便后续统计和可视化。
(Merge results: Run the emu combine-outputs command to merge all sample results at different taxonomic levels to generate a cross-sample abundance matrix for subsequent statistics and visualization.)

Emu 方法省去了 OTU 构建过程，直接得到定量的物种水平丰度。原理上，Emu 采用一种期望最大化算法处理全长长读数据，以获得更准确的物种丰度预测。
(The Emu method eliminates the OTU construction process and directly obtains quantitative species-level abundance. In principle, Emu uses an expectation-maximization algorithm to process full-length long read data to obtain more accurate species abundance predictions.)

References:

使用此脚本，请引用下文：

If used this script, please cited:

Yong-Xin Liu, Lei Chen, Tengfei Ma, Xiaofang Li, Maosheng Zheng, Xin Zhou, Liang Chen, Xubo Qian, Jiao Xi, Hongye Lu, Huiluo Cao, Xiaoya Ma, Bian Bian, Pengfan Zhang, Jiqiu Wu, Ren-You Gan, Baolei Jia, Linyang Sun, Zhicheng Ju, Yunyun Gao, Tao Wen, Tong Chen. 2023. EasyAmplicon: An easy-to-use, open-source, reproducible, and community-based pipeline for amplicon data analysis in microbiome research. iMeta 2: e83. https://doi.org/10.1002/imt2.83