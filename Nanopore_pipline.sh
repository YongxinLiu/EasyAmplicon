
[TOC]

# EasyAmplicon2 Nanopore Pipeline
# EasyAmplicon2 Nanopore 分析流程

    # 作者 Author: 刘永鑫（Yong-xin Liu）等
    # 更新时间 Update: 2025-07-24
    # 版本 Version: 2.1

    # 设置工作(work directory, wd)和软件数据库(database, db)目录
    # Set the working directory (wd) and the software/database directory (db)
    # 流程大部分在Git bash下运行，需要转换linux bash时会注明
    # Most of the pipeline runs in Git Bash. A note will be made when it is necessary to switch to Linux Bash.
    
    # 添加环境变量，并进入工作目录
    # Add environmental variables and enter the working directory
    # **每次打开Rstudio必须运行下面4行 Run it**，可选替换${db}为EasyMicrobiome安装位置
    # **The following 4 lines must be run every time you open RStudio**, you can optionally replace ${db} with your EasyMicrobiome installation location
    wd=/c/EasyAmplicon2
    db=/c/EasyMicrobiome
    PATH=$PATH:${db}/win
    cd ${wd}



## 1. 起始文件
## 1. Start files

    # 1. 分析流程pipeline.sh
    # 1. Analysis pipeline: pipeline.sh
    # 2. 样本元信息metadata.txt，保存于result目录
    # 2. Sample metadata: metadata.txt, saved in the result directory
    # 3. 测序数据fastq文件保存于seq目录，通常以`.fq.gz`结尾，每个样品一对文件
    # 3. Sequencing data: FASTQ files are saved in the seq directory, usually ending with `.fq.gz`, one file per sample.
    # 4. 进入对应工作目录，创建临时文件存储目录，分析结束可删除
    # 4. Enter the corresponding working directory, create a temporary directory for storing intermediate files, which can be deleted after the analysis is complete.
    cd Nanopore
    mkdir -p seq result temp 

### 1.1. 元数据/实验设计
### 1.1. Metadata/Experimental Design

    # 准备样本元数据result/metadata.txt
    # Prepare the sample metadata file: result/metadata.txt
    # csvtk统计表行(样本数，不含表头)列数，-t设置列分隔为制表符，默认为;
    # Use csvtk to count the rows (number of samples, excluding header) and columns of the table. -t sets the column separator to tab, default is comma.
    csvtk -t stat result/metadata.txt
    # 元数据至少3列，首列为样本ID(SampleID)，结尾列为描述(Description)
    # The metadata must have at least 3 columns: the first column is the sample ID (SampleID), and the last column is the description (Description).
    # cat查看文件，-A显示符号，"|"为管道符实现命令连用，head显示文件头，-n3控制范围前3行
    # Use 'cat -A' to view the file with non-printing characters. '|' is a pipe to chain commands. 'head -n3' displays the first 3 lines.
    cat -A result/metadata.txt | head -n3
    # windows用户结尾有^M，运行sed命令去除，再用cat -A检查
    # For Windows users, if there is a ^M at the end of the line, run the sed command to remove it, and then check with 'cat -A'.
    sed 's/\r//' result/metadata.txt
    cat -A result/metadata.txt | head -n3

### 1.2. 测序数据
### 1.2. Sequencing data
    # Nanopore 原始数据处理（以 FAST5 格式为例，可选）
    # Nanopore Raw Data Processing (FAST5 Example, Optional)
    ---
    ## Guppy 软件的下载安装
    ## Guppy Installation
    ### 第一步：下载 Guppy
    从 Oxford Nanopore 官网下载安装包：  
    Download the Guppy installer from the official Oxford Nanopore website:  
    [https://cdn.oxfordnanoportal.com/software/analysis/ont-guppy-cpu-6.5.7-win64.msi](https://cdn.oxfordnanoportal.com/software/analysis/ont-guppy-cpu-6.5.7-win64.msi)
    ---
    ### 第二步：安装 Guppy
    双击 `.msi` 文件并根据安装指引完成安装。安装完成后，将安装目录替换为你的实际路径，例如 `/d/guppy`：
    Double-click the .msi file and follow the installation instructions.
    After installation, replace /d/guppy with your actual installation path:
    ```bash
    mkdir -p /d/guppy
    第三步：处理 FAST5 格式文件
    使用 Guppy 对 FAST5 原始数据进行 basecalling，生成 .fastq 文件：
    Step 3: Process FAST5 Files
    Use Guppy to perform basecalling on FAST5 raw data, which will generate .fastq files:
    /d/guppy/bin/guppy_basecaller \
    -i fast5/ \
    -s output \
    --config dna_r9.4.1_450bps_fast.cfg \
    -r \
    --num_callers 4 \
    --cpu_threads_per_caller 2
    输出结果
    运行完成后，生成的 .fastq 文件将保存在 output/ 文件夹中。
    Output
    The resulting .fastq files will be saved in the output/ directory.
    ```

    # # 本段代码可在RStudio中Ctrl + Shift + C 取消注释“#”后运行
    # # This code block can be run in RStudio after uncommenting with Ctrl + Shift + C
    # # 测序数据下载百度网盘链接：https://pan.baidu.com/s/1Ikd_47HHODOqC3Rcx6eJ6Q?pwd=0315；文件路径：db/amplicon/Nanopore/seq
    # # Baidu Net Disk link for sequencing data download: https://pan.baidu.com/s/1Ikd_47HHODOqC3Rcx6eJ6Q?pwd=0315; File path: db/amplicon/Nanopore/seq
    # 公司返回的三代（Pacbio或者Nanopore）全长扩增子测序结果，通常为一个样品只有一个fq/fastq.gz格式压缩文件
    # The full-length amplicon sequencing results from third-generation sequencing (PacBio or Nanopore) usually consist of one compressed fq/fastq.gz file per sample.
    # 如果测序数据是.gz的压缩文件，有时需要使用gunzip解压后使用，vsearch通常可直接读取压缩文件
    # If the sequencing data is in a .gz compressed file, it sometimes needs to be decompressed with gunzip before use. Vsearch can usually read compressed files directly.
    gunzip seq/*.gz
    # zless按页查看压缩文件，空格翻页、q退出；head默认查看前10行，-n指定行
    # Use 'zless' to view compressed files page by page (space to scroll, q to quit); 'head' displays the first 10 lines by default, -n specifies the number of lines.
    ls -sh seq/
    head seq/E1.fastq
    # 每行太长，指定查看每行的1-60个字符
    # The lines are too long, so we view only the first 60 characters of each line.
    head seq/E1.fastq | cut -c 1-60
    # 统计测序数据，依赖seqkit程序
    # Use seqkit to get statistics of the sequencing data.
    seqkit stat seq/E1.fastq
    # 批量统计测序数据并汇总表
    # Batch process all sequencing files and summarize the statistics in a table.
    seqkit stat seq/*.fastq > result/seqkit.txt
    head result/seqkit.txt

    # ## (Optional)In order to perform quality filtering,  first assess sequence length and quality with nanoplot
    # ## （可选）为了进行质量过滤，首先使用nanoplot评估序列长度和质量，耗时较长
    # # mkdir -p nanoplot_reports  
    # ## 切换为linux bash，首次打开需要输入：bash
    # ## Switch to linux bash, type 'bash' when opening for the first time
    # conda activate easyamplicon2 
    # # Set directory paths (adjusted for your environment)
    # # 设置目录路径（根据您的环境进行调整）
    # INPUT_DIR="seq"
    # OUTPUT_DIR="result/nanoplot_reports"
    # 
    # # Create output directory
    # # 创建输出目录
    # mkdir -p "$OUTPUT_DIR"
    # 
    # # List of  FASTQ files
    # # FASTQ 文件列表
    # SAMPLES=("E1" "E2" "E3" "E4" "E5" "E6" "E7" "N1" "N2" "O1" "O2" "S1")
    # 
    # # Run NanoPlot
    # # 运行 NanoPlot
    # for SAMPLE in "${SAMPLES[@]}"; do
    #   echo "Processing $SAMPLE..."
    #   NanoPlot --fastq "$INPUT_DIR/${SAMPLE}.fastq" \
    #   --outdir "$OUTPUT_DIR/${SAMPLE}_report/" \
    #   --plots hex dot --loglength --N50
    # done

### 1.3. 流程和数据库
### 1.3. Pipeline & Database

    # 数据库第一次使用必须解压，以后可跳过此段
    # The database must be decompressed before the first use. This step can be skipped later.

    # usearch可用的16S/18S/ITS数据库：
    # Available 16S/18S/ITS databases for usearch:

    # 解压Silva数据库，需自行从官网或百度网盘 (https://pan.baidu.com/s/1Ikd_47HHODOqC3Rcx6eJ6Q?pwd=0315; db/amplicon/silva) 下载
    # Decompress the Silva database. You need to download it from the official website or the Baidu Net Disk (https://pan.baidu.com/s/1Ikd_47HHODOqC3Rcx6eJ6Q?pwd=0315; db/amplicon/silva).

    # 将SILVA数据库SILVA_modified.fasta保存在${db}/usearch/
    # Save SILVA_modified.fasta to ${db}/usearch/

    # 将SILVA数据库silva_nr99_v138.1_train_DADA2.fa.gz保存在${db}/DADA2/
    # Save silva_nr99_v138.1_train_DADA2.fa.gz to ${db}/DADA2/

    # 将其他可选数据库从百度网盘 (https://pan.baidu.com/s/1Ikd_47HHODOqC3Rcx6eJ6Q?pwd=0315; db/amplicon/usearch)保存到${db}/usearch/目录，例如：
    # Save other optional databases to ${db}/usearch/, for example:
    #   sintax_defalut_emu_database.fasta.gz
    #   sintax_ncbi_database.fasta.gz
    #   gtdb_sintax_database.fasta.gz
    seqkit stat ${db}/usearch/SILVA_modified.fasta
    # 将Silva_Emu/和GTDB_Emu/文件夹保存到{db}
    # Save Silva_Emu/、GTDB_Emu/ to {db}
    # gunzip -c ${db}/usearch/utax_reference_dataset_all_29.11.2022.fasta.gz >${db}/usearch/unite.fa
    # seqkit stat ${db}/usearch/unite.fa # 32.6万 / 326k sequences
    # Greengene数据库用于功能注释: ftp://greengenes.microbio.me/greengenes_release/gg_13_5/gg_13_8_otus.tar.gz
    # Greengenes database for functional annotation: ftp://greengenes.microbio.me/greengenes_release/gg_13_5/gg_13_8_otus.tar.gz
    # 默认解压会删除原文件，-c指定输出至屏幕，> 写入新文件(可改名)
    # Decompression deletes the original file by default. -c specifies output to the screen, > writes to a new file (can be renamed).
    gunzip -c ${db}/gg/97_otus.fasta.gz > ${db}/gg/97_otus.fa
    seqkit stat ${db}/gg/97_otus.fa


## 2. 序列合并和重命名
## 2. Reads merge and rename

      
### 2.1 文件改名
### 2.1 Reads rename

    # # 单个序列改名示例
    # # Example of renaming a single sequence file
    # i=E1
    # gunzip -c seq/${i}.fq.gz > seq/${i}.fq
    # usearch -fastx_relabel seq/${i}.fq -fastqout temp/${i}.merged.fq -prefix ${i}.

    # # 批量改名，解压(usearch不支持压缩格式)
    # # Batch rename, decompress first (usearch does not support compressed format)
    # gunzip seq/*.gz
        
    # usearch如果出现不兼容系统的情况，请解压${db}/win/usearch11.0.667_win32.exe.gz使用或官网https://www.drive5.com/usearch/download.html下载最新版本解压使用
    # If usearch is not compatible with your system, please decompress ${db}/win/usearch11.0.667_win32.exe.gz for use, or download the latest version from the official website: https://www.drive5.com/usearch/download.html and decompress it for use.
    time for i in `tail -n+2 result/metadata.txt|cut -f1`;do
     usearch -fastx_relabel seq/${i}.fastq -fastqout temp/${i}.merged.fastq -prefix ${i}.
     done &
    # # vsearch大数据方法参考“常见问题2”
    # # For large datasets, refer to "FAQ 2" for the vsearch method.

### 2.2 改名后序列整合
### 2.2 Integrate renamed reads

    #合并所有样品至同一文件
    # Merge all samples into a single file
    cat temp/*merged.fastq > temp/all.fq
    #查看文件大小8.4G，软件不同版本结果略有差异
    # Check the file size (e.g., 8.4G). Results may vary slightly with different software versions.
    ls -lsh temp/all.fq
    # 查看序列名，“.”之前是否为样本名，样本名绝不允许有点 (".")
    # Check the sequence names. The part before the "." should be the sample name. Sample names must not contain dots (".").
    # 样本名有点 (.) 的一个显著特征是生成的特征表会很大，特征表里面列很多，导致后面分析出现内存不足。
    # A significant feature of sample names with dots is that the generated feature table will be very large, with many columns, leading to memory shortages in subsequent analysis.
    # 后面分析获得特征表后要看一眼有没有问题，遇到内存不足问题，也要回头来排查。
    # After obtaining the feature table, you should check it for any issues. If you encounter memory problems, you should go back and investigate.
    head -n 6 temp/all.fq|cut -c1-60


## 3. 切除引物与质控
## 3. Cut primers and quality filter

    # 注意将gitbash切换为linux bash
    # Note: Switch from Git Bash to Linux Bash
    # conda activate easyAmplicon2
    # 过滤序列，去除接头："前向引物...反向引物",10s
    # Filter sequences, remove adapters: "forward_primer...reverse_primer", takes about 10s
    # cutadapt -g "AGAGTTTGATCCTGGCTCAG...TACGGYTACCTTGTTACGACTT" \
    # --error-rate=0.1 \
    # -j 10 \
    # --discard-untrimmed \
    # -o temp/allfilter.fastq temp/all.fq
    cutadapt -g "AGAGTTTGATCCTGGCTCAG...AAGTCSTAACAAGGTADCCSTA" \
    --error-rate=0.1 \
    --action=trim \
    --rc \
    -j 20 \
    --discard-untrimmed \
    -o temp/allfilter.fq \
    temp/all.fq
    
    # 方法1 Quality filtering with Vsearch
    # Method 1: Quality filtering with Vsearch

    # 细菌16S片段长度通常约为1500bp，为避免过长或过短序列干扰，使用vsearch进行片段长度筛选。筛选长度在1200~1800 bp之间的序列，输出为fasta格式
    # The length of bacterial 16S fragments is usually around 1500 bp. To avoid interference from excessively long or short sequences, use vsearch for fragment length filtering. Filter sequences between 1200 and 1800 bp and output in FASTA format.
    # 由于三代数据质量较高故设置--fastq_qmax 93参数（默认41，二代最高不超过40，三代最高不超过93），1m9s
    # Since the quality of third-generation data is high, set the --fastq_qmax parameter to 93 (default is 41; second-generation does not exceed 40, third-generation does not exceed 93). Takes about 1m 9s.
    # 在 .fastq 文件中，Q 值不是直接写数字，而是通过 ASCII 字符Phred+33编码
    # In a .fastq file, the Q score is not written as a number but is encoded by ASCII characters (Phred+33).
    # --fastq_maxee 控制每条序列的最大期望错误数，Q 值越高、允许的 maxEE 越低。对于二代Illumina数据常设为maxEE=1~2（对应Q≈30），而三代（如ONT/PacBio）因错误率高，常设为maxEE=50~100（对应Q≈12~15），以避免过度过滤。这里我们预设30(Q17)产出21781条序列，3m
    # --fastq_maxee controls the maximum expected error per sequence. The higher the Q score, the lower the allowed maxEE. For second-generation Illumina data, maxEE is often set to 1-2 (corresponding to Q≈30), while for third-generation (e.g., ONT/PacBio), due to high error rates, it is often set to maxEE=50-100 (corresponding to Q≈12-15) to avoid over-filtering. Here we preset 30 (Q17), which yields 21,781 sequences in 3 minutes.
    vsearch --fastx_filter temp/allfilter.fq --fastq_minlen 1200 --fastq_maxlen 1800 --fastaout temp/filtered.fa --fastq_maxee 30 --fastq_qmax 93
    # 另外为了后续筛序合适数量的otu序列，我们可以设置更严格的maxEE
    # Additionally, to screen for an appropriate number of OTU sequences later, we can set a stricter maxEE.
    # 查看文件了解fa文件格式,并对比质控前后变化
    # View the file to understand the FASTA format and compare the changes before and after quality control.
    head temp/filtered.fa

    # 方法2 Quality filtering with NanoFilt
    # Method 2: Quality filtering with NanoFilt
    # # 注意将gitbash切换为linux bash
    # # Note: Switch from Git Bash to Linux Bash
    # conda activate easyAmplicon2
    # # 对 allfilter.fq 文件进行筛选，保留长度 1200-1800bp、质量 Q>=20、去掉前10bp，输出为 filtered.fa
    # # Filter the allfilter.fq file, keeping sequences with a length of 1200-1800 bp, quality Q>=20, and removing the first 10 bp. Output to filtered.fa.
    # cat temp/allfilter.fq | NanoFilt -l 1200 -q 20 --headcrop 10 --maxlength 1800 > temp/filtered.fa
    # # 如果数据为双链，输出的序列头可能会带有 rc, 表示该序列是反向互补(reverse complement)，用以下命令删除
    # # If the data is double-stranded, the output sequence header may have "rc", indicating that the sequence is a reverse complement. Use the following command to remove it.
    # # sed 's/ rc$//' > temp/filtered.fa
    # ## For emu pipline: we should have Quality filtering with NanoFilt  and keep FASTQ
    # ## 对于emu流程：我们应该使用NanoFilt进行质量过滤并保留FASTQ格式
    # cat temp/allfilter.fq | NanoFilt -l 1200 -q 20 --headcrop 10 --maxlength 1800 > temp/filtered.fastq
    # 
    # # 查看文件了解fa文件格式,并对比质控前后变化
    # # View the file to understand the FASTA format and compare the changes before and after quality control.
    # head temp/filtered.fa
    
## 4. 去冗余挑选OTU/ASV
## 4. Dereplicate and cluster/denoise

### 4.1 序列去冗余
### 4.1 Dereplicate sequences

    #去冗余：完全一致的序列合并，统计size信息，输出唯一序列，-minsize表示保留序列的最小条数,-sizeout输出丰度,--fasta_width 0输出FASTA时不换行（单行序列）减少文件体积 --relabel必须加序列前缀更规范, 4s
    # Dereplicate: Merge identical sequences, count size information, output unique sequences. -minsize specifies the minimum number of sequences to keep. -sizeout outputs abundance. --fasta_width 0 outputs FASTA without line breaks (single-line sequence) to reduce file size. --relabel must be added with a sequence prefix for better standardization. Takes 4s.
    time vsearch --derep_fulllength temp/filtered.fa --fasta_width 0 --sizeout --relabel Uni_  --output temp/uniques.fa --minuniquesize 1 --threads 8
     
    #查看uniques.fa文件
    # View the uniques.fa file
    ls -lsh temp/uniques.fa
    # Uni_1;size=1  - 去冗余后序列的名字 Uni_1；该序列在所有样品测序数据中出现 1 次，三代长读长测序数据常会出现无冗余情况，用户根据自己的数据情况调整
    # Uni_1;size=1 - The name of the sequence after dereplication is Uni_1; this sequence appears once in all sample sequencing data. Third-generation long-read sequencing data often has no redundancy. Users should adjust according to their own data.

### 4.2 聚类OTU/去噪ASV
### 4.2 Cluster OTUs / Denoise ASVs

    #有两种方法：推荐unoise3去噪获得单碱基精度ASV，备选传统的97%聚类OTU(属水平精度)
    # There are two methods: unoise3 denoising to obtain single-base precision ASVs is recommended. The alternative is traditional 97% clustering for OTUs (genus-level precision).
    #usearch两种特征挑选方法均自带de novo去嵌合体
    # Both feature selection methods in usearch include de novo chimera removal.
    #-minsize二次过滤，这里因为数据size都是1故选方法1结合质控中maxEE参数，来控制OTU/ASV数量至1-5千，方便下游统计分析
    # -minsize for secondary filtering. Here, because the data size is all 1, method 1 is chosen in combination with the maxEE parameter in quality control to control the number of OTUs/ASVs to 1-5 thousand, which is convenient for downstream statistical analysis.

    #方法1. 97%聚类OTU，适合大数据/ASV规律不明显/reviewer要求
    # Method 1. 97% OTU clustering, suitable for big data / when ASV patterns are not obvious / required by reviewers.
    #结果耗时2m1s, 产生 4402 OTUs, 去除348 chimeras
    # Results took 2m 1s, produced 4402 OTUs, and removed 348 chimeras.
     usearch -cluster_otus temp/uniques.fa -minsize 1 \
      -otus temp/otus.fa \
      -relabel OTU_ 

    #方法2. ASV去噪 Denoise: predict biological sequences and filter chimeras
    # Method 2. ASV Denoise: predict biological sequences and filter chimeras
    #42s, 21783 good, 0 chimeras, 序列百万条可能需要几天/几周
    # 42s, 21,783 good, 0 chimeras. Millions of sequences may take days/weeks.
    # usearch -unoise3 temp/uniques.fa -minsize 1 \
    #   -zotus temp/zotus.fa
    ##三代测序数据size数会出现普遍较小的情况，根据实际情况考虑是否调整 -minsize 
    ## The size of third-generation sequencing data will generally be smaller. Consider whether to adjust -minsize according to the actual situation.
    #修改序列名：Zotu改为ASV方便识别
    # Rename sequences: Change Zotu to ASV for easy identification.
    # sed 's/Zotu/ASV_/g' temp/zotus.fa > temp/otus.fa
    head -n 2 temp/otus.fa

    #方法3. 数据过大无法使用usearch时，备选vsearch方法见"常见问题1"
    # Method 3. When the data is too large to use usearch, see "FAQ 1" for the alternative vsearch method.

### 4.3 （可选）Aligment with minimap2 and plosihing with meadaka 追求更精确的结果可选，PC耗时达好几天，建议高性能服务器尝试
### 4.3 (Optional) Alignment with minimap2 and polishing with Medaka. Optional for more accurate results. It can take several days on a PC, so it is recommended to try it on a high-performance server.

#     # 注意切换linux bash
#     # Note: Switch to Linux Bash
#     conda activate easyAmplicon2
# 
# #### Step1 : Alignment with minimap2: 
# #### 步骤1：使用minimap2进行比对
# 
#     # Set working directory
#     # 设置工作目录
#     cd /mnt/c/EasyAmplicon2/Nanopore
#     # Align reads to draft OTUs
#     # 将reads比对到草稿OTUs
#     minimap2 -ax map-ont -t 32 temp/otus.fa temp/all.fq > temp/aligned_output.sam
#     # Convert SAM to BAM, sort, and index
#     # 将SAM转换为BAM，排序并建立索引
#     samtools view -Sb temp/aligned_output.sam > temp/aligned_output.bam
#     samtools sort temp/aligned_output.bam -o temp/aligned_output_sorted.bam
#     samtools index temp/aligned_output_sorted.bam
# 
# #### Step 2: RACON Polishing
# #### 步骤2：RACON打磨
# 
#     # Apply RACON to polish OTUs
#     # 使用RACON打磨OTUs
#     racon -m 8 -x -6 -g -8 -w 500 -t 32 temp/all.fq temp/aligned_output.sam temp/otus.fa > temp/otus_racon.fa
#     ##Step 3: Prepare for MEDAKA
#     ##步骤3：准备MEDAKA
#     # Index FASTA (optional for medaka, helpful for inspection)
#     # 索引FASTA（medaka可选，有助于检查）
#     samtools faidx temp/otus_racon.fa
# 
#     # Run Medaka (replace 'r941_min_sup_g507' with correct model for your flowcell + basecaller)
#     # 运行Medaka（将'r941_min_sup_g507'替换为适合您的flowcell+basecaller的正确模型）
#     #Install Medaka
#     #安装Medaka
#     conda create -n medaka -c conda-forge -c nanoporetech -c bioconda medaka
#     #If already installed 
#     #如果已经安装
#     conda activate medaka
#     medaka_consensus -i temp/all.fq \
#     -d temp/otus_racon.fa \
#     -o temp/medaka_polished \
#     -t 14 \
#     -m r941_min_sup_g507
#     ##Then convert the consensus.fasta" into "otus.fa
#     ##然后将consensus.fasta转换为otus.fa
#     cd /d/EasyAmplicon_paper_materials/temp/medaka_polished
#     mv consensus.fasta temp/otus.fa
#     ##If size number present in consesus.fasta, remove using this code
#     ##如果consesus.fasta中存在size编号，请使用此代码将其删除
#     sed 's/;size=[0-9]*//' temp/otus.fa > temp/otus.fa
#     conda deactivate

### 4.4 基于参考去嵌合
### 4.4 Reference-based chimera detection

    # 不推荐，容易引起假阴性，因为参考数据库无丰度信息
    # Not recommended, as it can easily cause false negatives because the reference database lacks abundance information.
    # 而de novo时要求亲本丰度为嵌合体16倍以上防止假阴性
    # In de novo chimera detection, the abundance of the parent sequences is required to be more than 16 times that of the chimera to prevent false negatives.
    # 因为已知序列不会被去除，数据库选择越大越合理，假阴性率最低
    # Since known sequences will not be removed, the larger the database selected, the more reasonable it is, and the lower the false negative rate.
    mkdir -p result/raw

    # 方法1. vsearch+silva去嵌合
    # Method 1. Chimera removal with vsearch + SILVA
    vsearch --uchime_ref temp/otus.fa \
      -db ${db}/usearch/SILVA_modified.fasta \
      --nonchimeras result/raw/otus.fa \
      --threads 4
    # 4m41s, 325 chimeras
    # Win vsearch结果添加了windows换行符^M需删除，mac不要执行此命令
    # The results of vsearch on Windows have added a Windows newline character (^M), which needs to be deleted. Do not execute this command on a Mac.
    sed -i 's/\r//g' result/raw/otus.fa
    # 方法2. 不去嵌合
    # Method 2. Do not remove chimeras
    # cp -f temp/otus.fa result/raw/otus.fa


## 5. 特征表构建和筛选
## 5. Feature table creation and filtering

    # OTU和ASV统称为特征(Feature)，它们的区别是：
    # OTUs and ASVs are collectively referred to as Features. Their differences are:
    # OTU通常按97%聚类后挑选最高丰度或中心的代表性序列；
    # OTUs are usually representative sequences selected with the highest abundance or from the center after 97% clustering.
    # ASV是基于序列进行去噪(排除或校正错误序列，并挑选丰度较高的可信序列)作为代表性序列
    # ASVs are representative sequences obtained by denoising based on sequences (excluding or correcting erroneous sequences and selecting credible sequences with higher abundance).

### 5.1 生成特征表
### 5.1 Generate a Feature table
    # If you met this error "Fatal error: FASTA file expected, FASTQ file found (temp/filtered.fa)", than run this code 
    # 如果您遇到此错误“Fatal error: FASTA file expected, FASTQ file found (temp/filtered.fa)”，则运行此代码
    # vsearch --fastq_filter temp/filtered.fa \
    # --fastaout temp/filtered.fasta \
    # --fastq_qmax 50
    # 方法1. usearch生成特征表，小样本(<30)快；但大样本受限且多线程效率低，4971/21783(22.8%), 4核1m12s
    # Method 1. usearch generates a feature table, which is fast for small samples (<30); but it is limited for large samples and has low multi-threading efficiency. 4971/21783 (22.8%), 4 cores, 1m 12s.
    # time usearch -otutab temp/filtered.fa \
    #    -otus result/raw/otus.fa \
    #    -threads 4 \
    #    -otutabout result/raw/otutab.txt

    # 方法2. vsearch生成特征表
    # Method 2. vsearch generates a feature table
    # id(1)：100%相似度比对3797 of 21783 (17.43%) real time 3m14s
    # id(1): 100% similarity alignment, 3797 of 21783 (17.43%), real time 3m 14s.
    # 这里使用id(0.97)：(更高数据使用率，更快)Matching unique query sequences: 5147 of 21783 (23.63%),3m6s
    # Here, id(0.97) is used: (higher data usage, faster) Matching unique query sequences: 5147 of 21783 (23.63%), 3m 6s.
    time vsearch --usearch_global temp/filtered.fa \
      --db result/raw/otus.fa \
      --id 0.97 --threads 16 \
    	--otutabout result/raw/otutab.txt

    # vsearch结果windows用户删除换行符^M校正为标准Linux格式
    # For vsearch results, Windows users should delete the newline character ^M to correct it to the standard Linux format.
    sed -i 's/\r//' result/raw/otutab.txt
    head -n6 result/raw/otutab.txt | cut -f 1-6 |cat -A
    # csvtk统计表行列
    # Use csvtk to count the rows and columns of the table.
    # 这里一定看好列数，是不是等于你的样品数；如果不等，一般是样品命名存在问题，具体看上面解释
    # Be sure to check the number of columns here. Is it equal to your number of samples? If not, there is generally a problem with the sample naming. See the explanation above for details.
    csvtk -t stat result/raw/otutab.txt
    
### 5.2 物种注释，且/或去除质体和非细菌
### 5.2 Taxonomic annotation, and/or removal of plastids and non-Bacteria

    # 物种注释-去除质体和非细菌/古菌并统计比例(可选)
    # Taxonomic annotation - remove plastids and non-bacteria/archaea and calculate the proportion (optional).
    # SILVA数据库(SILVA_138.2_SSURef_NR99_tax_silva.fasta)更好注释真核、质体序列
    # The SILVA database (SILVA_138.2_SSURef_NR99_tax_silva.fasta) is better for annotating eukaryotic and plastid sequences.
    # 修改SILVA数据库(SILVA_138.2_SSURef_NR99_tax_silva.fasta)格式以适应注释代码，输出数据库文件SILVA_modified.fasta，包含种和株系的信息
    # Modify the format of the SILVA database (SILVA_138.2_SSURef_NR99_tax_silva.fasta) to adapt to the annotation code. The output database file is SILVA_modified.fasta, which contains species and strain information.
    # python database_silva.py
    # 置信阈值通常0.6/0.8，vserch最低0.1/usearch可选0输出最相似物种注释用于观察潜在分类,7m5s
    # The confidence threshold is usually 0.6/0.8. The minimum for vsearch is 0.1. usearch can optionally be set to 0 to output the most similar species annotation for observing potential classifications. Takes 7m 5s.
    vsearch --sintax result/raw/otus.fa \
      --db ${db}/usearch/SILVA_modified.fasta \
      --sintax_cutoff 0.1 \
      --tabbedout result/raw/otus.sintax 
    head result/raw/otus.sintax | cat -A
    sed -i 's/\r//' result/raw/otus.sintax


    # 方法1. 原始特征表行数
    # Method 1. Number of rows in the original feature table
    wc -l result/raw/otutab.txt
    #R脚本选择细菌古菌(真核)、去除叶绿体、线粒体并统计比例；输出筛选并排序的OTU表
    # R script to select bacteria/archaea (eukaryotes), remove chloroplasts, mitochondria, and calculate proportions; output a filtered and sorted OTU table.
    #输入为OTU表result/raw/otutab.txt和物种注释result/raw/otus.sintax
    # The input is the OTU table result/raw/otutab.txt and the taxonomic annotation result/raw/otus.sintax.
    #输出筛选并排序的特征表result/otutab.txt和
    # The output is the filtered and sorted feature table result/otutab.txt and
    #统计污染比例文件result/raw/otutab_nonBac.txt和过滤细节otus.sintax.discard
    # the contamination proportion file result/raw/otutab_nonBac.txt and the filtering details otus.sintax.discard.
    #真菌ITS数据，请改用otutab_filter_nonFungi.R脚本，只筛选真菌
    # For fungal ITS data, please use the otutab_filter_nonFungi.R script to filter only fungi.
    # Rscript ${db}/script/otutab_filter_nonBac.R -h # 显示参数说明 / Display parameter description
    Rscript ${db}/script/otutab_filter_nonBac.R \
      --input result/raw/otutab.txt \
      --taxonomy result/raw/otus.sintax \
      --output result/otutab.txt\
      --stat result/raw/otutab_nonBac.stat \
      --discard result/raw/otus.sintax.discard
    # 筛选后特征表行数
    # Number of rows in the feature table after filtering
    wc -l result/otutab.txt
    #过滤特征表对应序列
    # Filter the corresponding sequences in the feature table
    cut -f 1 result/otutab.txt | tail -n+2 > result/otutab.id
    usearch -fastx_getseqs result/raw/otus.fa \
        -labels result/otutab.id -fastaout result/otus.fa
    #过滤特征表对应序列注释
    # Filter the corresponding sequence annotations in the feature table
    awk 'NR==FNR{a[$1]=$0}NR>FNR{print a[$1]}'\
        result/raw/otus.sintax result/otutab.id \
        > result/otus.sintax

    # 方法2. 觉得筛选不合理可以不筛选，本次采用不筛选的方法
    # Method 2. If you think the filtering is unreasonable, you can skip it. This time, we will use the method without filtering.
     cp result/raw/otu* result/

    #可选统计方法：OTU表简单统计 Summary OTUs table
    # Optional statistical method: Simple summary of the OTU table
    usearch -otutab_stats result/otutab.txt \
      -output result/otutab.stat
    cat result/otutab.stat
    #注意最小值、分位数，或查看result/raw/otutab_nonBac.stat中样本详细数据量，用于重采样
    # Pay attention to the minimum value, quantiles, or view the detailed sample data volume in result/raw/otutab_nonBac.stat for resampling.

### 5.3 等量抽样标准化
### 5.3 Normalization by subsampling

    # Normlize by subsample
    # 通过子抽样进行标准化

    #使用vegan包进行等量重抽样，输入reads count格式Feature表result/otutab.txt
    # Use the vegan package for equal resampling. The input is the feature table in reads count format: result/otutab.txt.
    #可指定输入文件、抽样量和随机数，输出抽平表result/otutab_rare.txt和多样性alpha/vegan.txt
    # You can specify the input file, sampling depth, and random seed. The output is the rarefied table result/otutab_rare.txt and the diversity file alpha/vegan.txt.
    mkdir -p result/alpha
    Rscript ${db}/script/otutab_rare.R --input result/otutab.txt \
      --depth 139 --seed 1 \
      --normalize result/otutab_rare.txt \
      --output result/alpha/vegan.txt
    usearch -otutab_stats result/otutab_rare.txt \
      -output result/otutab_rare.stat
    cat result/otutab_rare.stat


## 6. α多样性
## 6. Alpha diversity

### 6.1. 计算α多样性
### 6.1. Calculate alpha diversity

    # 使用USEARCH计算14种alpha多样性指数(Chao1有错勿用)
    # Use USEARCH to calculate 14 alpha diversity indices (Chao1 has errors, do not use).
    #details in http://www.drive5.com/usearch/manual/alpha_metrics.html
    usearch -alpha_div result/otutab_rare.txt \
      -output result/alpha/alpha.txt

### 6.2. 计算稀释丰富度
### 6.2. Calculate rarefaction richness

    #稀释曲线：取1%-100%的序列中OTUs数量，每次无放回抽样
    # Rarefaction curve: Take the number of OTUs in 1%-100% of the sequences, sampling without replacement each time.
    #Rarefaction from 1%, 2% .. 100% in richness (observed OTUs)-method without_replacement https://drive5.com/usearch/manual/cmd_otutab_subsample.html
    usearch -alpha_div_rare result/otutab_rare.txt \
      -output result/alpha/alpha_rare.txt \
      -method without_replacement
    #预览结果
    # Preview the results
    head -n2 result/alpha/alpha_rare.txt
    #样本测序量低出现非数值"-"的处理，详见常见问题8
    # For samples with low sequencing depth, non-numeric values "-" may appear. See FAQ 8 for details on how to handle this.
    sed -i "s/-/\t0.0/g" result/alpha/alpha_rare.txt

### 6.3. 筛选高丰度菌
### 6.3. Filter by abundance

    #计算各特征的均值，有组再求分组均值，需根据实验设计metadata.txt修改组列名
    # Calculate the mean of each feature. If there are groups, calculate the group means. You need to modify the group column name according to the experimental design metadata.txt.
    #输入文件为feautre表result/otutab.txt，实验设计metadata.txt
    # The input files are the feature table result/otutab.txt and the experimental design metadata.txt.
    #输出为特征表按组的均值-一个实验可能有多种分组方式
    # The output is the mean of the feature table by group. An experiment may have multiple grouping methods.
    #-h显示脚本帮助(参数说明)
    # -h displays the script help (parameter description).
    Rscript ${db}/script/otu_mean.R -h
    #scale是否标准化，zoom标准化总和，all输出全部样本均值，type计算类型mean或sum
    # scale: whether to standardize; zoom: standardize the sum; all: output the mean of all samples; type: calculation type, mean or sum.
    Rscript ${db}/script/otu_mean.R --input result/otutab.txt \
      --metadata result/metadata.txt \
      --group Group --thre 0 \
      --scale TRUE --zoom 100 --all TRUE --type mean \
      --output result/otutab_mean.txt
    # 结果为全部和各组均值
    # The result is the mean of all samples and each group.
    head -n3 result/otutab_mean.txt

    #如以平均丰度>0.1%筛选，可选0.5或0.05，得到每个组的OTU组合
    # For example, filter by an average abundance of >0.1% (you can choose 0.5 or 0.05) to get the OTU combination for each group.
    awk 'BEGIN{OFS=FS="\t"}{if(FNR==1) {for(i=3;i<=NF;i++) a[i]=$i; print "OTU","Group";} \
        else {for(i=3;i<=NF;i++) if($i>0.1) print $1, a[i];}}' \
        result/otutab_mean.txt > result/alpha/otu_group_exist.txt
    head result/alpha/otu_group_exist.txt
    cut -f 2 result/alpha/otu_group_exist.txt | sort | uniq -c
    # 试一试：不同丰度下各组有多少OTU/ASV
    # Try it: How many OTUs/ASVs are there in each group at different abundances?
    # 可在 http://ehbio.com/test/venn/ 中绘图并显示各组共有和特有维恩或网络图
    # You can draw and display the common and unique Venn or network diagrams of each group at http://ehbio.com/test/venn/.
    # 也可在 http://www.ehbio.com/ImageGP 绘制Venn、upSetView和Sanky
    # You can also draw Venn, upSetView, and Sankey diagrams at http://www.ehbio.com/ImageGP.

## 7. β多样性
## 7. Beta diversity

    #结果有多个文件，需要目录
    # The results have multiple files, so a directory is needed.
    mkdir -p result/beta/
    #基于OTU构建进化树 Make OTU tree, 3min40s
    # Build a phylogenetic tree based on OTUs. Takes 3m 40s.
    usearch -cluster_agg result/otus.fa -treeout result/otus.tree
    #1s生成5种距离矩阵：bray_curtis, euclidean, jaccard, manhatten, unifrac
    # Generate 5 distance matrices in 1s: bray_curtis, euclidean, jaccard, manhattan, unifrac.
    usearch -beta_div result/otutab_rare.txt -tree result/otus.tree \
    -filename_prefix result/beta/


## 8. 物种注释分类汇总
## 8. Summary of taxonomic annotation

    # OTU对应物种注释2列格式：去除sintax中置信值，只保留物种注释，替换:为_，删除引号
    # 2-column format for OTU corresponding species annotation: remove the confidence value in sintax, keep only the species annotation, replace ":" with "_", and delete quotation marks.
    cut -f 1,4 result/otus.sintax \
      |sed 's/\t/\tk/;s/:/__/g;s/,/;/g;s/"//g' \
      > result/taxonomy2.txt
    head -n3 result/taxonomy2.txt

    #OTU对应物种8列格式：注意注释是非整齐
    # 8-column format for OTU corresponding species: Note that the annotation is not uniform.
    #生成物种表格OTU/ASV中空白补齐为Unassigned
    # In the generated species table, fill in the blanks in OTU/ASV with "Unassigned".
    awk 'BEGIN{OFS=FS="\t"}{delete a; a["k"]="Unassigned";a["p"]="Unassigned";a["c"]="Unassigned";a["o"]="Unassigned";a["f"]="Unassigned";a["g"]="Unassigned";a["s"]="Unassigned";\
      split($2,x,";");for(i in x){split(x[i],b,"__");a[b[1]]=b[2];} \
      print $1,a["k"],a["p"],a["c"],a["o"],a["f"],a["g"],a["s"];}' \
      result/taxonomy2.txt > temp/otus.tax
      
    sed 's/;/\t/g;s/.__//g;' temp/otus.tax|cut -f 1-8 | \
      sed '1 s/^/OTUID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n/' \
      > result/taxonomy.txt
    head -n3 result/taxonomy.txt

    # 补齐otu.sintax中未注释的行
    # Fill in the unannotated lines in otu.sintax
    awk -F"\t" 'BEGIN{OFS="\t"} {
    if ($2 == "") $2 = "d:(0.00)";
    if ($3 == "") $3 = "+";
    if ($4 == "") $4 = "d:Unassigned";
    print
    }' result/otus.sintax > result/otus.sintax.filled

    #统计门纲目科属，使用 rank参数 p c o f g s，为phylum, class, order, family, genus, species缩写
    # Count phylum, class, order, family, genus, species. Use the rank parameter p, c, o, f, g, s, which are abbreviations for phylum, class, order, family, genus, species.
    mkdir -p result/tax
    for i in p c o f g s;do
      usearch -sintax_summary result/otus.sintax.filled \
      -otutabin result/otutab_rare.txt -rank ${i} \
      -output result/tax/sum_${i}.txt
    done
    sed -i 's/(//g;s/)//g;s/"//g;s/#//g;s/\/Chloroplast//g' result/tax/sum_*.txt
    # 列出所有文件
    # List all files
    wc -l result/tax/sum_*.txt
    head -n3 result/tax/sum_g.txt

## 9. 有参定量特征表
## 9. Reference-based quantitative feature table

    # 比对Greengenes97% OTUs比对，用于PICRUSt/Bugbase功能预测
    # Align to Greengenes 97% OTUs for PICRUSt/Bugbase functional prediction.
    mkdir -p result/gg/

    #方法1. usearch比对更快，但文件超限报错选方法2
    # Method 1. usearch alignment is faster, but if the file size exceeds the limit and an error is reported, choose method 2.
    # 默认10核以下使用1核，10核以上使用10核
    # By default, 1 core is used for less than 10 cores, and 10 cores are used for more than 10 cores.
    usearch -otutab temp/filtered.fa -otus ${db}/gg/97_otus.fa \
    	-otutabout result/gg/otutab.txt -threads 10
    # 3707 / 21783 mapped to OTUs (17.0%)，10核41s, 内存使用731Mb
    # 3707 / 21783 mapped to OTUs (17.0%), 10 cores, 41s, memory usage 731Mb.
    head -n3 result/gg/otutab.txt

    # #方法2（本次使用）. vsearch比对，更准更慢，但并行24-96线程更强
    # # Method 2 (used this time). vsearch alignment is more accurate but slower, but it is more powerful with 24-96 parallel threads.
    vsearch --usearch_global temp/filtered.fa --db ${db}/gg/97_otus.fa \
       --otutabout result/gg/otutab.txt --id 0.97 --threads 12
    # 4063 of 21783 (18.65%)，12核4m24s
    # 4063 of 21783 (18.65%), 12 cores, 4m 24s.

    #统计
    # Statistics
    usearch -otutab_stats result/gg/otutab.txt -output result/gg/otutab.stat
    cat result/gg/otutab.stat


## 10. 空间清理及数据提交
## 10. Workspace cleanup and data submission

    #删除中间大文件
    # Delete large intermediate files
    rm -rf temp/*.fq

    # 统计md5值，用于数据提交
    # Calculate md5 values for data submission
    cd seq
    md5sum *.fastq > ../result/md5sum.txt
    cat ../result/md5sum.txt

# 除了上述传统的OTU/ASV分析路径，本脚本还提供了另一种基于Emu软件的快速物种注释和定量流程
# In addition to the traditional OTU/ASV analysis path described above, this script also provides another rapid species annotation and quantification process based on the Emu software.
# Emu通过直接将序列比对到参考数据库来进行物种鉴定和丰度计算，尤其适合需要快速获得物种组成
# Emu performs species identification and abundance calculation by directly aligning sequences to a reference database, which is particularly suitable for quickly obtaining species composition.
# 基于Emu软件的物种注释流程如下：
# The species annotation process based on the Emu software is as follows:

    ### 步骤1：为每个样本的Reads重新打上标签
    ### Step 1: Relabel Each Sample’s Reads
    # 为了便于后续追溯和合并，Emu流程的第一步是为每个样本的序列（Reads）添加一个唯一的样本标识符作为前缀。
    # To facilitate subsequent tracking and merging, the first step in the Emu workflow is to add a unique sample identifier as a prefix to the sequences (Reads) of each sample.
    # 这样可以确保在多样本混合分析时，每条序列都能准确地归属到其原始样本。
    # This ensures that each sequence can be accurately attributed to its original sample during multi-sample analysis.

    # 切换到Nanopore分析的主目录
    # Change to the main directory for Nanopore analysis
    cd /c/EasyAmplicon2/Nanopore/
    # 创建一个临时目录temp2来存放重命名后的文件
    # Create a temporary directory temp2 to store the relabeled files
    mkdir -p temp2

    # 使用一个for循环遍历seq/目录下的所有.fastq文件
    # Use a for loop to iterate through all .fastq files in the seq/ directory
    for fastq in seq/*.fastq; do
        # 从完整路径中提取不带扩展名的样本名
        # Extract the sample name without the extension from the full path
        sample=$(basename "$fastq" .fastq)
        # 使用usearch的-fastx_relabel命令为序列添加前缀
        # Use the -fastx_relabel command of usearch to add a prefix to the sequences
        # -prefix "${sample}." 会将如"read1"这样的序列ID变为"E1.read1"
        # -prefix "${sample}." will change a sequence ID like "read1" to "E1.read1"
        usearch -fastx_relabel "$fastq" -fastqout "temp2/${sample}.fq" -prefix "${sample}."
    done
    # 这个循环处理完后，temp2/目录下会包含所有样本的、且序列ID已更新的FASTQ文件。
    # After this loop is processed, the temp2/ directory will contain FASTQ files for all samples with updated sequence IDs.

    ### 步骤2：对每个样本进行引物修剪和质量过滤
    ### Step 2: Trim Primers and Quality Filter Each Sample
    # 这个步骤旨在去除测序读两端的引物序列，并根据序列长度和质量值进行过滤，以提高后续分析的准确性。
    # This step aims to remove primer sequences from both ends of the sequencing reads and filter them based on sequence length and quality values to improve the accuracy of subsequent analysis.
    
    # 切换到Linux Bash环境并激活conda环境
    # Switch to the Linux Bash environment and activate the conda environment
    bash
    conda activate easyamplicon2
    cd /mnt/c/EasyAmplicon2/Nanopore/

    # 使用cutadapt进行引物修剪
    # Primer trimming with cutadapt
    # 遍历temp2/中所有重命名后的.fq文件
    # Loop through all relabeled .fq files in temp2/
    for fq in temp2/*.fq; do
        sample=$(basename "$fq" .fq)
        # -g 指定5'端引物序列，...AAGTC... 指定3'端引物序列
        # -g specifies the 5' primer sequence, ...AAGTC... specifies the 3' primer sequence
        # --action=trim 保留匹配引物内部的序列
        # --action=trim keeps the sequence inside the matching primers
        # --rc 检查反向互补链的引物
        # --rc checks for primers on the reverse complement strand
        # -m 和 -M 分别设置序列的最小和最大长度
        # -m and -M set the minimum and maximum length of the sequences, respectively
        # --discard-untrimmed 丢弃未找到引物的序列
        # --discard-untrimmed discards sequences where primers were not found
        cutadapt -g "AGAGTTTGATCCTGGCTCAG...AAGTCSTAACAAGGTADCCSTA" \
            --error-rate=0.1 \
            --action=trim \
            --rc \
            -m 1000 \
            -M 1800 \
            -j 8 \
            --discard-untrimmed \
            -o "temp2/${sample}.filtered.fastq" \
            "$fq"
    done

    # 使用NanoFilt进行质量过滤
    # Quality filtering with NanoFilt
    # 遍历temp2/中所有经过引物修剪的.fastq文件
    # Loop through all primer-trimmed .fastq files in temp2/
    for fastq in temp2/*.filtered.fastq; do
        sample=$(basename "$fastq" .filtered.fastq)
        # -l 设置最小长度, -q 设置最低质量分数 (Q-score)
        # -l sets the minimum length, -q sets the minimum quality score (Q-score)
        # --headcrop 去除序列头部的碱基, --maxlength 设置最大长度
        # --headcrop removes bases from the head of the sequence, --maxlength sets the maximum length
        cat "$fastq" | NanoFilt -l 1000 -q 18 --headcrop 10 --maxlength 1800 > "temp2/${sample}.filtered.qc.fastq"
    done

    ### 步骤3：对每个样本运行Emu进行物种注释和丰度计算
    ### Step 3: Run Emu on Each Sample for Taxonomic Classification and Abundance Estimation
    # Emu通过将质控后的序列比对到参考数据库（如SILVA或GTDB），来快速确定物种组成和它们的相对丰度。
    # Emu rapidly determines species composition and their relative abundances by aligning quality-controlled sequences to a reference database (such as SILVA or GTDB).

    # 针对SILVA数据库运行Emu
    # Run Emu against the SILVA database
    for fastq in temp2/*.filtered.qc.fastq; do
        sample=$(basename "$fastq" .filtered.qc.fastq)
        # 为每个样本的输出创建一个单独的目录
        # Create a separate directory for the output of each sample
        mkdir -p "result/Emu/silva/$sample"
        # emu abundance 是核心命令
        # emu abundance is the core command
        # --type map-ont 指定使用为Nanopore数据优化的比对模式
        # --type map-ont specifies the alignment mode optimized for Nanopore data
        # --db 指向预先格式化好的Emu参考数据库路径
        # --db points to the path of the pre-formatted Emu reference database
        # --output-dir 指定结果输出目录
        # --output-dir specifies the output directory for the results
        emu abundance "$fastq" \
            --type map-ont \
            --db /mnt/c/EasyMicrobiome/Silva_Emu \
            --output-dir "result/Emu/silva/$sample" \
            --threads 4
    done

    # 针对GTDB数据库运行Emu
    # Run Emu against the GTDB database
    for fastq in temp2/*.filtered.qc.fastq; do
        sample=$(basename "$fastq" .filtered.qc.fastq)
        mkdir -p "result/Emu/gtdb/$sample"
        emu abundance "$fastq" \
            --type map-ont \
            --db /mnt/c/EasyMicrobiome/GTDB_Emu \
            --output-dir "result/Emu/gtdb/$sample" \
            --threads 8
    done

    ### 步骤4：合并与整理Emu输出结果
    ### Step 4: Combine and Tidy Emu Outputs
    # 单独运行Emu会为每个样本生成一个结果文件。为了进行跨样本比较，需要将这些结果合并成一个统一的物种丰度表。
    # Running Emu separately generates a result file for each sample. To perform cross-sample comparisons, these results need to be combined into a unified species abundance table.
    # 在合并之前，可能需要对Emu输出的TSV文件进行格式清理，以确保分类学列名和层级结构的一致性。
    # Before merging, it may be necessary to clean up the format of the Emu output TSV files to ensure consistency in taxonomic column names and hierarchical structure.

    # 首先，检查Emu输出的分类学列名是否规范
    # First, check if the taxonomic column names in the Emu output are standardized
    head -n 1 result/Emu/silva/E1/E1.filtered.qc_rel-abundance.tsv
    
    # 如果列名不规范，运行Python脚本进行重命名和清理
    # If the column names are not standard, run Python scripts for renaming and cleaning
    # 对SILVA结果进行清理
    # Clean up SILVA results
    python /mnt/c/EasyMicrobiome/script/renameSILVA_eachsample_columns.py
    # 对GTDB结果进行清理，并移除分类名前缀（如'g__'）
    # Clean up GTDB results and remove taxonomic prefixes (e.g., 'g__')
    python /mnt/c/EasyMicrobiome/script/renameGTDB_eachsample_columns.py
    python /mnt/c/EasyMicrobiome/script/Remove_gtdb_prefix.py

    ### 步骤5：使用Emu合并多样本结果
    ### Step 5: Combine Multi-Sample Results Using Emu
    # 清理格式后，使用`emu combine-outputs`命令在不同的分类学级别（物种、属、科等）上合并所有样本的丰度数据。
    # After cleaning the format, use the `emu combine-outputs` command to merge the abundance data of all samples at different taxonomic levels (species, genus, family, etc.).

    # 合并GTDB结果
    # Combine GTDB results
    # emu combine-outputs [输入目录] [分类学级别]
    # emu combine-outputs [input_directory] [taxonomic_rank]
    emu combine-outputs result/Emu/gtdb_cleaned_noprefix/ species
    emu combine-outputs result/Emu/gtdb_cleaned_noprefix/ genus
    emu combine-outputs result/Emu/gtdb_cleaned_noprefix/ family
    emu combine-outputs result/Emu/gtdb_cleaned_noprefix/ class
    emu combine-outputs result/Emu/gtdb_cleaned_noprefix/ order
    emu combine-outputs result/Emu/gtdb_cleaned_noprefix/ phylum

    # 合并SILVA结果
    # Combine SILVA results
    emu combine-outputs result/Emu/silva_cleaned/ species
    emu combine-outputs result/Emu/silva_cleaned/ genus
    emu combine-outputs result/Emu/silva_cleaned/ family
    emu combine-outputs result/Emu/silva_cleaned/ order
    emu combine-outputs result/Emu/silva_cleaned/ class
    emu combine-outputs result/Emu/silva_cleaned/ phylum



# PICRUSt2功能预测
# PICRUSt2 Functional Prediction

    # # 由于Nanopore数据质量较差，故不推荐做功能预测
    # # (可选)PICRUSt2(Linux/Windows下Linux子系统，要求>16GB内存)
    # # (Optional) PICRUSt2 (Linux/Linux subsystem on Windows, requires >16GB memory)
    # # 方法1. 直接安装
    # # Method 1. Direct installation
    # n=picrust2
    # conda create -n ${n} -c bioconda -c conda-forge ${n}=2.3.0_b
    # # 方法2. 导入安装环境(推荐)
    # # Method 2. Import installation environment (recommended)
    # n=picrust2
    # # 复制安装包，或下载我的环境打包
    # # Copy the installation package, or download my packaged environment
    # wget -c ftp://download.nmdc.cn/tools/conda/${n}.tar.gz
    # # 指定安装目录并解压
    # # Specify the installation directory and decompress
    # condapath=~/miniconda3
    # mkdir -p ${condapath}/envs/${n}
    # tar -xvzf ${n}.tar.gz -C ${condapath}/envs/${n}
    # # 激活环境并初始化
    # # Activate the environment and initialize
    # conda activate picrust2
    # conda unpack
    # # 加载环境
    # # Load the environment
    # conda activate picrust2
    # # 进入工作目录，服务器要修改工作目录
    # # Enter the working directory. The working directory needs to be modified on the server.
    # mkdir -p /mnt/c/EasyAmplicon2/Nanopore/result/picrust2
    # cd /mnt/c/EasyAmplicon2/Nanopore/result/picrust2
    # # 运行流程，内存15.7GB，耗时12m；内存不足会导致程序中断，out目录里只有一个intermediate/文件夹，需使用大内存设备运行
    # # Run the pipeline. Memory: 15.7GB, time: 12m. Insufficient memory will cause the program to be interrupted. The out directory will only have an intermediate/ folder. A device with large memory is required.
    # picrust2_pipeline.py -s ../otus.fa -i ../otutab.txt -o ./out -p 8
    # # 添加EC/KO/Pathway注释
    # # Add EC/KO/Pathway annotations
    # cd out
    # add_descriptions.py -i pathways_out/path_abun_unstrat.tsv.gz -m METACYC \
    # -o pathways_out/path_abun_unstrat_descrip.tsv.gz
    # add_descriptions.py -i EC_metagenome_out/pred_metagenome_unstrat.tsv.gz -m EC \
    #   -o EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz
    # add_descriptions.py -i KO_metagenome_out/pred_metagenome_unstrat.tsv.gz -m KO \
    #   -o KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz 
    # # KEGG按层级合并
    # # Merge KEGG by hierarchy
    # db=/mnt/c/EasyMicrobiome/
    # zcat KO_metagenome_out/pred_metagenome_unstrat.tsv.gz > KEGG.KO.txt
    # python3 ${db}/script/summarizeAbundance.py \
    #   -i KEGG.KO.txt \
    #     -m ${db}/kegg/KO1-4.txt \
    #     -c 2,3,4 -s ',+,+,' -n raw \
    #     -o KEGG
    # # 统计各层级特征数量
    # # Count the number of features at each level
    # wc -l KEGG*


# Evolution进化树
# Evolution Tree

    #切换gitbash
    # Switch to gitbash
    cd ${wd}/Nanopore
    mkdir -p result/tree
    cd ${wd}/Nanopore/result/tree

## 筛选高丰度/指定的特征
## Filter high-abundance/specified features

    #方法1. 按丰度筛选特征，一般选0.001或0.005，且OTU数量在30-150个范围内
    # Method 1. Filter features by abundance, usually 0.001 or 0.005, and the number of OTUs is in the range of 30-150.
    #统计特征表中ASV数量，如总计3794个
    # Count the number of ASVs in the feature table, e.g., a total of 3794.
    tail -n+2 ../otutab_rare.txt | wc -l
    #按相对丰度0.2%筛选高丰度OTU
    # Filter high-abundance OTUs by a relative abundance of 0.2%.
    usearch -otutab_trim ../otutab_rare.txt \
        -min_otu_freq 0.002 \
        -output otutab.txt
    #统计筛选OTU表特征数量，总计~19个
    # Count the number of features in the filtered OTU table, a total of ~19.
    tail -n+2 otutab.txt | wc -l

    #方法2. 按数量筛选
    # Method 2. Filter by quantity
    # #按丰度排序，默认由大到小
    # # Sort by abundance, from largest to smallest by default.
    # usearch -otutab_sortotus ../otutab_rare.txt  \
    #     -output otutab_sort.txt
    # #提取高丰度中指定Top数量的OTU ID，如Top100,
    # # Extract the specified number of top OTU IDs from the high-abundance OTUs, e.g., Top100.
    # sed '1 s/#OTU ID/OTUID/' otutab_sort.txt \
    #     | head -n101 > otutab.txt

    #修改特征ID列名
    # Modify the feature ID column name
    sed -i '1 s/#OTU ID/OTUID/' otutab.txt
    #提取ID用于提取序列
    # Extract IDs for sequence extraction
    cut -f 1 otutab.txt > otutab_high.id

    # 筛选高丰度菌/指定差异菌对应OTU序列
    # Filter the OTU sequences corresponding to high-abundance bacteria/specified differential bacteria.
    usearch -fastx_getseqs ../otus.fa -labels otutab_high.id \
        -fastaout otus.fa
    head -n 2 otus.fa

    ## 筛选OTU对物种注释
    ## Filter OTU for species annotation
    awk 'NR==FNR{a[$1]=$0} NR>FNR{print a[$1]}' ../taxonomy.txt \
        otutab_high.id > otutab_high.tax

    #获得OTU对应组均值，用于样本热图
    # Get the mean of the OTU corresponding group for the sample heatmap.
    #依赖之前otu_mean.R计算过按Group分组的均值
    # Depends on the mean calculated by group in the previous otu_mean.R.
    awk 'NR==FNR{a[$1]=$0} NR>FNR{print a[$1]}' ../otutab_mean.txt otutab_high.id \
        | sed 's/#OTU ID/OTUID/' > otutab_high.mean
    head -n3 otutab_high.mean

    #合并物种注释和丰度为注释文件
    # Merge species annotation and abundance into an annotation file.
    cut -f 2- otutab_high.mean > temp
    paste otutab_high.tax temp > annotation.txt
    head -n 3 annotation.txt

## 构建进化树
## Build a phylogenetic tree

    # 起始文件为 result/tree目录中 otus.fa(序列)、annotation.txt(物种和相对丰度)文件
    # The starting files are otus.fa (sequences) and annotation.txt (species and relative abundance) in the result/tree directory.
    # Muscle软件进行序列对齐，3s
    # Use Muscle for sequence alignment, takes 3s.
    muscle -in otus.fa -out otus_aligned.fas

    ### 方法1. 利用IQ-TREE快速构建ML进化树，2m
    ### Method 1. Use IQ-TREE to quickly build an ML phylogenetic tree, takes 2m.
    # rm -rf iqtree
    # mkdir -p iqtree
    # iqtree -s otus_aligned.fas \
    #  -bb 1000 -redo -alrt 1000 -nt AUTO \
    #  -pre iqtree/otus

    ### 方法2. FastTree快速建树(Linux)
    ### Method 2. FastTree for fast tree building (Linux)
    # 注意FastTree软件输入文件为fasta格式的文件，而不是通常用的Phylip格式。输出文件是Newick格式。
    # Note that the input file for FastTree is in FASTA format, not the commonly used Phylip format. The output file is in Newick format.
    # 该方法适合于大数据，例如几百个OTUs的系统发育树！
    # This method is suitable for big data, such as phylogenetic trees with hundreds of OTUs!
    # Ubuntu上安装fasttree可以使用`sudo apt install fasttree`
    # On Ubuntu, you can use `sudo apt install fasttree` to install fasttree.
    cd result/tree
    fasttree -gtr -nt otus_aligned.fas > otus.nwk
    
